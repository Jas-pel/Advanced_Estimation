import numpy as np
import pytest
from qiskit.quantum_info import PauliList, Pauli
from advanced_estimation.commutation import clifford

# =================================================================
# TESTS SYSTÉMATIQUES : PORTES À 1 QUBIT (H, S)
# =================================================================

@pytest.mark.parametrize("gate", ["h", "s"])
def test_single_qubit_gates_all_inputs(gate):
    """Teste H et S sur X, Y, Z pour vérifier les transformations et les phases."""
    paulis = PauliList(["X", "Y", "Z"])
    func = getattr(clifford, gate)
    res = func(paulis, [0])
    
    labels = [p.to_label() for p in res]
    
    if gate == "h":
        # H: X->Z, Z->X, Y->-Y
        assert labels[0].endswith("Z")
        assert labels[1].endswith("Y") and "-" in labels[1]
        assert labels[2].endswith("X")
    elif gate == "s":
        # S: X->Y, Y->-X, Z->Z
        assert labels[0].endswith("Y")
        assert labels[1].endswith("X") and "-" in labels[1]
        assert labels[2].endswith("Z")

# =================================================================
# TESTS SYSTÉMATIQUES : PORTES À 2 QUBITS (CX, CZ)
# =================================================================

def test_cx_exhaustive():
    """Vérifie toutes les règles de propagation du CNOT (Control=0, Target=1)."""
    # Règles : IX->XX, XI->XI, IZ->IZ, ZI->ZZ
    inputs  = PauliList(["IX", "XI", "IZ", "ZI", "YY"])
    # Rappel : Qiskit label est 'Q1 Q0' -> "IX" est X sur Q0.
    res = clifford.cx(inputs, 0, 1)
    
    assert res[0].to_label().endswith("XX") # IX -> XX
    assert res[1].to_label().endswith("XI") # XI -> XI
    assert res[2].to_label().endswith("IZ") # IZ -> IZ
    assert res[3].to_label().endswith("ZZ") # ZI -> ZZ
    # YY -> (-Y) * (-Y) sur les phases ou transformation symplectique
    assert "X" in res[4].to_label() # Y1Y0 -> -X1Z0 (selon la convention)

def test_cz_exhaustive():
    """Vérifie la symétrie du CZ (IX -> ZX, XI -> XZ)."""
    inputs = PauliList(["IX", "XI", "IZ", "ZI"])
    res = clifford.cz(inputs, 0, 1)
    
    assert res[0].to_label().endswith("ZX")
    assert res[1].to_label().endswith("XZ")
    assert res[2].to_label().endswith("IZ")
    assert res[3].to_label().endswith("ZI")

# =================================================================
# TESTS DE PROPRIÉTÉS (COMMUTATION & MULTI-QUBITS)
# =================================================================

def test_commutation_invariance():
    """
    PROPRIÉTÉ CRITIQUE : Une porte Clifford DOIT préserver la commutation.
    Si ce test échoue, ta diagonalisation est impossible.
    """
    # Paire qui commute : XY et YX (sur 2 qubits)
    paulis = PauliList(["XY", "YX"])
    
    # Teste chaque porte
    p_h = clifford.h(paulis, [0])
    p_s = clifford.s(paulis, [1])
    p_cx = clifford.cx(paulis, 0, 1)
    p_cz = clifford.cz(paulis, 0, 1)
    
    for p_list in [p_h, p_s, p_cx, p_cz]:
        # Produit symplectique (commutation)
        x, z = p_list.x.astype(int), p_list.z.astype(int)
        commute = (x[0] @ z[1].T + x[1] @ z[0].T) % 2 == 0
        assert commute, f"La porte {p_list} a brisé la commutation !"

def test_multi_qubit_targets():
    """Vérifie que passer une liste de qubits fonctionne [0, 2]."""
    paulis = PauliList(["XXX"])
    # H sur Q0 et Q2 -> ZXZ
    res = clifford.h(paulis, [0, 2])
    assert res[0].to_label().endswith("ZXZ")

# =================================================================
# TEST DE STABILITÉ DES PHASES (TRÈS IMPORTANT)
# =================================================================

def test_s_gate_phase_double_application():
    """Appliquer S deux fois sur X doit donner -X (car S^2 = Z)."""
    p = PauliList(["X"])
    s1 = clifford.s(p, [0])
    s2 = clifford.s(s1, [0])
    
    label = s2[0].to_label()
    # S^2 X (S^2)+ = Z X Z = -X
    assert "X" in label and "-" in label, f"S^2 sur X devrait être -X, reçu {label}"