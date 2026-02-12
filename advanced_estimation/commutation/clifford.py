from typing import Union
import numpy as np
from qiskit.quantum_info import PauliList

def h(paulis: PauliList, qubits: Union[int, list[int]]) -> PauliList:
    if isinstance(qubits, int):
        qubits = [qubits]

    new_z = paulis.z.copy()
    new_x = paulis.x.copy()

    z_table = new_z[:, qubits].astype(int)
    x_table = new_x[:, qubits].astype(int)

    new_phase = paulis.phase + 2 * np.einsum('pq, pq -> p', z_table, x_table)
    new_phase = np.mod(new_phase, 4)

    new_z[:, qubits] = x_table
    new_x[:, qubits] = z_table

    return PauliList.from_symplectic(new_z, new_x, new_phase)


def s(paulis: PauliList, qubits: Union[int, list[int]]) -> PauliList:
    if isinstance(qubits, int):
        qubits = [qubits]

    new_z = paulis.z.copy()
    new_x = paulis.x.copy()

    z_table = new_z[:, qubits].astype(int)
    x_table = new_x[:, qubits].astype(int)

    new_phase = paulis.phase + 2 * np.einsum('pq, pq -> p', z_table, x_table)
    new_phase = np.mod(new_phase, 4)

    new_z[:, qubits] = np.logical_xor(z_table, x_table)

    return PauliList.from_symplectic(new_z, new_x, new_phase)


def cx(paulis: PauliList, control_qubits: Union[int, list[int]], target_qubits: Union[int, list[int]]) -> PauliList:
    if isinstance(control_qubits, int):
        control_qubits = [control_qubits]
    if isinstance(target_qubits, int):
        target_qubits = [target_qubits]
    assert len(control_qubits) == len(target_qubits)

    new_z = paulis.z.copy()
    new_x = paulis.x.copy()
    new_phase = paulis.phase.copy()

    for c, t in zip(control_qubits, target_qubits):
        zc = new_z[:, c]
        xc = new_x[:, c]
        zt = new_z[:, t]
        xt = new_x[:, t]

        term_phase = xc.astype(int) * zt.astype(int) * (1 - np.logical_xor(zc, xt).astype(int))
        add_phase = 2 * term_phase
        new_phase = np.mod(new_phase + add_phase, 4)

        new_z[:, c] = np.logical_xor(zc, zt)
        new_x[:, t] = np.logical_xor(xt, xc)

    return PauliList.from_symplectic(new_z, new_x, new_phase)


def cz(paulis: PauliList, control_qubits: Union[int, list[int]], target_qubits: Union[int, list[int]]) -> PauliList:
    # NOTE : Cette implémentation utilise une boucle plutôt qu'une approche vectorisée.
    # La raison : lorsque plusieurs paires control-target partagent des qubits 
    # (ex: CZ(0,1) puis CZ(1,2)), chaque transformation CZ doit voir les modifications 
    # des transformations précédentes. Une approche vectorisée qui extrait toutes les 
    # valeurs au début appliquerait les transformations en parallèle au lieu de 
    # séquentiellement, ce qui donnerait des résultats incorrects.
    if isinstance(control_qubits, int):
        control_qubits = [control_qubits]
    if isinstance(target_qubits, int):
        target_qubits = [target_qubits]
    assert len(control_qubits) == len(target_qubits)

    new_z = paulis.z.copy()
    new_x = paulis.x.copy()
    new_phase = paulis.phase.copy()

    for c, t in zip(control_qubits, target_qubits):
        zc = new_z[:, c]
        xc = new_x[:, c]
        zt = new_z[:, t]
        xt = new_x[:, t]

        add_phase = 2 * xc.astype(int) * xt.astype(int) * (zc ^ zt).astype(int)
        new_phase = np.mod(new_phase + add_phase, 4)

        new_z[:, c] = np.logical_xor(zc, xt)
        new_z[:, t] = np.logical_xor(zt, xc)

    return PauliList.from_symplectic(new_z, new_x, new_phase)