from qiskit.quantum_info import PauliList

from advanced_estimation.commutation import clifford


def test_h():

    paulis = PauliList(["XX", "ZZ", "YY", "XZ"])
    ref_paulis = PauliList(["ZX", "XZ", "-YY", "ZZ"])

    new_paulis = clifford.h(paulis, [1])

    assert new_paulis == ref_paulis


def test_s():

    paulis = PauliList(["XX", "ZZ", "YY", "XY"])
    ref_paulis = PauliList(["YY", "ZZ", "XX", "-YX"])

    new_paulis = clifford.s(paulis, [0, 1])

    assert new_paulis == ref_paulis


def test_cx():

    paulis = PauliList(["XX", "ZZ", "YY", "XZ"])
    ref_paulis = PauliList(["IX", "ZI", "-ZX", "XZ"])

    new_paulis = clifford.cx(paulis, 0, 1)

    print(new_paulis)

    assert new_paulis == ref_paulis
