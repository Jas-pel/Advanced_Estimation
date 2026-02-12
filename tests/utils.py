import numpy as np
from numpy.typing import NDArray
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator, PauliList


def check_diag_transformation(paulis: PauliList, diag_paulis: PauliList, circuit: QuantumCircuit) -> bool:
    """
    Check if the circuit actually diagonalize the Paulis into the given diagonal Pauli

    Args:
        paulis (PauliList): _description_
        diag_paulis (PauliList): _description_
        circuit (QuantumCircuit): _description_

    Returns:
        bool: _description_
    """

    assert np.all(~diag_paulis.x)

    circuit_unitary = Operator(circuit).to_matrix()
    paulis_matrices = paulis.to_matrix(array=True)
    diag_paulis_matrices = diag_paulis.to_matrix(array=True)

    ref_diag_paulis_matrices = np.einsum("ij,lk,pjk->pil", circuit_unitary, circuit_unitary.conj(), paulis_matrices)

    return np.allclose(diag_paulis_matrices, ref_diag_paulis_matrices)


def check_expectation_values_within_range(
    exp_values: NDArray[np.float64],
    ref_exp_values: NDArray[np.float64],
    variances: NDArray[np.float64],
    sigma: int = 4,
    min_relative_distance=0.1,
):

    within_range = (exp_values - ref_exp_values) ** 2 < sigma * variances
    within_min_rel = (exp_values - ref_exp_values) ** 2 < (min_relative_distance * ref_exp_values) ** 2

    return np.sum(np.logical_or(within_range, within_min_rel).astype(int)) / within_range.size > 0.9
