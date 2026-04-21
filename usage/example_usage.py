import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp
from qiskit_aer import AerSimulator
from qiskit_ibm_runtime import SamplerV2 as Sampler

from advanced_estimation.commutation.base_commutation import BaseCommutation
from advanced_estimation.commutation.general_commuting import GeneralCommutation
from advanced_estimation.estimation.observable_estimation import (
    iterative_estimate_sparse_pauli_op_expectation_value,
)


def main() -> None:
    """
    Simple example of how to use the iterative estimation for a SparsePauliOp observable.
    """
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cx(0, 1)

    observable = SparsePauliOp.from_list(
        [
            ("ZZ", 0.7),
            ("XX", 0.2),
            ("ZI", 0.1),
        ]
    )

    simulator = AerSimulator()
    sampler = Sampler(mode=simulator)

    exp_vals, variances = iterative_estimate_sparse_pauli_op_expectation_value(
        observable=observable,
        state_circuit=qc,
        sampler=sampler,
        commutation_module=GeneralCommutation(),
        shots_budget=2000,
        num_iterations=5,
    )

    print("Iter | Expectation value | Variance")
    for i, (ev, var) in enumerate(zip(exp_vals, variances), start=1):
        print(f"{i:>4} | {ev:>17.8f} | {var:>10.8f}")

    print("\nDernière estimation:")
    print(f"  <O>  = {exp_vals[-1]:.8f}")
    print(f"  Var  = {variances[-1]:.8f}")
    print(f"  Std  = {np.sqrt(max(variances[-1], 0.0)):.8f}")


if __name__ == "__main__":
    main()
