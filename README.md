# Advanced Estimator

**Jasmin**

Un framework pour estimer efficacement la valeur moyenne d'observables quantiques en regroupant les chaînes de Pauli commutantes via la théorie des graphes, la diagonalisation de Clifford et l'allocation adaptative de shots.

---

## Table des matières

1. [Problème et motivation](#problème-et-motivation)
2. [Architecture du projet](#architecture-du-projet)
3. [Pipeline d'estimation](#pipeline-destimation)
4. [Module commutation](#module-commutation)
5. [Module clifford.py](#module-cliffordpy)
6. [Module estimation](#module-estimation)
7. [Scripts d'utilisation](#scripts-dutilisation)
8. [Installation](#installation)
9. [Utilisation rapide](#utilisation-rapide)
10. [Tests](#tests)
11. [Références](#références)

---

## Problème et motivation

En informatique quantique, **mesurer un état le détruit**. Pour estimer la valeur moyenne d'une observable

$$\langle \hat{A} \rangle = \sum_i c_i \langle P_i \rangle$$

il faut répéter les mesures (*shots*) un grand nombre de fois. Mesurer chaque chaîne de Pauli séparément est très coûteux (ex: 100 Pauli × 1000 shots = 100 000 circuits).

**L'idée clé** : regrouper les chaînes de Pauli qui **commutent** en **cliques**, les mesurer simultanément avec un seul circuit, et ainsi réduire drastiquement le nombre total de mesures nécessaires.

---

## Architecture du projet

```
advanced_estimator/
├── commutation/                  # Stratégies de groupement
│   ├── base_commutation.py       # Classe abstraite + recherche de cliques (NetworkX)
│   ├── no_commuting.py           # Aucun groupement (1 Pauli = 1 circuit)
│   ├── bitwise_commuting.py      # Commutation qubit-par-qubit
│   ├── general_commuting.py      # Commutation générale + diagonalisation Clifford
│   └── clifford.py               # Transformations symplectiques H, S, CX, CZ
├── estimation/
│   ├── pauli_estimation.py       # Estimation des ⟨Pᵢ⟩ et matrice de covariance
│   └── observable_estimation.py  # Estimation itérative de ⟨Â⟩ avec réallocation de shots
tests/                            # Tests unitaires
usage/                            # Scripts de comparaison et visualisation
```

---

## Pipeline d'estimation

```
Observable Â = Σᵢ cᵢ Pᵢ
        │
        ▼
┌──────────────────────────────┐
│  1. Graphe de commutation    │  Nœud = Pauli, Arête = commutent
│     commutation_table()      │
└──────────┬───────────────────┘
           ▼
┌──────────────────────────────┐
│  2. Cliques maximales        │  NetworkX : find_cliques()
│     find_commuting_cliques() │  Chaque clique = mesurable ensemble
└──────────┬───────────────────┘
           ▼
┌──────────────────────────────┐
│  3. Diagonalisation          │  Circuit Clifford (H, S, CX, CZ)
│     diagonalize_paulis_      │  Résultat : tous les Pauli → Z/I
│     with_circuit()           │
└──────────┬───────────────────┘
           ▼
┌──────────────────────────────┐
│  4. Mesure                   │  |ψ⟩ → circuit diag → mesure Z
│     Sampler Qiskit           │  Simulateur ou hardware réel
└──────────┬───────────────────┘
           ▼
┌──────────────────────────────┐
│  5. Statistiques             │  Valeurs moyennes ⟨Pᵢ⟩
│     + covariances            │  Matrice de covariance Cov(Pᵢ, Pⱼ)
└──────────┬───────────────────┘
           ▼
┌──────────────────────────────┐
│  6. Réallocation itérative   │  Plus de shots aux cliques
│     des shots                │  avec plus grande variance
└──────────────────────────────┘
           │
           ▼
        ⟨Â⟩ ± σ
```

---

## Module commutation

Trois stratégies de commutation, toutes héritant de `BaseCommutation` :

### `BaseCommutation` (classe abstraite)

- `find_commuting_cliques(paulis)` : construit le graphe de commutation et trouve les cliques maximales via `networkx.find_cliques()`
- Méthodes abstraites : `commutation_table()` et `diagonalize_paulis_with_circuit()`

### `NoCommutation`

- **Table de commutation** : matrice identité (rien ne commute sauf avec soi-même)
- **Diagonalisation** : H sur les X, S†H sur les Y, rien sur les Z
- **Résultat** : autant de circuits que de chaînes de Pauli
- **Usage** : baseline naïve pour comparer

### `BitwiseCommutation`

- **Condition** : Pᵢ et Pⱼ commutent **qubit par qubit** (*bitwise*)
- **Table** : `¬(z₁·x₂ ⊕ x₁·z₂)` sur chaque qubit, `all()` sur les qubits
- **Diagonalisation** : portes H et S† uniquement (single-qubit)
- Groupement conservateur mais simple

### `GeneralCommutation`

- **Condition** : commutation quantique standard

$$\sum_q \left( z_i^{(q)} x_j^{(q)} + x_i^{(q)} z_j^{(q)} \right) \equiv 0 \pmod{2}$$

- **Diagonalisation** : algorithme de Gokhale et al. (arXiv:1907.13623)
  1. Réduction en générateurs (forme échelon binaire)
  2. Packing diagonal sur les blocs X et Z
  3. Élimination avec H, CX, CZ, S
- **Le plus efficace** : produit le moins de cliques = le moins de circuits

---

## Module clifford.py

Transformations symplectiques des chaînes de Pauli dans le **tableau de Gottesman-Knill**. Opèrent directement sur la représentation binaire (Z | X | φ) en O(n) au lieu de manipuler des matrices 2ⁿ × 2ⁿ.

| Porte | Transformation Z | Transformation X | Mise à jour de la phase |
|-------|------------------|------------------|------------------------|
| **H** | Z ↔ X | X ↔ Z | φ + 2·z·x |
| **S** | Z ← Z ⊕ X | X inchangé | φ + 2·z·x |
| **CX** | Zc ← Zc ⊕ Zt | Xt ← Xt ⊕ Xc | φ + 2·xc·zt·(1 − zc ⊕ xt) |
| **CZ** | Zc ← Zc ⊕ Xt, Zt ← Zt ⊕ Xc | inchangé | φ + 2·xc·xt·(zc ⊕ zt) |

Ces fonctions sont utilisées par `GeneralCommutation.diagonalize_paulis_with_circuit()` pour transformer les Pauli en Pauli diagonaux tout en construisant le circuit correspondant.

---

## Module estimation

### `pauli_estimation.py`

| Fonction | Rôle |
|----------|------|
| `bitstrings_to_bits()` | Convertit les bitstrings de mesure (`"011"`) en matrice booléenne numpy |
| `diag_paulis_expectation_values_and_covariances()` | Calcule ⟨Pᵢ⟩ et Cov(Pᵢ, Pⱼ) à partir des counts de mesure |
| `estimate_cliques_expectation_values_and_covariances()` | Orchestre mesure + statistiques pour chaque clique séparément |
| `overall_paulis_expectation_values_and_covariances()` | Combine les résultats de toutes les cliques (moyenne pondérée par shots) |
| `get_paulis_shots()` | Calcule le nombre total de shots alloués à chaque Pauli |

### `observable_estimation.py`

| Fonction | Rôle |
|----------|------|
| `iterative_estimate_sparse_pauli_op_expectation_value()` | Estimation itérative de ⟨Â⟩ avec réallocation adaptative des shots selon la variance pondérée par clique |
| `compute_weighted_cliques_variances()` | Calcule Σⱼₖ cⱼ · Cov(Pⱼ, Pₖ) · cₖ pour chaque clique |

**Algorithme itératif** :
1. Distribuer les shots uniformément entre les cliques
2. Estimer ⟨Pᵢ⟩ et les covariances
3. Calculer la variance pondérée de chaque clique
4. Réallouer les shots proportionnellement à √(variance) de chaque clique
5. Répéter → convergence vers l'allocation optimale

---

## Scripts d'utilisation

| Script | Description |
|--------|-------------|
| `usage/compare_paulis_estimations.py` | Compare les 3 stratégies (No / Bitwise / General) sur l'estimation de ⟨Pᵢ⟩ individuels. Produit des barres d'erreur et des heatmaps de covariance. |
| `usage/compare_observable_estimations.py` | Compare les 3 stratégies sur l'estimation itérative de ⟨Â⟩. Montre la convergence avec barres d'erreur par itération. |
| `usage/compare_paulis_bitwise_general.py` | Comparaison directe Bitwise vs General sur un état aléatoire. |

---

## Installation

```bash
# Cloner le dépôt
git clone <url-du-repo>
cd advanced-estimation-jasmin_thomas

# Créer et activer l'environnement virtuel
python -m venv venv

# Windows PowerShell :
venv\Scripts\Activate.ps1

# Linux / Mac :
source venv/bin/activate

# Installer les dépendances
pip install -r requirements.txt
```

### Dépendances

| Package | Rôle |
|---------|------|
| `numpy` | Calcul numérique et opérations symplectiques |
| `networkx` | Graphes de commutation et recherche de cliques |
| `qiskit` | Circuits quantiques et représentation des Pauli |
| `qiskit-aer` | Simulateur quantique local |
| `qiskit-ibm-runtime` | Interface Sampler V2 |

---

## Utilisation rapide

```python
import numpy as np
from qiskit import transpile
from qiskit.circuit.random import random_circuit
from qiskit.quantum_info import SparsePauliOp, pauli_basis
from qiskit_aer import AerSimulator
from qiskit_ibm_runtime import SamplerV2 as Sampler

from advanced_estimator.commutation import GeneralCommutation
from advanced_estimator.estimation.observable_estimation import (
    iterative_estimate_sparse_pauli_op_expectation_value,
)

# Setup
simulator = AerSimulator()
sampler = Sampler(mode=simulator)

# Observable aléatoire sur 3 qubits
paulis = pauli_basis(3)
coeffs = 2 * np.random.random(paulis.size) - 1
observable = SparsePauliOp(paulis, coeffs)

# Circuit d'état aléatoire
state_circuit = transpile(random_circuit(3, depth=4), simulator)

# Estimation avec commutation générale
exp_values, variances = iterative_estimate_sparse_pauli_op_expectation_value(
    observable, state_circuit, sampler,
    commutation_module=GeneralCommutation(force_single_qubit_generators=True),
    shots_budget=2000,
    num_iterations=5,
)

print(f"<A> = {exp_values[-1]:.4f} +/- {np.sqrt(variances[-1]):.4f}")
```

---

## Tests

```bash
# Tous les tests
pytest tests/

# Un fichier spécifique
pytest tests/test_clifford.py
pytest tests/test_general.py
pytest tests/test_bitwise.py
pytest tests/test_estimation.py
```

---

## Références

- Gokhale et al., *Minimizing State Preparations in Variational Quantum Eigensolver by Partitioning into Commuting Families*, arXiv:1907.13623 (2019)
- Aaronson & Gottesman, *Improved Simulation of Stabilizer Circuits*, Phys. Rev. A 70, 052328 (2004)
- [Documentation Qiskit](https://qiskit.org)
