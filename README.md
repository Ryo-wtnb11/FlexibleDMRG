# FlexibleDMRG

DMRG implementation including the following features:
 - search of the optimal site ordering by swapping local sites.

## How to use the package
1. Assuming you have Julia installed. You move to the repo directory and call `julia --project` to start the REPL.

2. You first need to run the code in the REPL:

```
using Pkg
Pkg.instantiate()
```

## Examples

### Basic DMRG with site ordering optimization
You can try the example for the dimer covering in the `examples/dimer_covering.jl` file.

```
> julia --project
(FlexibleDMRG) julia> include("examples/dimer_covering.jl")
```

### XXZ Model Hamiltonians
The package now includes convenient functions for creating XXZ model Hamiltonians:

```julia
using FlexibleDMRG

# 1. Nearest neighbor XXZ model
H = create_nearest_neighbor_xxz(1/2, 10, 1.0, 0.5; periodic=false)

# 2. XXZ model with custom couplings from data
J_data = [(1, 2, 1.0), (2, 3, 1.0), (3, 4, 1.0)]
Δ_data = [(1, 2, 0.5), (2, 3, 0.5), (3, 4, 0.5)]
H = create_xxz_hamiltonian_from_data(1/2, 4, J_data, Δ_data)

# 3. XXZ model with full coupling matrices
sites = siteinds("S=1/2", 4)
J_ij = [0.0 1.0 0.0 0.0;
        1.0 0.0 1.0 0.0;
        0.0 1.0 0.0 1.0;
        0.0 0.0 1.0 0.0]
Δ_ij = [0.0 0.5 0.0 0.0;
        0.5 0.0 0.5 0.0;
        0.0 0.5 0.0 0.5;
        0.0 0.0 0.5 0.0]
H = create_xxz_hamiltonian(sites, J_ij, Δ_ij)

# 4. XXZ model from data file
H = create_xxz_hamiltonian_from_file(1/2, 6, "couplings.dat")
```

#### Data File Format
For the file-based approach, create a text file with the format:
```
i j J_{i,j} Δ_{i,j}
```
where `i` and `j` are site indices (1-based), and `J_{i,j}` and `Δ_{i,j}` are the coupling parameters.

Example file (`couplings.dat`):
```
# Nearest neighbor couplings
1 2 1.0 0.5
2 3 1.0 0.5
3 4 1.0 0.5
4 5 1.0 0.5
5 6 1.0 0.5

# Next-nearest neighbor couplings
1 3 0.3 0.2
2 4 0.3 0.2
3 5 0.3 0.2
4 6 0.3 0.2
```

Try the examples:
```
(FlexibleDMRG) julia> include("examples/xxz_example.jl")
(FlexibleDMRG) julia> include("examples/xxz_from_file_example.jl")
(FlexibleDMRG) julia> include("examples/simple_xxz_example.jl")
```

### Parameter File System
For complex calculations, you can use parameter files to configure all aspects of the DMRG calculation:

```toml
# parameters.toml
[system]
spin_size = 0.5
system_size = 6

[hamiltonian]
coupling_file = "couplings.dat"

[dmrg]
max_bond_dim = 50
energy_tolerance = 1e-8
convergence_sets = 3
max_sweeps = 10
max_sets = 20

[output]
filename = "results"
folder = "output"
```

Run with parameter file:
```bash
julia run_with_parameters.jl parameters.toml
```

Below is the minimal code base for augumented DMRG method based on the site ordering optimization.
Since relying on the [ITensorMPS.jl](https://github.com/ITensor/ITensorMPS.jl), you should take a deep read of the docs at least of `dmrg()` and `OpSum()` functions.

```julia
using FlexibleDMRG

H = ... # Construct the MPO with `OpSum()`.
psi0 = ... # Construct the MPS with `MPS()`.
sweeps = ... # Construct the `Sweeps()` object.
dmrg_kwargs = ... # Construct the `kwargs` dictionary of `dmrg()`.

kwargs = (
    nswapsites = 3, # Number of sites to swap at each step.
    proposer = MinimumEESiteOrderProposer(), # Proposer of the site ordering.
)

energy, psi, whole_order = siteordering(H, psi0, sweeps; kwargs..., dmrg_kwargs=dmrg_kwargs)
```
If you put the `nswapsites` as 2, an order is optimized just by adjacent swaps as Ref. [1].
You can choose `proposer` from `src/siteorderproposer.jl`. The default one is `MinimumEESiteOrderProposer()`, then the order with low entanglement entropy (EE) is proposed as Ref. [2, 3].

## Upcoming features
- Augumented by Clifford + Matchgate in Ref. [4, 5].

## References
1. [W. Li et al. J. Phys.: Condens. Matter 34 254003 (2022)](https://iopscience.iop.org/article/10.1088/1361-648X/ac640e/meta)
2. [T. Hikihara et al. Phys. Rev. Research 5, 013031 (2023)](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.5.013031)
3. [R. Watanabe et al. arXiv:2505.05908 (2025)](https://arxiv.org/abs/2505.05908)
4. [X. Qian et al. Phys. Rev. Lett. 133, 190402 (2024)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.133.190402)
5. [J. Huang et al. arXiv:2505.08635 (2025)](https://arxiv.org/abs/2505.08635)

