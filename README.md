# FlexibleDMRG

DMRG implementation including the following features:
 - reconfiguration of the site ordering.

## How to use the package
1. Assuming you have Julia installed. You move to the repo directory and call `julia --project` to start the REPL.

2. You first need to run the code in the REPL:

```
using Pkg
Pkg.instantiate()
```

## Example
You can try the example that is the Heisenberg dimer covering problem in the `examples/dimer_covering.jl` file.

```
> julia --project
(FlexibleDMRG) julia> include("examples/dimer_covering.jl")
```

Below is the minimal code for the usage.
Since we rely on the [ITensorMPS.jl](https://github.com/ITensor/ITensorMPS.jl), you should deep read the documentations at least of `dmrg()` and `OpSum()` functions.

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

