# Lithium inventory checks for PyBaMM issue #5700 and PR #5524

Supporting material for a comment on https://github.com/pybamm-team/PyBaMM/issues/5700.
It contains the inventory code, a wrapper that reproduces the headline numbers, the
per-run outcome table for all 84 configurations, and the raw time series for the
configurations quoted in the comment.

## Revisions checked

| Label | Source | Commit |
|---|---|---|
| baseline | tag v26.4.1 | `d9e5afbd1c5f0b05b6f5ea521c2aab3a4c3c3c81` |
| main | main branch, 2026-09-05 | `258fdc8c125869d2598757c2bc4ff1e87b60cb3b` |
| pr | PR #5524 head, 2026-09-03 | `719763b6925dd7a9b2c69932df37ed4057799379` |

Each revision was installed from the GitHub source archive for that commit. The
installed package reports version 0.0.0 because it was built from an archive, so
identity is the commit, not the version string.

## Environment

Python 3.12.13, CasADi 3.7.2, NumPy 2.5.2, SciPy 1.18.1, Pint 0.25.3, SymPy 1.14.0.
Native solver: pybammsolvers 0.6.0 for v26.4.1 (its declared range) and 0.9.0 for
main and the PR. Exact pins are in `requirements-baseline.lock` and
`requirements-main.lock`; `requirements-common.txt` lists the shared packages.
Because the native solver differs between v26.4.1 and the current revisions,
baseline-to-current comparisons carry a solver-version difference. Main and PR
share an identical environment.

## Cases

All cases run the protocol from the issue description: rest 10 min, 1C discharge
20 min, rest 10 min, 1C charge 20 min, rest 20 min, 5 s period, initial SOC 0.45.

| Case | Model | Parameters |
|---|---|---|
| basic_variable | `BasicDFN` | ORegan2022, negative particle diffusivity 3.3e-14 |
| basic_constant | `BasicDFN` | as above, transference number fixed at its value for the initial concentration and ambient temperature |
| half | `DFN({"working electrode": "positive"})` | Xu2019 |
| dfn_control | `DFN()` | ORegan2022, negative particle diffusivity 3.3e-14 |

Solver: `IDAKLUSolver(rtol=tol, atol=tol, root_method="casadi", root_tol=1e-10,
options={"num_threads": 1})`. Mesh: `nx` cells per x-domain, `nr` shells per
particle. Schedule: (nx, nr, tol) in (20,20,1e-6), (20,20,1e-8), (20,20,1e-10),
(40,40,1e-10), (80,80,1e-10), (80,20,1e-10), (20,80,1e-10), for each case and
revision. 84 runs; 73 completed; 11 half-cell runs failed inside the solver and
are listed with their error strings in `summary.csv`.

## Inventory definition

N(t) = A [ sum over electrolyte domains of eps_k * sum_i c_e,k,i * dx_i
        + sum over electrodes of eps_s,e * sum_i dx_i * sum_j w_j * c_s,e,j,i ]

where A is electrode width times height times the number of parallel electrodes,
dx_i are the finite-volume cell widths from the mesh edges, and
w_j = (r_{j+1}^3 - r_j^3) / R^3 are spherical-shell volume fractions from the
particle mesh edges. Concentrations are read at the mesh nodes. The sum uses
PyBaMM's mesh, parameters and solution; it does not use `x_average`, `r_average`,
`Integral` or `Total lithium [mol]`. A second implementation with explicit loops
(`referee.inventory_loop`) is compared against the vectorised one on every run.

Full cell: R(t) = N(t) - N(0). Positive half-cell: R(t) = N(t) - N(0) - q(t) with
q(t) the integral of I/F over each constant-current step (metal is outside the
inventory; positive discharge current delivers lithium into the electrolyte).

`referee.controls()` runs the inventory on synthetic arrays with known answers
and confirms that an altered inventory, a reversed boundary sign, an omitted
source and a wrong volume factor are each detected. These are arithmetic checks
of the inventory code, not battery results.

## Reproducing the headline numbers

With one of the three PyBaMM revisions importable (for example a checkout of the
commit installed with `pip install -e`), and this directory as the working
directory:

```sh
PYBAMM_DISABLE_TELEMETRY=true python minimal_reproducer.py --output reviewer-output
```

This runs basic_variable, basic_constant and half at nx = nr = 20, tol = 1e-8 and
prints, per case, the maximum absolute residual, the residual relative to N(0),
and the final residual. Expected on main or the PR head:

| Case | max abs R (mol) | max abs R / N(0) | final R (mol) |
|---|---:|---:|---:|
| basic_variable | 6.518391162e-3 | 2.441742013e-2 | +6.518391162e-3 |
| basic_constant | 4.996003611e-15 | 1.871466687e-14 | +4.940492460e-15 |
| half | 1.637627307e-11 | 1.678862993e-7 | -2.187574618e-12 |

Each case writes `series.csv` plus one NPZ file per protocol step with the raw
arrays. Runtime is under a minute on a laptop.

## Files

- `measurement.py`: builds the simulation, reads fields at the mesh nodes, sums the inventory, integrates the current, and writes the per-step arrays and `series.csv`.
- `referee.py`: the two inventory implementations and the synthetic controls.
- `minimal_reproducer.py`: wrapper for the three 20/20, 1e-8 cases.
- `summary.csv`: one row per run (84), with residuals, rates, closure checks, runtime, and the error string for failed runs. `summary.json` has the full metadata per run, including solver and model options.
- `common-results/<run>/series.csv`: time series for the 20/20, 1e-8 runs on all three revisions for basic_variable and half, plus the PR basic_constant and dfn_control controls. Columns: step, time_s, inventory_mol, input_mol, residual_mol, inventory_rate_mol_s, input_rate_mol_s, balance_rate_mol_s, predicted_rate_mol_s, predicted_integral_mol, loop_inventory_mol, supplied_inventory, current_A, capacity_Ah. `supplied_inventory` is PyBaMM's `Total lithium [mol]` for the half-cell and the issue's fallback expression for BasicDFN.
- `common-results/pr_*_x20_r20_t1e-08/step[0-4].npz`: raw arrays for the two PR runs quoted in the comment. `ce_*` are mol/m3 at x cell centres; `cs_*` are mol/m3 with axes (radial shell, x cell, time); `dx_*` and `r_*` are metres; `electrolyte_flux` is mol/m2/s at the x cell edges; `area_m2` includes parallel electrodes.
- `main-pr-series-comparison.json`: byte-identity of `series.csv` for the 24 main/PR pairs that completed on both.
- `boundary-type-check.json`: count of boundary-condition types declared by each model at the PR head.
- `source-comparison.json`: SHA-256 of the model files that define the equations, per revision.
- `controls.json`: output of `referee.controls()`.
- `conservation.png`: residual time courses at the PR head, 20/20, 1e-8.
- `PACKAGE-SHA256.json`: SHA-256 of every other file here.

Do not integrate `balance_rate_mol_s` across a step boundary with a single
trapezoid; the current is discontinuous there. Integrate within each step.
