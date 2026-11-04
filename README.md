# Composite Laminate Analysis (CLT) in MATLAB

A MATLAB port of my Python toolkit
[composite-fuselage-laminate-analysis](https://github.com/Sanjanakamboj/composite-fuselage-laminate-analysis):
Classical Lamination Theory for a composite fuselage skin panel, from lamina stiffness to
ply-level failure assessment.

```
material -> Q / Qbar -> laminate ABD -> mid-plane strain + curvature -> ply stress recovery
         -> Max-Stress / Tsai-Hill / Tsai-Wu failure indices -> baseline laminate comparison
```

> **Scope.** This port covers the analysis chain above and the baseline laminate-family
> comparison. The bounded design-space search and the mass-reduction study stay in the
> Python repo. Loads and strengths are illustrative, preliminary-design-level values, not
> certified allowables.

## Contents

| File | Purpose |
|---|---|
| `src/reduced_stiffness.m` | Plane-stress reduced stiffness `Q` (ply 1-2 axes) |
| `src/transformed_reduced_stiffness.m` | `Qbar` in global x-y axes (Jones/Reddy closed form) |
| `src/stress_transform_matrix.m`, `src/transform_strain_to_material.m` | Stress/strain rotation with the Reuter matrix (engineering shear strain) |
| `src/laminate_abd.m` | A, B, D assembly for a bottom-to-top stacking sequence, areal mass |
| `src/solve_deformation.m` | `[eps0; kappa] = ABD \ [N; M]`, with a condition-number guard |
| `src/recover_ply_states.m` | Strain/stress at the bottom and top of every ply |
| `src/failure_indices.m` | Maximum Stress, Tsai-Hill, Tsai-Wu (default `F12 = -0.5*sqrt(F11*F22)`) |
| `src/run_baseline_study.m` | Four symmetric laminate families under a membrane load |
| `examples/fuselage_baseline_study.m` | Reproduces the baseline comparison table |
| `tests/run_tests.m` | Verification against the Python implementation |

Conventions follow the Python original: ply angle in degrees from global +x to the fibre axis,
counter-clockwise positive; stress order `[sigma_x sigma_y tau_xy]`; engineering shear strain.

## Run

```matlab
run_tests                      % from the tests/ folder, or add tests/ to the path
run('examples/fuselage_baseline_study.m')
```

Plain MATLAB, no toolboxes. It also runs unchanged in GNU Octave.

## Verification

`tests/run_tests.m` checks the port against reference values produced by the Python package
(same material, strengths and load case):

* `Q` and `Qbar(30 deg)` match to 1e-12 (relative)
* A11, A22, A66 for all four laminate families match to 1e-9
* Max-Stress, Tsai-Hill and Tsai-Wu failure indices match to 1e-9
* closed-form identities: `Qbar(0) = Q`, +45/-45 coupling-term symmetry, strain round trip,
  `B = 0` for a symmetric laminate, rejection of a singular ABD

## Example output

Representative carbon/epoxy ply (E1 = 161 GPa, E2 = 11.4 GPa, G12 = 5.17 GPa, nu12 = 0.32),
membrane load `N = [98 000; 116 000; 20 000] N/m`:

| Laminate | A11 (MN/m) | A22 (MN/m) | A66 (MN/m) | Max-Stress FI | Tsai-Hill FI | Tsai-Wu FI |
|---|---:|---:|---:|---:|---:|---:|
| quasi-isotropic [0/45/-45/90]s | 101.02 | 101.02 | 34.41 | 0.241 | 0.059 | 0.155 |
| 0-biased [0/0/45/-45]s | 156.47 | 45.56 | 34.41 | 0.478 | 0.230 | 0.392 |
| 90-biased [90/90/45/-45]s | 45.56 | 156.47 | 34.41 | 0.386 | 0.150 | 0.291 |
| +-45-biased [45/-45/45/-45]s | 74.22 | 74.22 | 61.20 | 0.224 | 0.056 | 0.140 |

## License

MIT, see `LICENSE`.
