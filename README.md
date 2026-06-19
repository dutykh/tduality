# tduality

Companion code and data for the paper

> **Finite-size T-duality electrodynamics in (2+1) dimensions**
> Davide Batic and Denys Dutykh, Mathematics Department, Khalifa University of
> Science and Technology, Abu Dhabi, United Arab Emirates (2026).

This repository contains the Julia scripts, the numerical output, and the
manuscript for a finite-size formulation of T-duality-deformed electrostatics
in 2+1 dimensions: corrected Green functions, self-energy-subtracted pair
potentials, boundary effects, and the running interaction dimension. The
manuscript is included in [`pdf/T_Duality.pdf`](pdf/T_Duality.pdf).

## Background

In two spatial dimensions the static Green function of the Laplacian is
logarithmic, so the ordinary Coulomb interaction carries an unavoidable
infrared (large-distance) ambiguity and a short-distance divergence. The
T-duality (zero-point-length) deformation replaces the Maxwell kinetic operator
by a nonlocal form factor

$$f(k) = \ell_0\,k\,K_1(\ell_0 k),$$

where $K_1$ is a modified Bessel function and $\ell_0$ is the zero-point length.
This factor smooths the ultraviolet (short-distance) behaviour while leaving the
infrared sector to be fixed by the *geometry*. The paper makes that split
explicit:

- **Corrected infinite plane.** The T-duality kernel is equivalent to sourcing a
  point charge by a finite, normalised charge cloud
  $\rho_{\ell_0}(r) = (\pi\ell_0^2)^{-1}(1 + r^2/\ell_0^2)^{-2}$. The
  self-energy-subtracted neutral-pair energy is finite at the origin and recovers
  the logarithmic Coulomb law at large distance:

$$\Delta V_{\infty,\ell_0}(r) = \frac{Q^2}{4\pi}\,\ln\!\left(1 + \frac{r^2}{\ell_0^2}\right).$$

- **Bounded / finite domains.** The Green function is defined through the
  *spectral calculus* of the Laplacian, $\mathcal{G} = F_{\ell_0}(A_B)\,A_B^{-1}P_\perp$,
  so that the boundary condition is imposed *before* the nonlocal operator is
  formed. The paper works out three canonical realisations:
  - **square torus** — charge neutrality via zero-mode subtraction;
  - **grounded half-plane** — Dirichlet boundary, image-type self-energy;
  - **grounded strip** — two Dirichlet boundaries, closed image-summed form.

- **Running interaction dimension.** An additive-constant-free diagnostic
  $\mathbb{D}_\Delta(r) = 3 - \partial \ln \mathcal{V}(r)/\partial \ln r$ built from
  the neutral-pair energy. The ultraviolet flow to one dimension survives in
  finite volume, while the infrared behaviour is geometry-dependent (the torus
  saturates; grounded boundaries screen the logarithm).

The two scripts in this repository reproduce **Figure 1** (the finite-volume
T-duality pair potential on the square torus) and **Figure 2** (the running
interaction dimension on the square torus), together with the data and
convergence checks behind them.

## Repository layout

```
.
├── codes/                                       Julia scripts (self-contained)
│   ├── plot_torus_potential_pdfjl               Fig. 1 — axial torus pair potential
│   └── plot_running_dimension_torus_noqt_fixed.jl   Fig. 2 — running interaction dimension
├── convergence/                                 Reference output and convergence checks
│   ├── torus_potential_a_0p020.csv
│   ├── torus_potential_a_0p050.csv
│   ├── torus_potential_a_0p100.csv
│   ├── convergence_report.txt
│   ├── running_dimension_a_0p020.csv
│   ├── running_dimension_a_0p050.csv
│   ├── running_dimension_a_0p100.csv
│   └── running_dimension_convergence_report.txt
├── pdf/
│   └── T_Duality.pdf                            The manuscript
├── LICENSE                                      GNU LGPL v2.1
└── README.md
```

> **Note on a filename.** The first script is committed as
> `codes/plot_torus_potential_pdfjl` — the dot before the extension is missing
> (it is conceptually `plot_torus_potential_pdf.jl`). Julia runs a file
> regardless of its name, so the commands below work as written.

## Requirements

- [Julia](https://julialang.org/) (developed and tested with 1.12).
- Internet access on first run, so the scripts can fetch their one dependency.

The scripts are deliberately **self-contained and headless**. They avoid
`Plots.jl`, `GR.jl`, Qt, Makie, GUI backends, and `pdflatex`. Each script
creates an isolated Julia environment in a hidden folder next to itself
(`.torus_potential_env_noqt` / `.running_dimension_env_noqt`), installs only
[`SpecialFunctions.jl`](https://github.com/JuliaMath/SpecialFunctions.jl) for the
modified Bessel functions $K_0, K_1$, and writes a hand-rolled vector PDF
directly. Nothing else in your Julia setup is touched.

## Running

```bash
cd codes

# Figure 1 — finite-volume T-duality potential on the square torus
julia plot_torus_potential_pdfjl

# Figure 2 — potential-based running interaction dimension on the torus
julia plot_running_dimension_torus_noqt_fixed.jl
```

Each script prints its progress, then writes its results to an output directory
created alongside it:

| Script | Output directory | Files |
| --- | --- | --- |
| `plot_torus_potential_pdfjl` | `codes/torus_potential_output_pdf/` | `torus_potential.pdf`, three `torus_potential_a_*.csv`, `convergence_report.txt` |
| `plot_running_dimension_torus_noqt_fixed.jl` | `codes/running_dimension_output_noqt/` | `running_dimension_torus.pdf`, three `running_dimension_a_*.csv`, `running_dimension_convergence_report.txt` |

The CSV files and convergence reports archived under [`convergence/`](convergence/)
are the reference outputs of these runs and can be compared against your own.

### Parameters

Both scripts sweep the dimensionless ratio $\mathfrak{a} = \ell_0/L$ over
`A_VALUES = [0.02, 0.05, 0.10]` and the dimensionless separation
$\xi = r/L \in [0, 0.5]$ on a grid of 501 points (the axial path on the torus).
These, and the convergence tolerances, are `const`s near the top of each script.

## Output formats

**`torus_potential_a_*.csv`** — columns

```
r_over_L, torus_deltaV_over_Q2, infinite_plane_deltaV_over_Q2
```

i.e. the separation $\xi = r/L$, the zero-mode-subtracted toroidal pair energy
$\Delta V_{L,\ell_0}/Q^2$, and the corrected infinite-plane reference
$\frac{1}{4\pi}\ln(1 + \xi^2/\mathfrak{a}^2)$.

**`running_dimension_a_*.csv`** — columns

```
xi, D_torus, D_infinite
```

i.e. the separation $\xi$, the toroidal running dimension
$\mathbb{D}_{\Delta,\mathbb{T}^2}(\xi;\mathfrak{a})$, and the infinite-plane
reference $\mathbb{D}_\Delta^{(\infty)}(\xi/\mathfrak{a})$.

**Convergence reports** record, for each $\mathfrak{a}$, the symmetric lattice
cutoff $N$ ($-N \le n_1,n_2 \le N$, zero mode omitted), the analytic tail
estimate $K_0(2\pi\mathfrak{a}N)/\pi$, and the direct $N$-vs-$2N$ difference on
the plotting grid. The tail tolerance is $10^{-10}$ for the potential; the
dimension curve additionally targets $10^{-8}$ on the $N$-vs-$2N$ difference.
In the archived runs the direct differences are at or below the printed
double-precision floor, confirming that the visible deviations from the
infinite-plane reference are genuine finite-size effects rather than truncation
artefacts.

## Key formulas implemented

Square-torus axial pair potential (Eq. 81 in the paper), with
$\mathfrak{a} = \ell_0/L$ and $\xi = r/L$:

$$\mathcal{U}_{\mathrm{ax}}(\xi;\mathfrak{a}) = \sum_{\mathbf{n}\neq\mathbf{0}}
\frac{\mathfrak{a}\,K_1(2\pi\mathfrak{a}|\mathbf{n}|)}{2\pi|\mathbf{n}|}
\bigl[1 - \cos(2\pi n_1 \xi)\bigr].$$

Potential-based running interaction dimension (Eq. 123 / 135):

$$\mathbb{D}_\Delta(\xi;\mathfrak{a}) = 3 - \frac{\xi\,\partial_\xi \mathcal{U}}{\mathcal{U}},$$

evaluated analytically from the truncated Fourier series rather than by finite
differences. Infinite-plane references:
$\mathcal{U}_\infty = \frac{1}{4\pi}\ln(1 + \xi^2/\mathfrak{a}^2)$ and
$\mathbb{D}_\Delta^{(\infty)}(s) = 3 - 2s^2/[(1+s^2)\ln(1+s^2)]$ with
$s = \xi/\mathfrak{a}$.

## Citing

If you use this code or data, please cite the paper:

> D. Batic and D. Dutykh, *Finite-size T-duality electrodynamics in (2+1)
> dimensions* (2026).

## License

This project is released under the **GNU Lesser General Public License v2.1**.
See [`LICENSE`](LICENSE) for the full text.

## Authors

- Davide Batic — `davide.batic@ku.ac.ae`
- Denys Dutykh — `denys.dutykh@ku.ac.ae`

Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi,
United Arab Emirates.
