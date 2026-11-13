# Original Dycore (from GFDL)

This branch keeps the **original dry spectral dynamical core** that the moist
model on [`main`](https://github.com/PeterChang1023/Moist_Dycore) was built
from. It is here as a reference point: no moisture, no precipitation, only the
dynamics and Held–Suarez forcing.

The code is `JGCM` (`IdealizedSpectralGCM.jl`), a Julia implementation by
Daniel Zhengyu Huang, originally published under CliMA, that follows the
classical idealized spectral atmosphere model developed at the Geophysical Fluid
Dynamics Laboratory (GFDL).

## What is in this branch

```
IdealizedSpectralGCM.jl/
├── Project.toml            Julia package definition (package name: JGCM)
├── src/
│   ├── JGCM.jl             module entry point
│   ├── Atmos_Spectral/     spectral transforms, dynamics, time stepping, output
│   └── Atmos_Param/        Held–Suarez forcing
└── exp/
    ├── HSt21/              Held–Suarez test at T21, 10 levels
    ├── HSt42/              Held–Suarez test at T42, 20 levels
    ├── Barotropic/         barotropic vorticity equation
    └── Shallow_Water/      shallow-water equations
```

## Requirements

- **Julia.** Tested with Julia 1.10.4.
- **Python with matplotlib.** The model draws its figures through `PyPlot`. If
  Julia cannot find a Python with matplotlib, `PyPlot` installs a private one
  the first time it is built.

## Setup

Get this branch and install the Julia dependencies (once):

```shell
git clone -b Original_Dycore_from_GFDL https://github.com/PeterChang1023/Moist_Dycore.git
cd Moist_Dycore/IdealizedSpectralGCM.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

This downloads FFTW, JLD2, MAT, PyPlot and Revise, and takes a few minutes.

## Running the experiments

Each experiment is run from its own folder, because output files are written to
the current directory. `--project=../..` points Julia at the package two levels
up.

### Held–Suarez, T21 (quick test)

```shell
cd exp/HSt21
julia --project=../.. Run_HS.jl
```

As shipped this integrates 2 days and finishes in about 20 seconds. It prints
one line per model day:

```
[ Info: ("Day: ", 1, " Max |U|,|V|,|P|,|T| : ", 1.205, 0.769, 95150.2, 272.69)
[ Info: ("Day: ", 2, " Max |U|,|V|,|P|,|T| : ", 1.782, 0.841, 95285.9, 279.68)
```

and writes:

| File | Content |
|---|---|
| `Contourf_T.png`, `Contourf_Teq.png`, `Contourf_U.png`, `Contourf_V.png` | Zonal-mean temperature, equilibrium temperature, and winds |
| `HS_mean.dat` | Zonal-mean `T`, `u`, `v` (JLD2 format) |
| `HS_OpM.dat` | The full output manager (JLD2 format) |

### Held–Suarez, T42 (climate run)

```shell
cd exp/HSt42
julia --project=../.. Run_HS.jl
```

As shipped this integrates 1200 days and discards the first 200 as spin-up, the
standard Held–Suarez protocol. It takes far longer than the T21 test. The
`Contourf_*.png` files already in `exp/HSt42/` are the result of this run.

### Barotropic model

```shell
cd exp/Barotropic
julia --project=../.. Run_Barotropic.jl
```

Integrates 8 days (about 25 seconds) and writes `Barotropic_vel_u.png` and
`Barotropic_vor.png`.

### Shallow-water model

```shell
cd exp/Shallow_Water
julia --project=../.. Run_SW.jl
```

Integrates 8 days (about 1 minute) and writes `Shallow_Water_h.png`.

## Changing a run

**Length.** Edit `end_day` and `spinup_day` at the top of `Run_HS.jl`. Output is
averaged over the days after `spinup_day`.

**Resolution.** Edit this line in `HS.jl`:

```julia
num_fourier, nθ, nd = 21, 32, 10   # T21, 32 latitudes, 10 levels
num_fourier, nθ, nd = 42, 64, 20   # T42, 64 latitudes, 20 levels
```

The time step `Δt` (1200 s) is set a few lines below it and should be reduced
when the resolution is increased.

**Held–Suarez forcing.** The parameters are passed in from `Run_HS.jl`:

```julia
physics_params = Dict{String,Float64}(
    "σ_b"  => 0.7,        # top of the boundary layer in sigma
    "k_f"  => 1.0,        # Rayleigh friction rate (1/day)
    "k_a"  => 1.0/40.0,   # free-atmosphere thermal relaxation rate (1/day)
    "k_s"  => 1.0/4.0,    # surface thermal relaxation rate (1/day)
    "ΔT_y" => 60.0,       # equator-to-pole temperature difference (K)
    "Δθ_z" => 10.0,       # vertical potential temperature difference (K)
)
```

## Troubleshooting

- **No display, or a remote machine.** Tell matplotlib not to open windows:

  ```shell
  MPLBACKEND=Agg julia --project=../.. Run_HS.jl
  ```

- **`Package JGCM not found`.** The `--project` path is wrong. It must point at
  the folder that contains `Project.toml`.

## References

1. GFDL idealized spectral models quickstart:
   https://www.gfdl.noaa.gov/idealized-spectral-models-quickstart/
2. Held, I. M., and M. J. Suarez, 1994: A proposal for the intercomparison of
   the dynamical cores of atmospheric general circulation models. *Bull. Amer.
   Meteor. Soc.*, **75**, 1825–1830.
3. Ehrendorfer, M., 2011: *Spectral Numerical Weather Prediction Models*. SIAM.
4. Durran, D. R., 2013: *Numerical Methods for Wave Equations in Geophysical
   Fluid Dynamics*. Springer.
5. Lauritzen, P. H., et al., eds., 2011: *Numerical Techniques for Global
   Atmospheric Models*. Springer.

The original package README is in
[`IdealizedSpectralGCM.jl/README.md`](IdealizedSpectralGCM.jl/README.md).
