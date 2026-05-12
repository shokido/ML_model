# Fortran Tutorial Code

This directory contains the Fortran version of the one-dimensional mixed
layer model tutorial.

## Source Files

- `do_mlmodel_nn.f90`: main program. Edit this file to change the experiment
  period, grid spacing, initial condition, surface forcing, and output file.
- `ml_param.f90`: common parameters and constants.
- `ml_utils.f90`: utility routines, including density, stratification, shear,
  and time conversion helpers.
- `nnf_sub.f90`: turbulence closure routines used to diagnose vertical
  diffusivities.
- `solve_diag.f90`: solvers for the vertically discretized diffusion
  equations.
- `Makefile`: build recipe for the executable.

For ordinary tutorial experiments, start from `do_mlmodel_nn.f90`. The other
files provide the model machinery and usually do not need to be edited.

## Requirements

- A Fortran compiler. The provided `Makefile` uses `gfortran`.
- `make`.

If your compiler has a different command name, edit the `FC` variable in the
`Makefile`.

## Compile

From this directory:

```sh
cd /Users/skido/WORK/ML_model/TUTORIAL/Fortran
make
```

This compiles the source files in the required order and creates:

```text
exec_1d_nn.out
```

You can remove compiled objects and the executable with:

```sh
make clean
```

## Run

After compilation:

```sh
./exec_1d_nn.out
```

The default run writes:

```text
out_case_nn.txt
```

Check that the output was created:

```sh
ls out_case_nn.txt
```

## Basic Experiment Settings

The main settings are near the top of `do_mlmodel_nn.f90`:

- `dt`: time step in seconds. The default tutorial value is `240.0`.
- `dt_start`, `dt_end`: start and end date of the integration.
- `dt_output`: output interval. The default is one day.
- `bottom_depth`: total water-column depth, default `1000.0`.
- `dz`: vertical grid spacing, default `5`.
- `lat`: latitude used to compute the Coriolis parameter, default `45.0`.
- `fname_out`: output file name, default `out_case_nn.txt`.

The initial temperature, salinity, and velocity profiles are set in the
`Set initial condition` block.

The surface forcing is set inside the main time-stepping loop. The active
default is Case 1:

```fortran
hflx_nosolar=0.0_idx;hflx_solar=0.0_idx
sflx=0.0_idx
uflx=0.2;vflx=0.0
```

Comment or uncomment the nearby case blocks to try surface cooling, surface
heating, combined wind and cooling, or freshwater forcing.

## Notes on Signs and Units

- `uflx`, `vflx`: surface momentum fluxes in `N m-2`.
- `hflx_nosolar`: non-solar surface heat flux in `W m-2`. In the tutorial
  cases, negative values represent heat loss from the ocean.
- `hflx_solar`: shortwave heat flux in `W m-2`.
- `sflx`: surface salinity flux in `psu m s-1`; negative values can be used
  for freshwater input in the sample cases.

