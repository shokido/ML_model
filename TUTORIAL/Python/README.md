# Python Tutorial Code

This directory contains the Python version of the one-dimensional mixed layer
model tutorial. It mirrors the Fortran version closely, but does not require
compilation.

## Source Files

- `do_mlmodel_nn.py`: main script. Edit this file to change the experiment
  period, grid spacing, initial condition, surface forcing, and output file.
- `ml_param.py`: common parameters and constants.
- `ml_util.py`: utility functions, including density, stratification, shear,
  and shortwave absorption helpers.
- `nnf_sub.py`: turbulence closure functions used to diagnose vertical
  diffusivities.
- `solve_diag.py`: solvers for the vertically discretized diffusion
  equations.

For ordinary tutorial experiments, start from `do_mlmodel_nn.py`. The other
files are imported automatically.

## Requirements

Use Python 3 with:

- `numpy`
- `scipy`

The standard-library modules `datetime` and `time` are also used.

Install the external packages if needed:

```sh
python -m pip install numpy scipy
```

Depending on your environment, the Python command may be `python3` instead of
`python`.

## Run

From this directory:

```sh
cd /Users/skido/WORK/ML_model/TUTORIAL/Python
python do_mlmodel_nn.py
```

The default run writes:

```text
out_case_nn.txt
```

Check that the output was created:

```sh
ls out_case_nn.txt
```

No compilation step is required. The helper files are loaded through Python
imports.

## Basic Experiment Settings

The main settings are near the top of `do_mlmodel_nn.py`:

- `fname_out`: output file name, default `out_case_nn.txt`.
- `dt_start`, `dt_end`: start and end date of the integration.
- `dz`: vertical grid spacing in meters, default `5.0`.
- `lat`: latitude used to compute the Coriolis parameter, default `45.0`.
- `dt`: time step in seconds. The default tutorial value is `240.0`.
- `dt_output`: output interval. The default is one day.
- `bottom_depth`: total water-column depth, default `1000.0`.

The initial temperature, salinity, and velocity profiles are set in the block
beginning with:

```python
# Set initial conditions for temperature, salinity, and velocity
```

The surface forcing is set inside the main time-stepping loop. The active
default is Case 1:

```python
hflx_nosolar = 0.0; hflx_solar = 0.0
sflx = 0.0
uflx = 0.2; vflx = 0.0
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

## Plotting

To plot the output with the sample script:

```sh
cd /Users/skido/WORK/ML_model/TUTORIAL/Gallery
cp ../Python/out_case_nn.txt .
python plot_results.py
```

The script draws time-depth sections of temperature, salinity, zonal velocity,
and meridional velocity.

