# Tutorial: 1-D Ocean Mixed Layer Model

This directory contains a compact tutorial version of a one-dimensional
ocean mixed layer model written in both Fortran and Python. The model is
designed for small experiments that can be run on a desktop computer, while
retaining the essential physics of vertical mixing in the upper ocean.

The tutorial code is based on the mixed layer model framework of Mellor and
Yamada (1982), Furuichi et al. (2012), and Kido and Tozuka (2017). It solves
the time evolution of temperature, salinity, and horizontal velocity in a
single water column. Horizontal advection can be prescribed as an external
tendency, but the central process represented explicitly here is vertical
turbulent transport.

## Directory Layout

- `Fortran/`: Fortran source code and a `Makefile`.
- `Python/`: Python source code with the same basic experiment setup.
- `Gallery/`: sample output, plotting script, and an example figure.

The Fortran and Python versions are intentionally similar. In most tutorial
experiments, you only need to edit the main driver:

- `Fortran/do_mlmodel_nn.f90`
- `Python/do_mlmodel_nn.py`

## Physical Problem

The tutorial experiments consider a water column from the sea surface to
1000 m depth at 45 degrees north. The default vertical grid spacing is 5 m,
so tracer and velocity values are defined at layer centers from 2.5 m to
997.5 m.

The default experiment integrates from 1900-01-01 00:00:00 to
1900-01-21 00:00:00 with a 240 s time step. Output is written once per day.

The default initial condition is:

- temperature is 15 degC from the surface to 20 m depth;
- temperature decreases linearly from 15 degC at 20 m to 12 degC at 50 m;
- temperature is 12 degC below 50 m;
- salinity is vertically uniform at 35 psu;
- zonal and meridional velocities are initially zero.

The default forcing is the tutorial wind-mixing case:

- zonal surface momentum flux: `0.2 N m-2`;
- meridional surface momentum flux: `0.0 N m-2`;
- surface heat flux: `0.0 W m-2`;
- surface salinity flux: `0.0 psu m s-1`.

This represents a westerly wind stress that mixes warm surface water with
colder water below.

## Governing Equations

The model solves simplified one-dimensional forms of the ocean model
equations. Only the vertical turbulent transport terms are solved explicitly
inside the column, while other processes may be supplied as external
tendencies.

For potential temperature `T`,
```math
\frac{\partial T}{\partial t} = \frac{\partial}{\partial z}( \kappa_{T}\frac{\partial T}{\partial z}  ) + g_T(z) + \frac{1}{\rho_{0} C_p} \frac{\partial I(z)}{\partial z}
```
```math
\frac{\partial S}{\partial t} = \frac{\partial}{\partial z}( \kappa_{S}\frac{\partial S}{\partial z}  ) + g_S(z)
```
For zonal and meridional velocity `u` and `v`,
```math
\frac{\partial u}{\partial t} = \frac{\partial}{\partial z}( \kappa_{M}\frac{\partial u}{\partial z})+fv + g_u(z)
```
```math
\frac{\partial v}{\partial t} = \frac{\partial}{\partial z}( \kappa_{M}\frac{\partial v}{\partial z})-fu + g_u(z)
```
where:

- $z$ is the vertical coordinate;
- $\kappa_{M}$ is vertical eddy viscosity for momentum;
- $\kappa_{T}$ and $\kappa_{S}$ are vertical eddy diffusivities for temperature and salinity;
- $g_T$, $g_S$, $g_u$, and $g_v$ are prescribed external tendencies, such as
  advection;
- $I(z)$ is penetrating shortwave radiation;
- $rho_{0}$ is reference seawater density;
- $C_{p}$ is seawater heat capacity;
- $f$ is the Coriolis parameter.

The turbulent diffusivities are diagnosed from turbulent quantities:

```math
\kappa_{M} = S_{M}ql
\kappa_{T}=\kappa_{S} = S_{H}ql
```

Here $q$ is related to turbulent kinetic energy, $l$ is the turbulent length
scale, and $S_{M}$ and $S_{H}$ are stability functions. The strength of vertical mixing depends mainly on stratification, vertical shear, and surface forcing.

## Changing Experiments

The easiest way to modify the tutorial case is to edit the forcing block in
the main driver. Several cases are already included as commented examples:

- Case 1: wind mixing with `uflx = 0.2`, no heat or salinity flux.
- Case 2: surface cooling with `hflx_nosolar = -30.0`, no wind stress.
- Case 3: surface heating with `hflx_nosolar = 30.0`, no wind stress.
- Case 4: wind mixing plus surface cooling.
- Case 5: wind mixing plus freshwater input through negative salinity flux.

Useful parameters near the top of the main driver include:

- `dt`: time step in seconds;
- `dt_start`, `dt_end`: integration period;
- `dt_output`: output interval;
- `bottom_depth`: bottom depth in meters;
- `dz`: vertical grid spacing in meters;
- `lat`: latitude used for the Coriolis parameter;
- `fname_out`: output file name.

## Output Format

Both versions write results to `out_case_nn.txt` by default.

The file structure is:

1. number of output time steps;
2. output depths in meters;
3. date and time for the first output record;
4. temperature profile at that time;
5. salinity profile at that time;
6. zonal velocity profile at that time;
7. meridional velocity profile at that time;
8. the next output date and time, followed by the same four profiles;
9. repeated until the end of the simulation.

The profile values are instantaneous values at each output time, not temporal
averages.

## Plotting Results

After generating `out_case_nn.txt`, you can visualize it with the sample
script in `Gallery/`:

```sh
cd /Users/skido/WORK/ML_model/TUTORIAL/Gallery
cp ../Python/out_case_nn.txt .
python plot_results.py
```

If you ran the Fortran version instead, copy `../Fortran/out_case_nn.txt`.
The plotting script reads the text output and draws time-depth sections of
temperature, salinity, zonal velocity, and meridional velocity.