# spaceseries
 Estimation of space series from the early Swarm mission pearls-on-a-string orbital formation.

 ## Overview

For an overview and an application to auroral physics, see Burchill (2023), "How variable are Birkeland currents?", submitted to Earth, Planets and Space.

The software is written in C mainly to make processing of the several GB of Swarm data fast.

One way to compile the source code is to install Docker and run, from a terminal, the docker-build.sh and docker-setup-dev-container.sh scripts from the Docker folder, in that order, then open the source folder in Visual Studio Code, connect to the remote container SpaceSeries, install the Microsoft C/C++ Extension Pack and CMake Tools extensions on the container, select a build kit, run a terminal in VS Code, cd into the build folder, type make, and press the "return" key.

For details on the space series calculations, see spaceseries.c. The source compiles to the binary `swarm_space_series`.

## Dependencies

- [Quasidipole magnetic coordinates library](https://github.com/JohnathanBurchill/QuasiDipole)
- An _apexsh.dat_ binary file containing unformatted spherical harmonics coefficients for the calculation of quasi-dipole magnetic latitude and local times is included in this distribution. This file is different than the one generated by the vanilla Emmert and Richmond Fortran code; it does not contain superfluous Fortran binary formatting. Place that file in the working directory of the `swarm_space_series` program, or provide its location as an option to the program.

- [NASA CDF C library](https://cdf.gsfc.nasa.gov/html/sw_and_docs.html)
- [GNU Scientific Library](https://www.gnu.org/software/gsl/)

## Inputs

Processing requires data files for two or three satellites. Only one file should be present for each satellite. The behaviour is not well defined when multiple files for a given day, with different time ranges, are present.

When processing ion drift from the TII cross-track ion drift dataset (version 0302), no additional files are needed.

If processing Langmuir probe parameters or magnetic field residuals (or the vertical component of electric current density), ephemeres will be loaded from TII cross-track ion drift files, if present, otherwise Swarm Medium Orbit Determination (MOD) files must be present covering the requested time range.

Magnetic residual files are produced by the [`chaos`](https://github.com/JohnathanBurchill/chaos) program.

## Output

The software stores space series either as unformatted binary (default) or formatted ASCII output files. Each series consists of a header with the timestamp, number of space points, and satellite ephemeres at that timestamp, followed by the space series data points. Each point includes ephemeres, the derived measurement, its time derivative, and for each satellite the difference between the satellite's measurement at that location and the derived measurement (i.e., the error at the given time associated with assuming a steady-state environment). 

## Analysis

For statistical analysis of space series generated with `swarm_space_series` based on quasi-dipole latitude and magnetic local time, have a look at _analyze_space_series.c_ and _statistics.c_, which compile to `analyze_space_series`.


