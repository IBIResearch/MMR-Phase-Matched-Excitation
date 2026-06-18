#  Phase-Matched Chirps for Energy-Efficient Excitation of Magneto-Mechanical Resonators
This repository contains code to simulate different excitation waveforms for magneto-mechanical resonators. The different waveforms are sinusoidal and rectangular with different instantaneous phase models -- fixed frequency and linear chirp-

The method corresponding to this code is described in the associated publication (see below).

## Installation
In order to use this code, one first has to download [Julia](https://julialang.org/) (version 1.11 or later) and clone this repository.

## Execution
After installation the example code can be executed by navigating to the folder, running `julia` and entering
```
include("src/main.jl")
```
to simulate experiments with excitation of all proposed excitation signals. The example script automatically activates the environment and installs all necessary packages. This will take several minutes before the actual code is run, since all packages are precompiled during installation. Resulting images will be saved in the folder results.

## Citation
If you use this code in your research, please cite the following paper:
TODO
