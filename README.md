# PythonMEEP
Script to simulate a resonant chamber simulated in python MEEP

## Installation
The software and packages needed are specified in [installation](https://meep.readthedocs.io/en/latest/Installation/).

## What the script does
The script resonant_chamber.py performs the FDTD electromagnetic simulation in MEEP of a rectangular resonant chamber. The simulation outputs files with dielectric and electric field (Ez) information in selected places.
Finally it plots the permittivity maps, Ez fields and performs the FFT of the Ez signal in the receiver antennas.
