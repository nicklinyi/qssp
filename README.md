# QSSP - Complete synthetic seismograms for a spherical earth

Notes: This is a modified version of the original code. 

## What I have done?
- add a Makefile to compile the code
- modify qpmain.f to qpmain.f90 so we can pass the command-line arguments to the code

Code for calculating complete synthetic seismograms of a spherical earth using the normal mode theory.

## Highlights
(1) all-in-one code for body waves, surface waves, free oscillations, tsunami for uniform ocean, infrasound waves for a standard atmosphere and static deformation
(2) generating Green’s function database or simulating complete seismograms for any given kinematic source model
(3) hybrid algorithm (numerical integration for low frequency / small harmonic degrees and analytical propagator algorithm for high frequency / large harmonic degrees)
(4) complex frequency technique for supressing the time-domain aliasing problem
(5) differential filter technique for suppressing numerical phases

## Related codes
QSSPSTATIC - Co- and post-seismic viscoelastic deformation based on a spherical visco-elastic-gravitational earth model.

## Downloads
- ftp://ftp.gfz-potsdam.de/pub/home/turk/wang/qssp2017-code+input.rar

- ftp://ftp.gfz-potsdam.de/pub/home/turk/wang/qsspstatic-code+input.rar

## References
- Wang, R., S. Heimann, Y. Zhang, H. Wang, and T. Dahm (2017). Complete synthetic seismograms based on a spherical self-gravitating Earth model with an atmos-phere-ocean-mantle-core structure. Geophysical Journal International, doi: 10.1093/gji/ggx259.