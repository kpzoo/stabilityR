# stabilityR


We present Matlab code to reproduce all analyses and figures from the paper "When is the R = 1 epidemic threshold meaningful?" by Kris V Parag, Anne Cori and Uri Obolski. 

This work focuses on showing that R=1 is not a good measure of epidemic stability in (realistic) heterogeneous setting because it can hide resurgent (local) groups. We propose a previously developed statistic E (the risk averse reproduction number) as an improvement for guaranteeing local stability given a global statistic.

System Requirements

Should work with Matlab v2023a and above. Tested on macOS v13.6.9.
The implementation in R has been tested on R version 3.6.2
Dependence on the linespecer and boundedline package (license included in main folder).
Reproduction number estimates use local versions of the EpiFilter package: https://github.com/kpzoo/EpiFilter.


Instructions and installation

Should work with any standard Matlab installation and generate figures.
Run FigX.m where X is the figure in the manuscript to be reproduced.
Code is self contained so no external installations required.
Run times of all scripts are of the order of minutes or faster.
