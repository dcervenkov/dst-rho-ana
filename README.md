# dst-rho-ana

[![Build Status](https://travis-ci.org/dcervenkov/dst-rho-ana.svg?branch=master)](https://travis-ci.org/dcervenkov/dst-rho-ana)

This is a collection of software tools used to perform a time-dependent angular CP violation study of the B → D* + ρ decay. 

### DSRho
*Disclaimer: A lot of the code for this module was stolen from various Belle members and used as it was. Because of this I have given up on any code and documentation standards, sorry.*

A BASF module that goes through the (hopefully skimmed) data and creates an ntuple with many variables to be used for further analysis. It does:
 - particle reconstruction
 - various cuts
 - computation of transversity parameters
 - mass and/or vertex constrained fits for D* and D
 - MC truth matching
 - continuum suppression
 - vertexing
 - tagging

### DSRhoEfficiency
The goal of this program is to test various efficiency models by fitting them to signal MC, pick a satisfactory one and then use the model together with its parameters gained from the fit in further steps of the analysis.

### DSRhoCPFit
The complete time-dependent angular fitter. Incorporates angular efficiency and Tatami. It can fit:
 - transversity amplitudes 
 - x and y Cartesian parameters 
 - lifetime
 - Δm (mixing parameter)

### DSRhoFit
Time-independent angular fitter that can be used to extract the transversity amplitudes. It incorporates efficiency parametrization.

### DSRhoLifetime
Time-dependent non-angular fitter that can be used to extract B lifetime. It incorporates Tatami.

### DSRhoPeek
This small tool reads in the ROOT file created by the DSRho module for (usually generic) MC and creates plots of all the variables with overlaid components (signal, good D & bad ρ, etc.). Several python scripts that can be used to generate publication quality plots are also included. 

### DSRhoSkim
A BASF module that goes through the whole Belle dataset and preselects data for the main analysis module. It does reconstruction and a few very simple cuts.

### DSRhoYield
This program performs a 1D Mbc fit to extract signal yield from data. It first fixes signal, self-cross-feed and background shapes and some ratios from generic MC. We introduce a smearing parameter when fitting data to account for MC vs. data differences in Mbc resolution
