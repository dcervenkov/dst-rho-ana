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

### DSRhoBackground
Model the background angular and temporal distributions. The angular result is a 1D x 1D x 1D functional form or a RooHistPdf when correlations are important. Delta t can be either an empirical PDF or a per-event dependent physics-inspired Tatami-based PDF.

### DSRhoEfficiency
The goal of this program is to test various efficiency models by fitting them to signal MC, pick a satisfactory one and then use the model together with its parameters gained from the fit in further steps of the analysis.

### DSRhoCPFit
The complete fitter time-dependent angular fitter. Incorporates angular efficiency and Tatami. It can fit:
 - transversity amplitudes
 - x and y Cartesian parameters
 - lifetime
 - Δm (mixing parameter)
It can also be switched to use a time-independent PDF which gives access to the transversity amplitudes only, but is much faster and simpler then the time-dependent PDF.

### DSRhoLifetime
Time-dependent non-angular fitter that can be used to extract B lifetime. It incorporates Tatami.

### DSRhoPeek
This small tool reads in the ROOT file created by the DSRho module for (usually generic) MC and creates plots of all the variables with overlaid components (signal, good D & bad ρ, etc.). Several python scripts that can be used to generate publication quality plots are also included.

### DSRhoSkim
A BASF module that goes through the whole Belle dataset and preselects data for the main analysis module. It does reconstruction and a few very simple cuts.

### DSRhoSidebands
A program to model and fit the shape difference of background in the signal window and sidebands. It's output is part of the DSRhoCPFit config and improves the correspondence between our model and data.

### DSRhoTMVATraining
Adapted TMVA tutorial code used to train continuum suppression classifiers.

### DSRhoYield
This program performs a 1D Mbc fit to extract signal yield from data. It first fixes signal, self-cross-feed and background shapes and some ratios from generic MC. We introduce a smearing parameter when fitting data to account for MC vs. data differences in Mbc resolution

### panther2TTree
A tool to convert `hbook` files to `ROOT` files. It's a one to one conversion, ntuple field names are preserved.

### printTree
Allows you to print the decay tree of chosen event(s) in a human-friendly format.
