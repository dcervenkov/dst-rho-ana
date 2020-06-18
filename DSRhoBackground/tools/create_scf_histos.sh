#!/bin/bash

DATE=$(date +%y%m%d)

cd ..
./DSRhoBackground --histo --plot-dir=plots/scf results/scf_Kpi_${DATE}.root    ../data/basf2_190307_Kpi/*
./DSRhoBackground --histo --plot-dir=plots/scf results/scf_Kpipi0_${DATE}.root ../data/basf2_190528_Kpipi0/*
./DSRhoBackground --histo --plot-dir=plots/scf results/scf_K3pi_${DATE}.root   ../data/basf2_190529_K3pi/*
