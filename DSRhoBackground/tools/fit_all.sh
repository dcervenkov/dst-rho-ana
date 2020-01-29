#!/bin/bash

mkdir -p results
mkdir -p plots
mkdir -p logs

for CHANNEL in "Kpi" "Kpipi0" "K3pi"; do
    ./DSRhoBackground --physics --plot-dir=plots/${CHANNEL}_mc_scf results/${CHANNEL}_mc_scf.json ../data/${CHANNEL}/signal_mc/*.root &> logs/${CHANNEL}_mc_scf &
    ./DSRhoBackground --physics --plot-dir=plots/${CHANNEL}_mc_bkg results/${CHANNEL}_mc_bkg.json ../data/${CHANNEL}/mc_wo_signal/*.root &> logs/${CHANNEL}_mc_bkg &
    ./DSRhoBackground --physics --plot-dir=plots/${CHANNEL}_data_sidebands_bkg results/${CHANNEL}_data_sidebands_bkg.json ../data/${CHANNEL}/sidebands/*.root &> logs/${CHANNEL}_data_sidebands_bkg &
done

./DSRhoBackground --physics --plot-dir=plots/mc_scf results/mc_scf.json ../data/*/signal_mc/*.root &> logs/mc_scf &
./DSRhoBackground --physics --plot-dir=plots/mc_bkg results/mc_bkg.json ../data/*/mc_wo_signal/*.root &> logs/mc_bkg &
./DSRhoBackground --physics --plot-dir=plots/data_sidebands_bkg results/data_sidebands_bkg.json ../data/*/sidebands/*.root &> logs/data_sidebands_bkg &
