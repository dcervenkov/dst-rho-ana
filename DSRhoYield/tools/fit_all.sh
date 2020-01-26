#!/bin/bash

mkdir -p logs
mkdir -p plots
for CHANNEL in "Kpi" "Kpipi0" "K3pi"; do
    for STREAM in $(seq 0 5); do
        ./DSRhoYield --MC ../data/${CHANNEL}/realistic_mc ../data/${CHANNEL}/realistic_mc/stream${STREAM} plots/${CHANNEL}_stream${STREAM} &> logs/${CHANNEL}_stream${STREAM} &
    done
    ./DSRhoYield ../data/${CHANNEL}/realistic_mc ../data/${CHANNEL} plots/${CHANNEL} &> logs/${CHANNEL} &
done
