#!/bin/bash

for CHANNEL in Kpi Kpipi0 K3pi; do
    sed "s/Kpi/${CHANNEL}/g" plots.json > plots_${CHANNEL}.json
    ./publication_plot.py --format=.pdf --dir=plots_${CHANNEL} plots_${CHANNEL}.json
done
