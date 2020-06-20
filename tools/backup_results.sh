#!/bin/bash
DATE=$(date +%y%m%d)
BAKDIR="old/$DATE"

for APP in DSRhoBackground DSRhoCPFit DSRhoLifetime DSRhoYield; do
    echo "Backing up $APP logs, plots, and results to $BAKDIR"
    cd ../$APP/ && mkdir $BAKDIR && mv logs plots results $BAKDIR/. && mkdir logs plots results && cd ../tools
done
