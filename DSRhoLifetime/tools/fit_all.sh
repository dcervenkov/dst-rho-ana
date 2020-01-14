#!/bin/bash

for CHANNEL in Kpi Kpipi0 K3pi; do
    for TYPE in lifetime mixing; do
        # for COMPONENT in CR CRSCF; do
        #     seq -w 00 98 | parallel --will-cite "./DSRhoLifetime -c 1 --components=${COMPONENT} --config=config.json --channel=${CHANNEL} --${TYPE} results/${CHANNEL}/${COMPONENT}_${TYPE}_{} ../data/${CHANNEL}/signal_mc/DSRho-mdst_${CHANNEL}_basf2_{}_svd*.root"
        # done
        for COMPONENT in all; do
            seq 0 5 | parallel --will-cite "./DSRhoLifetime -c 1 --components=${COMPONENT} --config=config.json --channel=${CHANNEL} --${TYPE} results/${CHANNEL}/${COMPONENT}_${TYPE}_{} ../data/${CHANNEL}/realistic_mc/stream{}/*.root"
        done
    done
done
