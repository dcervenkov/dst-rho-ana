#!/bin/bash

# for CHANNEL in Kpi Kpipi0 K3pi; do
#     for TYPE in lifetime mixing; do
#         # for COMPONENT in CR CRSCF; do
#         #     seq -w 00 98 | parallel --will-cite "./DSRhoLifetime -c 1 --components=${COMPONENT} --config=configs/config_mc.json --channel=${CHANNEL} --${TYPE} results/${CHANNEL}/${COMPONENT}_${TYPE}_{} ../data/${CHANNEL}/signal_mc/DSRho-mdst_${CHANNEL}_basf2_{}_svd*.root"
#         # done
#         for COMPONENT in all; do
#             seq 0 5 | parallel --will-cite "./DSRhoLifetime -c 1 --components=${COMPONENT} --config=configs/config_mc.json --channel=${CHANNEL} --${TYPE} results/${CHANNEL}/${COMPONENT}_${TYPE}_{} ../data/${CHANNEL}/realistic_mc/stream{}/*.root"
#         done
#     done
# done

mkdir -p logs
mkdir -p results
mkdir -p plots

for CHANNEL in Kpi Kpipi0 K3pi; do
    ./DSRhoLifetime -c 1 --physics --components=all --config=configs/config_mc.json --channel=${CHANNEL} --lifetime --plot-dir=plots/${CHANNEL}_mc results/${CHANNEL}_mc ../data/${CHANNEL}/realistic_mc/stream0/*.root &> logs/${CHANNEL}_mc &
    ./DSRhoLifetime -c 1 --physics --components=all --config=configs/config_data.json --channel=${CHANNEL} --lifetime --plot-dir=plots/${CHANNEL}_data results/${CHANNEL}_data ../data/${CHANNEL}/*.root &> logs/${CHANNEL}_data &
done

./DSRhoLifetime -c 1 --physics --components=all --config=configs/config_mc_all.json --channel=Kpi --lifetime --plot-dir=plots/mc results/mc ../data/K*/realistic_mc/stream0/*.root &> logs/mc &
./DSRhoLifetime -c 1 --physics --components=all --config=configs/config_data_all.json --channel=Kpi --lifetime --plot-dir=plots/data results/data ../data/K*/*.root &> logs/data &

./DSRhoLifetime -c 1 --physics --components=all --config=configs/config_mc_all_sidebands.json --channel=Kpi --lifetime --plot-dir=plots/mc_sidebands results/mc_sidebands ../data/K*/realistic_mc/stream0/DSRho-*.root ../data/K*/sidebands/*.root &> logs/mc_sidebands &

for CHANNEL in Kpi Kpipi0 K3pi; do
    ./DSRhoLifetime -c 1 --physics --components=all --config=configs/config_mc_sidebands.json --channel=${CHANNEL} --lifetime --plot-dir=plots/${CHANNEL}_mc_sidebands results/${CHANNEL}_mc_sidebands ../data/${CHANNEL}/realistic_mc/stream0/DSRho-*.root ../data/${CHANNEL}/sidebands/*.root &> logs/${CHANNEL}_mc_sidebands &
done
