#!/bin/bash

CMD_FILE="tools/cmds.tmp"

cd $(dirname $0)
cd ..

rm -f ${CMD_FILE}
touch ${CMD_FILE}
mkdir -p logs
mkdir -p results
mkdir -p plots

# Plotless fits of all types on all available MC (99 streams for CR and CRSCF)
for CHANNEL in Kpi Kpipi0 K3pi; do
    mkdir -p results/${CHANNEL}
    mkdir -p logs/${CHANNEL}
    for TYPE in lifetime mixing; do
        for COMPONENT in CR CRSCF; do
            for STREAM in $(seq -w 00 98); do
                echo "./DSRhoLifetime -c 1 --physics --components=${COMPONENT} --config=configs/config_mc.json --channel=${CHANNEL} --${TYPE} results/${CHANNEL}/${COMPONENT}_${TYPE}_${STREAM} ../data/${CHANNEL}/signal_mc/DSRho-mdst_${CHANNEL}_basf2_${STREAM}_svd?.root &> logs/${CHANNEL}/${COMPONENT}_${TYPE}_${STREAM}" >> ${CMD_FILE}
            done
        done
        for COMPONENT in all; do
            for STREAM in $(seq -w 0 5); do
                echo "./DSRhoLifetime -c 1 --physics --components=${COMPONENT} --config=configs/config_mc.json --channel=${CHANNEL} --${TYPE} results/${CHANNEL}/${COMPONENT}_${TYPE}_${STREAM} ../data/${CHANNEL}/realistic_mc/stream${STREAM}/*.root &> logs/${CHANNEL}/${COMPONENT}_${TYPE}_${STREAM}" >> ${CMD_FILE}
            done
        done
    done
done

# Fits with all components, with plots, on all available MC (6 streams) and data
for CHANNEL in Kpi Kpipi0 K3pi; do
    for STREAM in $(seq 0 5); do
        echo "./DSRhoLifetime -c 1 --physics --components=all --config=configs/config_mc.json --channel=${CHANNEL} --lifetime --plot-dir=plots/${CHANNEL}_mc_stream${STREAM} results/${CHANNEL}_mc_stream${STREAM} ../data/${CHANNEL}/realistic_mc/stream${STREAM}/*.root &> logs/${CHANNEL}_mc_stream${STREAM}" >> ${CMD_FILE}
    done
    echo "./DSRhoLifetime -c 1 --physics --components=all --config=configs/config_data.json --channel=${CHANNEL} --lifetime --plot-dir=plots/${CHANNEL}_data results/${CHANNEL}_data ../data/${CHANNEL}/*.root &> logs/${CHANNEL}_data" >> ${CMD_FILE}
done

# As above, but with all channels together in a single fit
for STREAM in $(seq 0 5); do
    echo "./DSRhoLifetime -c 1 --physics --components=all --config=configs/config_mc_all.json --channel=Kpi --lifetime --plot-dir=plots/mc_stream${STREAM} results/mc_stream${STREAM} ../data/K*/realistic_mc/stream${STREAM}/*.root &> logs/mc_stream${STREAM}" >> ${CMD_FILE}
done
echo "./DSRhoLifetime -c 1 --physics --components=all --config=configs/config_data_all.json --channel=Kpi --lifetime --plot-dir=plots/data results/data ../data/K*/*.root &> logs/data" >> ${CMD_FILE}


# As above but with sidebands substituted for MC BKG
for CHANNEL in Kpi Kpipi0 K3pi; do
    echo "./DSRhoLifetime -c 1 --physics --components=all --config=configs/config_mc_sidebands.json --channel=${CHANNEL} --lifetime --plot-dir=plots/${CHANNEL}_mc_sidebands results/${CHANNEL}_mc_sidebands ../data/${CHANNEL}/realistic_mc/stream0/DSRho-*.root ../data/${CHANNEL}/sidebands/*.root &> logs/${CHANNEL}_mc_sidebands" >> ${CMD_FILE}
done
echo "./DSRhoLifetime -c 1 --physics --components=all --config=configs/config_mc_all_sidebands.json --channel=Kpi --lifetime --plot-dir=plots/mc_sidebands results/mc_sidebands ../data/K*/realistic_mc/stream0/DSRho-*.root ../data/K*/sidebands/*.root &> logs/mc_sidebands" >> ${CMD_FILE}

echo "Running $(wc -l ${CMD_FILE} | cut -d' ' -f1) fits..."
parallel --will-cite --bar < ${CMD_FILE}

if [ "$?" -eq 0 ]; then
    rm ${CMD_FILE}
    exit 0
else
    exit 1
fi
