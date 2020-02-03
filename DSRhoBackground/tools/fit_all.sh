#!/bin/bash

CMD_FILE="tools/cmds.tmp"

cd $(dirname $0)
cd ..

rm -f ${CMD_FILE}
touch ${CMD_FILE}
mkdir -p results
mkdir -p plots
mkdir -p logs

for CHANNEL in "Kpi" "Kpipi0" "K3pi"; do
    echo "./DSRhoBackground --physics --plot-dir=plots/${CHANNEL}_mc_scf results/${CHANNEL}_mc_scf.json ../data/${CHANNEL}/signal_mc/*.root &> logs/${CHANNEL}_mc_scf" >> ${CMD_FILE}
    echo "./DSRhoBackground --physics --plot-dir=plots/${CHANNEL}_mc_bkg results/${CHANNEL}_mc_bkg.json ../data/${CHANNEL}/mc_wo_signal/*.root &> logs/${CHANNEL}_mc_bkg" >> ${CMD_FILE}
    echo "./DSRhoBackground --physics --plot-dir=plots/${CHANNEL}_data_sidebands_bkg results/${CHANNEL}_data_sidebands_bkg.json ../data/${CHANNEL}/sidebands/*.root &> logs/${CHANNEL}_data_sidebands_bkg" >> ${CMD_FILE}
done

echo "./DSRhoBackground --physics --plot-dir=plots/mc_scf results/mc_scf.json ../data/*/signal_mc/*.root &> logs/mc_scf" >> ${CMD_FILE}
echo "./DSRhoBackground --physics --plot-dir=plots/mc_bkg results/mc_bkg.json ../data/*/mc_wo_signal/*.root &> logs/mc_bkg" >> ${CMD_FILE}
echo "./DSRhoBackground --physics --plot-dir=plots/data_sidebands_bkg results/data_sidebands_bkg.json ../data/*/sidebands/*.root &> logs/data_sidebands_bkg" >> ${CMD_FILE}

echo "Running $(wc -l ${CMD_FILE} | cut -d' ' -f1) fits..."
parallel --will-cite --bar --nice 10 < ${CMD_FILE}

if [ "$?" -eq 0 ]; then
    rm ${CMD_FILE}
    exit 0
else
    exit 1
fi
