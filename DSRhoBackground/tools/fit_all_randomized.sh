#!/bin/bash

CMD_FILE="tools/cmds.tmp"

cd $(dirname $0)
cd ..

rm -f ${CMD_FILE}
touch ${CMD_FILE}
mkdir -p results
mkdir -p plots
mkdir -p logs

NUM=100

echo "./DSRhoBackground --randomize=$NUM --notail --nodelta --physics results/randomized/together_mc_scf.json         ../data/K*/signal_mc/*_0[0-5]_*.root &> logs/together_mc_scf_randomized" >> ${CMD_FILE}
echo "./DSRhoBackground --randomize=$NUM --notail           --physics results/randomized/together_data_sidebands.json ../data/K*/sidebands/*.root          &> logs/together_data_sidebands_randomized" >> ${CMD_FILE}

# Non-physics-based fits
for CHANNEL in "Kpi" "Kpipi0" "K3pi"; do
    echo "./DSRhoBackground --randomize=$NUM --angular results/randomized/nonphys_${CHANNEL}_mc_bkg.json ../data/${CHANNEL}/mc_wo_signal/*.root &> logs/nonphys_${CHANNEL}_mc_bkg" >> ${CMD_FILE}
done

echo "Running $(wc -l ${CMD_FILE} | cut -d' ' -f1) fits..."
parallel --will-cite --bar --nice 10 < ${CMD_FILE}

if [ "$?" -eq 0 ]; then
    rm ${CMD_FILE}
    exit 0
else
    exit 1
fi
