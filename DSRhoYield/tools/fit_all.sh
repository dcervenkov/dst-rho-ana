#!/bin/bash

CMD_FILE="tools/cmds.tmp"

cd $(dirname $0)
cd ..

rm -f ${CMD_FILE}
touch ${CMD_FILE}
mkdir -p logs
mkdir -p plots

for CHANNEL in "Kpi" "Kpipi0" "K3pi"; do
    for STREAM in $(seq 0 5); do
        echo "./DSRhoYield --MC ../data/${CHANNEL}/realistic_mc ../data/${CHANNEL}/realistic_mc/stream${STREAM} plots/${CHANNEL}_stream${STREAM} &> logs/${CHANNEL}_stream${STREAM}" >> ${CMD_FILE}
    done
    echo "./DSRhoYield ../data/${CHANNEL}/realistic_mc ../data/${CHANNEL} plots/${CHANNEL} &> logs/${CHANNEL}" >> ${CMD_FILE}
done

echo "Running $(wc -l ${CMD_FILE} | cut -d' ' -f1) fits..."
parallel --will-cite --bar --nice 10 < ${CMD_FILE}

if [ "$?" -eq 0 ]; then
    rm ${CMD_FILE}
    exit 0
else
    exit 1
fi
