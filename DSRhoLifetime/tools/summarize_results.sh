#!/bin/bash

echo "Lifetime"
for CHANNEL in Kpi Kpipi0 K3pi together; do
    echo "${CHANNEL}"
    for TYPE in CR CRSCF all; do
        printf "%6s " ${TYPE}
        awk '{life += $1; err += $2; count++} END {printf "%.3f +- %.3f\n",life/count,err/count}' ../results/${CHANNEL}/${TYPE}_lifetime_*[0-6]
    done
    echo
done

echo "Mixing"
for CHANNEL in Kpi Kpipi0 K3pi together; do
    echo "${CHANNEL}"
    for TYPE in CR CRSCF all; do
        printf "%6s " ${TYPE}
        awk '{tau += $1; tau_err += $2; dm += $3; dm_err += $4; count++} END {printf "%.3f +- %.3f, %.3f +- %.3f\n",tau/count, tau_err/count, dm/count, dm_err/count}' ../results/${CHANNEL}/${TYPE}_mixing_*[0-6]
    done
    echo
done
