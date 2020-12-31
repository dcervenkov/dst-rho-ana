#!/bin/bash

for NUM in $(seq -f "%.5f" 0.00000 0.00002 0.00030); do
    sed 's/"bkg_thetab_corr_f":/"bkg_thetab_corr_f": '"$NUM"',/' DSRhoCPFit/configs/nominal_rbin/config_data_sidebkg.json > config.json

    ./DSRhoCPFit/DSRhoCPFit --MC=0 --config=config.json --components=all --output=DSRhoCPFit/results/test_Kpi_$NUM --cpus=4 --exclude-channels=Kpipi0,K3pi --time-independent --plot-dir=DSRhoCPFit/plots/test/Kpi_$NUM &&\
    ./DSRhoCPFit/tools/get_chi2.py DSRhoCPFit/results/test_Kpi_${NUM}.root Kpi | tail -n1 > test_Kpi_${NUM}.chi2

    ./DSRhoCPFit/DSRhoCPFit --MC=0 --config=config.json --components=all --output=DSRhoCPFit/results/test_Kpipi0_$NUM --cpus=4 --exclude-channels=Kpi,K3pi --time-independent --plot-dir=DSRhoCPFit/plots/test/Kpipi0_$NUM &&\
    ./DSRhoCPFit/tools/get_chi2.py DSRhoCPFit/results/test_Kpipi0_${NUM}.root Kpipi0 | tail -n1 > test_Kpipi0_${NUM}.chi2

    ./DSRhoCPFit/DSRhoCPFit --MC=0 --config=config.json --components=all --output=DSRhoCPFit/results/test_K3pi_$NUM --cpus=4 --exclude-channels=Kpi,Kpipi0 --time-independent --plot-dir=DSRhoCPFit/plots/test/K3pi_$NUM &&\
    ./DSRhoCPFit/tools/get_chi2.py DSRhoCPFit/results/test_K3pi_${NUM}.root K3pi | tail -n1 > test_K3pi_${NUM}.chi2
done
