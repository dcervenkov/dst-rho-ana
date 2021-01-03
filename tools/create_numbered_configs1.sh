#!/bin/bash

# Copy `DSRhoCPFit/efficiencies/randomized` (`DSRhoEfficiency` itself can generate these randomized efficiency models)
# Copy `DSRhoCPFit/scfmodels/randomized` (`DSRhoBackground` itself can generate these randomized SCF models (maybe))
# Copy `DSRhoBackground/results/randomized_corr` (`DSRhoBackground` itself can generate these randomized BKG models (maybe))

cd tools
./create_randomized_configs.py ../DSRhoCPFit/configs/randomized_eff/templates/config_data_mcbkg.template.json -d ../DSRhoCPFit/configs/randomized_eff/numbered -r 200729,200729_rnd_{} -r eff_,randomized/eff_ -n100
./create_randomized_configs.py ../DSRhoCPFit/configs/randomized_scf/templates/config_data_mcbkg.template.json -d ../DSRhoCPFit/configs/randomized_scf/numbered -r 200619,200802_rnd_{} -r scfmodels,scfmodels/randomized -n100
./create_randomized_configs.py ../DSRhoCPFit/configs/randomized_bkg_corr/templates/config_data_mcbkg.template.json -d ../DSRhoCPFit/configs/randomized_bkg_corr/numbered -r Kpi_mc_bkg.json,randomized_corr/Kpi_mc_bkg_rnd_{}.json -r Kpipi0_mc_bkg.json,randomized_corr/Kpipi0_mc_bkg_rnd_{}.json -r K3pi_mc_bkg.json,randomized_corr/K3pi_mc_bkg_rnd_{}.json -n100
./create_randomized_configs.py ../DSRhoCPFit/configs/randomized_bkg_dt_corr/templates/config_data_mcbkg.template.json -d ../DSRhoCPFit/configs/randomized_bkg_dt_corr/numbered -r together_data_sidebands.json,randomized_corr/together_data_sidebands_rnd_{}.json -n100
cd -

cd tools
./config_from_template.py ../DSRhoCPFit/configs/sigma_tau_dm/templates/config_data_mcbkg.template.json > ../DSRhoCPFit/configs/sigma_tau_dm/numbered/config_data_mcbkg_0.json
cd -
cd DSRhoCPFit/configs/sigma_tau_dm/numbered
sed 's/"tau": 1.519/"tau": 1.515/' config_data_mcbkg_0.json > config_data_mcbkg_1.json
sed 's/"tau": 1.519/"tau": 1.523/' config_data_mcbkg_0.json > config_data_mcbkg_2.json
sed 's/"dm": 0.5065/"dm": 0.5046/' config_data_mcbkg_0.json > config_data_mcbkg_3.json
sed 's/"dm": 0.5065/"dm": 0.5084/' config_data_mcbkg_0.json > config_data_mcbkg_4.json
cd -

cd tools
./config_from_template.py ../DSRhoCPFit/configs/measured_scf/templates/config_data_mcbkg.template.json > ../DSRhoCPFit/configs/measured_scf/numbered/config_data_mcbkg_0.json
cd -
cd DSRhoCPFit/configs/measured_scf/numbered
sed 's/scf_Kpipi0_200619.root/scf_Kpipi0_200816_measured.root/' config_data_mcbkg_0.json > config_data_mcbkg_1.json
cd -

cd tools
./config_from_template.py ../DSRhoCPFit/configs/scf_swap/templates/config_data_mcbkg.template.json > ../DSRhoCPFit/configs/scf_swap/numbered/config_data_mcbkg_0.json
cd -
