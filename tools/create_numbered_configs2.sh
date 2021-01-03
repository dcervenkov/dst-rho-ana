#!/bin/bash

cd tools
./create_randomized_configs.py ../DSRhoCPFit/configs/randomized_yield_rbin_corr/templates/config_data_mcbkg.template.json -d ../DSRhoCPFit/configs/randomized_yield_rbin_corr/numbered -r Yield/results/,Yield/results_randomized_rbin_corr/rnd_{}/ -n100
cd -
