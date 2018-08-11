#!/usr/bin/env python

"""Run fits with all possible combinations of the fitted parameters.

The results are saved to files named par1_par2...root
"""

from itertools import combinations
from subprocess import run

pars = ["ap", "apa", "a0", "ata"]

par_combinations = []
for i in range(4):
    par_combinations += combinations(pars, i)

for par_combination in par_combinations:
    print(par_combination)
    run(["./DSRhoCPFit", "--efficiency-model=0", "--time-independent",
         "--perfect-tag", "--cpus=4", "--fit=CR",
         "--fix=" + ",".join(par_combination),
         "../data/evtgen/basf2_mod_real_unmod_1.root",
         "results/" + "_".join(par_combination)],
        cwd="..")
