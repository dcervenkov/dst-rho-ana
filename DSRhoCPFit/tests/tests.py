"""A module that runs fit tests (end-to-end tests).

This module runs predefined fit tests and compares their results to reference
ones. It outputs a summary table and quits with the exit code 0 when all tests
pass and 1 when any of them fail.
"""
import os
import subprocess
import sys

# Define tests to be run
test_configs = {
    "td_cr_fit": ["--efficiency-model=6",
                  "--efficiency-file=eff_Kpi_190531.root",
                  "--fit=CR",
                  "--mixing",
                  "--fix=apa,a0,ata,x0,xt,yp,y0,yt,xpb,x0b,xtb,ypb,y0b,ytb",
                  "tests/current_result",
                  "tests/data/signalMC.root"],

    "ti_cr_fit": ["--efficiency-model=6",
                  "--efficiency-file=eff_Kpi_190531.root",
                  "--fit=CR",
                  "--time-independent",
                  "tests/current_result",
                  "tests/data/signalMC.root"],

    "td_crscf_fit": ["--efficiency-model=6",
                     "--efficiency-file=eff_Kpi_190531.root",
                     "--fit=CRSCF",
                     "--scf-histo=scf_Kpi_190531.root",
                     "--mixing",
                     "--fix=apa,a0,ata,x0,xt,yp,y0,yt,xpb,x0b,xtb,ypb,y0b,ytb",
                     "tests/current_result",
                     "tests/data/signalMC.root"],

    "ti_crscf_fit": ["--efficiency-model=6",
                     "--efficiency-file=eff_Kpi_190531.root",
                     "--fit=CRSCF",
                     "--scf-histo=scf_Kpi_190531.root",
                     "--time-independent",
                     "tests/current_result",
                     "tests/data/signalMC.root"],

    "td_all_fit": ["--efficiency-model=6",
                   "--efficiency-file=eff_Kpi_190531.root",
                   "--config=config_Kpi.json",
                   "--fit=all",
                   "--scf-histo=scf_Kpi_190531.root",
                   "--mixing",
                   "--fix=apa,a0,ata,x0,xt,yp,y0,yt,xpb,x0b,xtb,ypb,y0b,ytb",
                   "tests/current_result",
                   "tests/data/signalMC.root",
                   "tests/data/chargedMC.root",
                   "tests/data/charmMC.root",
                   "tests/data/mixedMC.root",
                   "tests/data/udsMC.root"],

    "ti_all_fit": ["--efficiency-model=6",
                   "--efficiency-file=eff_Kpi_190531.root",
                   "--config=config_Kpi.json",
                   "--fit=all",
                   "--scf-histo=scf_Kpi_190531.root",
                   "--time-independent",
                   "tests/current_result",
                   "tests/data/signalMC.root",
                   "tests/data/chargedMC.root",
                   "tests/data/charmMC.root",
                   "tests/data/mixedMC.root",
                   "tests/data/udsMC.root"],
}


def run_test(name, config):
    """Run a single test defined by its name and config."""
    print("Running " + name)
    config.insert(0, "./DSRhoCPFit")
    return_code = subprocess.call(config, cwd="../.")

    reference_result = ""
    current_result = ""

    with open("references/" + name + ".reference", "r") as f:
        reference_result = f.readline()

    with open("current_result", "r") as f:
        current_result = f.readline()

    if return_code == 0 and reference_result == current_result:
        print("The result matches the reference result.")
        os.remove("current_result")
        os.remove("current_result.root")
        return 0, name
    else:
        print("The result does NOT match the reference result!")
        os.rename("current_result", name + ".result")
        os.rename("current_result.root", name + ".result.root")
        return 1, name


def print_results_table(results):
    """Print a summary of the results in a table."""
    print("\nSummary")
    print("-" * 27)

    for result in test_results:
        result_str = ""
        if result[0]:
            result_str = "FAIL"
        else:
            result_str = "pass"

        print("{:<20} | {}".format(result[1], result_str))

    num_failed = sum(x[0] for x in test_results)

    print("-" * 27)
    print("{:<20} | {:>4}".format("Total Failed", num_failed))


test_results = []

for test_name, test_config in test_configs.items():
    # Skip tests that weren't explicitely requested if some were
    if len(sys.argv) > 1 and test_name not in sys.argv[1:]:
        continue

    test_results.append(run_test(test_name, test_config))

print_results_table(test_results)


# Return code for Travis automation
if any([result[0] for result in test_results]):
    sys.exit(1)
else:
    sys.exit(0)
