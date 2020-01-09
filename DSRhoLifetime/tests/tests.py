"""A module that runs fit tests (end-to-end tests).

This module runs predefined fit tests and compares their results to reference
ones. It outputs a summary table and quits with the exit code 0 when all tests
pass and 1 when any of them fail.
"""
import os
import subprocess
import sys
import time

# Define command to be run
cmd = "./DSRhoLifetime"

# Define tests to be run
test_configs = {
    "lifetime_fit_CR": ["--config=tests/config.json",
                        "--channel=Kpi",
                        "--components=CR",
                        "--lifetime",
                        "tests/current_result",
                        "tests/data/Kpi.root"],
    "lifetime_fit_CRSCF": ["--config=tests/config.json",
                           "--channel=Kpi",
                           "--components=CRSCF",
                           "--lifetime",
                           "tests/current_result",
                           "tests/data/Kpi.root"],
    "lifetime_fit_all": ["--config=tests/config.json",
                         "--channel=Kpi",
                         "--components=all",
                         "--lifetime",
                         "tests/current_result",
                         "tests/data/Kpi.root"],
    "mixing_fit_CR": ["--config=tests/config.json",
                      "--channel=Kpi",
                      "--components=CR",
                      "--mixing",
                      "tests/current_result",
                      "tests/data/Kpi.root"],
    "mixing_fit_CRSCF": ["--config=tests/config.json",
                         "--channel=Kpi",
                         "--components=CRSCF",
                         "--mixing",
                         "tests/current_result",
                         "tests/data/Kpi.root"],
    "mixing_fit_all": ["--config=tests/config.json",
                       "--channel=Kpi",
                       "--components=all",
                       "--mixing",
                       "tests/current_result",
                       "tests/data/Kpi.root"],
    "physpdf_lifetime_all": ["--config=tests/config.json",
                                      "--channel=Kpi",
                                      "--components=all",
                                      "--lifetime",
                                      "tests/current_result",
                                      "tests/data/Kpi.root"],
}

red_code = "\033[31m"
green_code = "\033[32m"
yellow_code = "\033[33m"
reset_code = "\033[0m"

def run_test(name, config):
    """Run a single test defined by its name and config."""
    print("Running " + name)
    config.insert(0, cmd)
    start = time.time()
    return_code = subprocess.call(config, cwd="../.")
    elapsed = (time.time() - start)

    reference_result = ""
    current_result = ""

    if (os.path.exists("references/" + name + ".reference")):
        with open("references/" + name + ".reference", "r") as f:
            reference_result = f.readline()
    else:
        print(yellow_code + "The result does NOT have a reference result!" + reset_code)
        os.rename("current_result", name + ".result")
        return 2, name, elapsed

    with open("current_result", "r") as f:
        current_result = f.readline()

    if return_code == 0 and reference_result == current_result:
        print(green_code + "The result matches the reference result." + reset_code)
        os.remove("current_result")
        # Delete results from previous failed runs if it now works
        if (os.path.exists(name + ".result")):
            os.remove(name + ".result")
        return 0, name, elapsed
    else:
        print(red_code + "The result does NOT match the reference result!" + reset_code)
        os.rename("current_result", name + ".result")
        return 1, name, elapsed


def print_results_table(results):
    """Print a summary of the results in a table."""
    print("\nSummary")
    print("-" * 34)

    num_failed = 0
    for result in test_results:
        result_str = ""
        if result[0] == 0:
            result_str = green_code + "pass" + reset_code
        elif result[0] == 2:
            result_str = yellow_code + "N/A" + reset_code
            num_failed += 1
        else:
            result_str = red_code + "FAIL" + reset_code
            num_failed += 1

        # The padding is 13 because python counts the ANSI escape sequences as chars
        print("{:<20} | {:>13} | {:>4.1f}".format(result[1], result_str, result[2]))

    print("-" * 34)
    print("{:<20} | {:>4} | {:>4.1f}".format("Total Failed & Time", num_failed, sum([x[2] for x in results])))


test_results = []

for test_name, test_config in test_configs.items():
    # Skip tests that weren't explicitely requested if some were
    if len(sys.argv) > 1 and test_name not in sys.argv[1:]:
        continue

    test_results.append(run_test(test_name, test_config))

test_results.sort()
print_results_table(test_results)


# Return code for Travis automation
if any([result[0] for result in test_results]):
    sys.exit(1)
else:
    sys.exit(0)
