#!/usr/bin/env python3
"""A module that runs fit tests (end-to-end tests).

This module runs predefined fit tests and compares their results to reference
ones. It outputs a summary table and quits with the exit code 0 when all tests
pass and 1 when any of them fail.

Tests are supplied via JSON files. Each JSON test config file must have the
following structure:
{
    "command": "name_of_program_to_test",
    "comparison_file": "file_that_will_be_compared_to_reference",
    "temporary_paths": [
        "files_that_should_be_deleted",
        "after_each_test"
    ],
    "post_commands": [
        [
            "a_command",
            "with_arguments",
            "to_be_run_after_running",
            "the_test_itself",
            "to_ie_extract_the_results_from_an_interstage_file"
        ]
    ],
    "tests": {
        "name_of_test": [
            "arguments",
            "to_the_command"
        ]
    }
}

temporary_paths and post_commands are optional.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time

RED_CODE = "\033[31m"
GREEN_CODE = "\033[32m"
YELLOW_CODE = "\033[33m"
RESET_CODE = "\033[0m"


def decode_arguments():
    """Decode CLI arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("tests", nargs="*")
    parser.add_argument("-f", "--config", action="append",
                        help="test config in JSON format [default: tests.json]",
                        default=[])
    args = parser.parse_args()

    if not args.config:
        args.config.append("tests.json")

    return args.config, args.tests


def run_test(config):
    """Run a single test defined by its name and config."""
    name = config['name']
    test_config = config['config']

    print("Running " + name)
    test_config.insert(0, config['command'])
    start = time.time()
    return_code = subprocess.call(test_config, cwd="../.")
    elapsed = (time.time() - start)

    if 'post_commands' in config:
        for post_command in config['post_commands']:
            subprocess.call(post_command)

    if 'temporary_paths' not in config:
        config['temporary_paths'] = []

    reference_result = ""
    current_result = ""

    if not os.path.exists(config['comparison_file']):
        print(RED_CODE + "Test failed to create comparison file!" + RESET_CODE)
        delete_paths(config['temporary_paths'])
        return 3, name, elapsed

    if os.path.exists("references/" + name + ".reference"):
        with open("references/" + name + ".reference", "r") as f:
            reference_result = f.read()
    else:
        print(YELLOW_CODE + "Test does NOT have a reference result!" + RESET_CODE)
        os.rename(config['comparison_file'], name + ".result")
        delete_paths(config['temporary_paths'])
        return 2, name, elapsed

    with open(config['comparison_file'], "r") as f:
        current_result = f.read()

    if return_code == 0 and reference_result == current_result:
        print(GREEN_CODE + "Result matches reference result." + RESET_CODE)
        os.remove(config['comparison_file'])
        # Delete results from previous failed runs if it now works
        if os.path.exists(name + ".result"):
            os.remove(name + ".result")
        delete_paths(config['temporary_paths'])
        return 0, name, elapsed
    else:
        print(RED_CODE + "Result does NOT match reference result!" + RESET_CODE)
        os.rename(config['comparison_file'], name + ".result")
        delete_paths(config['temporary_paths'])
        return 1, name, elapsed


def delete_paths(paths):
    for path in paths:
        if os.path.isfile(path):
            os.remove(path)
        elif os.path.isdir(path):
            shutil.rmtree(path)
        else:
            print("Can't delete " + path)


def print_results_table(results):
    """Print a summary of the results in a table."""
    print("\nSummary")
    print("-" * 34)

    num_failed = 0
    for result in test_results:
        result_str = ""
        if result[0] == 0:
            result_str = GREEN_CODE + "pass" + RESET_CODE
        elif result[0] == 2:
            result_str = YELLOW_CODE + "N/A" + RESET_CODE
            num_failed += 1
        elif result[0] == 1:
            result_str = RED_CODE + "FAIL" + RESET_CODE
            num_failed += 1
        else:
            result_str = RED_CODE + "ERR " + RESET_CODE
            num_failed += 1

        # The padding is 13 because python counts the ANSI escape sequences as chars
        print("{:<20} | {:>13} | {:>4.1f}".format(
            result[1], result_str, result[2]))

    print("-" * 34)
    print("{:<20} | {:>4} | {:>4.1f}".format(
        "Total Failed & Time", num_failed, sum([x[2] for x in results])))


def load_tests_to_run(config_files, tests_to_run):
    """Create a dictionary of tests from JSON input files.

    If tests_to_run are specified, only those tests will be run.
    """
    all_test_configs = []
    for config_file in config_files:
        with open(config_file, 'r') as file:
            test_configs = json.load(file)
            # Skip tests that weren't explicitely requested if some were
            for test in test_configs['tests'].keys():
                if tests and test not in tests:
                    continue
                test_entry = {}
                test_entry['name'] = test
                test_entry['command'] = test_configs['command']
                test_entry['comparison_file'] = test_configs['comparison_file']
                test_entry['config'] = test_configs['tests'][test]
                if 'temporary_paths' in test_configs:
                    test_entry['temporary_paths'] = test_configs['temporary_paths']
                if 'post_commands' in test_configs:
                    test_entry['post_commands'] = test_configs['post_commands']
                all_test_configs.append(test_entry)

    return all_test_configs


config_files, tests = decode_arguments()

test_configs = load_tests_to_run(config_files, tests)
test_results = []
for test_config in test_configs:
    test_results.append(run_test(test_config))

test_results.sort()
print_results_table(test_results)

# Return code for Travis automation
if any([result[0] for result in test_results]):
    sys.exit(1)
else:
    sys.exit(0)
