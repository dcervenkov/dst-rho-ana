#!/usr/bin/env python3

import copy
import json
import sys


def create_modified_values(variations_path):
    multiplier = 3
    modified_values = []
    with open(variations_path, "r") as f:
        variations = json.load(f)
        for par in variations['parameters'].items():
            par_name = par[0]
            par_value = par[1][0]
            par_error = par[1][1]
            modified_values.append({par_name: par_value + multiplier * par_error})
            modified_values.append({par_name: par_value - multiplier * par_error})

    return modified_values

def create_modified_config(org_config, modified_value):
    new_config = copy.deepcopy(org_config)
    for key, value in modified_value.items():
        new_config['modelParameters'][key] = value
    return new_config


def main():
    original_config_path = sys.argv[1]
    variations_path = sys.argv[2]

    with open(original_config_path, "r") as org_f:
        original_config = json.load(org_f)
    
        modified_values = create_modified_values(variations_path)
        for i, modified_value in enumerate(modified_values):
            modified_config = create_modified_config(original_config, modified_value)
            with open("var_config_{}.json".format(i), "w") as new_f:
                new_f.write(json.dumps(modified_config, indent=4))




main()