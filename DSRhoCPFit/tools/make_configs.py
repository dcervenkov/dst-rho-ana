#!/usr/bin/env python3

import json


def generate_ranges(begin, end, steps):
    ranges = []
    interval = (end - begin)/steps
    for step in range(steps):
        if step == 0:
            ranges.append([[begin + interval, end]])
        elif step == steps - 1:
            ranges.append([[begin, end - interval]])
        else:
            ranges.append([[begin, begin + step * interval],
                           [begin + (step + 1) * interval, end]])

    return ranges


def make_configs(var_name, begin, end, steps, config_template, filename_base):
    with open(config_template, "r") as f:
        json_template = json.load(f)

    for i, var_ranges in enumerate(generate_ranges(begin, end, steps)):
        config = json_template
        config["fitRanges"][var_name] = var_ranges

        # print(json.dumps(config, indent=4))
        file_name = filename_base + str(i + 1) + ".json"
        with open(file_name, "w") as f:
            json.dump(config, f, indent=4)


NUM_STEPS = 5
CONFIG_TEMPLATE = "config_template.json"

make_configs("thetat", 0, 3.1415, NUM_STEPS, CONFIG_TEMPLATE, "cfg_exc_thetat_")
make_configs("thetab", 0.5, 2.95, NUM_STEPS, CONFIG_TEMPLATE, "cfg_exc_thetab_")
make_configs("phit", -3.1415, 3.1415, NUM_STEPS, CONFIG_TEMPLATE, "cfg_exc_phit_")