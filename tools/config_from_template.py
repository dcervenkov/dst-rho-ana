#!/usr/bin/env python3
"""A script that takes a JSON config template and creates a final config.

The scripts looks for 'substitute' keys and replaces them with values from the
referenced files. Each 'substitute' key is an array of one or more JSON
objects. Each such object must have a 'file' key, which holds the name of a
file to be read-in. An object might also have a 'replace' key which specifies
that the newly read-in key names should be modified by replacing a substring of
each key name with a new string.

An example template:
{
    "modelParameters": {
        "substitute": [
            {
                "file": "Kpi_mc_scf.json",
                "replace": ["old", "new"]
            }
        ],
        "cr_f": 0.764,
        "scf_f": 0.094
    }
}
"""
import json
import os
import sys


def read_in_json(file):
    """Read in a JSON file; fix certain invalid JSON problems."""
    try:
        with open(file, 'r') as f:
            config = json.load(f)
        return config
    except json.decoder.JSONDecodeError:
        stripped_file = '{'
        with open(file, 'r') as f:
            raw_file = f.readlines()
            for line in raw_file:
                if not line.startswith('===') and not line.isspace():
                    stripped_file += line

        stripped_file = stripped_file.strip('\n,')
        stripped_file += '}'
        return json.loads(stripped_file)


def remove_extra_keys(dictionary):
    """Remove certain keys, that are relicts from the fit programs."""
    unwanted_keys = ['vrerr6', 'vterr6']
    for unwanted_key in unwanted_keys:
        dictionary.pop(unwanted_key, None)


def replace_key_substring(dictionary, old, new):
    """Recursively replace substrings in all keys in a dictionary."""
    if type(dictionary) is dict:
        for key in dictionary.keys():
            dictionary[key.replace(old, new)] = dictionary.pop(key)
            if type(dictionary[key.replace(old, new)]) is dict:
                dictionary[key.replace(old, new)] = replace_key_substring(
                    dictionary, old, new)

    return dictionary


def substitute_files(dictionary, from_dir):
    """Replace all 'substitute' keys in a dictionary.

    The keys are replaced based on the 'file' key. The filenames should be
    relative to from_dir.
    """
    for key, value in dictionary.copy().items():
        if isinstance(value, dict):
            substitute_files(value, from_dir)
        if key == 'substitute':
            del dictionary[key]
            for substitution in value:
                sub_dict = read_in_json(os.path.join(
                    from_dir, substitution['file']))
                remove_extra_keys(sub_dict)
                if 'replace' in substitution:
                    text_to_replace = substitution['replace'][0]
                    text_to_use = substitution['replace'][1]
                    sub_dict = replace_key_substring(
                        sub_dict, text_to_replace, text_to_use)

                dictionary.update(sub_dict)

    return dictionary


def main():
    template_path = sys.argv[1]
    sub_path = sys.argv[2]
    config_template = read_in_json(template_path)
    config_full = substitute_files(config_template, sub_path)
    print(json.dumps(config_full, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
