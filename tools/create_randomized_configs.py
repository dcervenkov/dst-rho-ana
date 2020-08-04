#!/usr/bin/env python3

import argparse
import os
import subprocess


def decode_arguments():
    """Decode CLI arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--config-dir",
        required=True,
        help="where the generated configs should be stored",
    )
    parser.add_argument(
        "-r",
        "--search-and-replace",
        required=True,
        help="text to search for and replace",
        action="append",
    )
    parser.add_argument(
        "-n",
        "--number",
        required=True,
        type=int,
        help="number of configs to be generated",
    )
    parser.add_argument("template")
    args = parser.parse_args()
    return args


def read_in_template(file):
    """Read in a template file."""
    with open(file, "r") as f:
        template = f.readlines()
    return template


def get_replacement_line(line, to_replace, replace_by, number):
    start = line.find(to_replace)
    if start != -1:
        new_line = (
            line[:start] + replace_by.format(number) + line[start + len(to_replace) :]
        )
        return new_line
    else:
        return line


def create_config_template_path(template_path, config_dir, count):
    basename = os.path.basename(template_path)
    return os.path.join(
        config_dir, basename.replace(".template", "_" + str(count) + ".template")
    )


def main():
    args = decode_arguments()
    template = read_in_template(args.template)

    for i in range(args.number):
        config_template_path = create_config_template_path(
            args.template, args.config_dir, i
        )
        if not os.path.isdir(args.config_dir):
            os.mkdir(args.config_dir)
        with open(config_template_path, "w") as f:
            for line in template:
                for sar in args.search_and_replace:
                    search, replace = sar.split(',')
                    line = get_replacement_line(line, search, replace, i)
                f.write(line)
        with open(config_template_path.replace(".template", ""), "w") as f:
            subprocess.call(
                ["./config_from_template.py", config_template_path], stdout=f
            )
        os.remove(config_template_path)


if __name__ == "__main__":
    main()
