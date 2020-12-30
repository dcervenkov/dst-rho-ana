#!/usr/bin/env python3
"""This script prints a TNamed's Title field.

It's intended use is to easily print strings stored in ROOT files in TNamed
objects. (The recommended way to store text in ROOT files is to create a TNamed
object and store the text in its title member.)
"""

import argparse
import sys

import uproot


def decode_arguments():
    """Decode CLI arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("files", metavar="FILE", nargs="+")
    parser.add_argument(
        "-n",
        "--name",
        type=str,
        action="append",
        help="name of the object whose title to print (can be used multiple times)",
    )
    parser.add_argument(
        "-l", "--list", help="list objects in file", action="store_true", default=False
    )
    args = parser.parse_args()

    if args.name is None:
        args.name = ["pull_table"]

    return args.files, args.name, args.list


def print_titles(filename, names):
    for name in names:
        value = None
        try:
            with uproot.open(filename) as root_file:
                value = root_file[name].member("fTitle")
        except KeyError:
            print(f"ERROR: '{name}' not present in file '{filename}'", file=sys.stderr)
            print_objects(filename)
            sys.exit(1)

        if len(names) > 1:
            if "\n" in value:
                print(name + ":")
            else:
                print(name + ":", end=" ")

        print(value)


def print_objects(filename):
    print("The following objects are present in the file:\n")
    with uproot.open(filename) as root_file:
        objects = [
            (name, classname)
            for name, classname in root_file.classnames().items()
        ]

    longest_name = len(max(objects, key=lambda obj: len(obj[0]))[0])
    longest_class = len(max(objects, key=lambda obj: len(obj[0]))[1])

    print(f"{'Name':{longest_name+1}}Class")
    print("-" * (longest_name + longest_class + 1))
    for name, classname in objects:
        print(f"{name:{longest_name+1}}{classname}")


def main():
    filenames, names, showlist = decode_arguments()

    for filename in filenames:
        if len(filenames) > 1:
            print("=" * 40)
            print(filename)
            print("=" * 40)
        if showlist:
            print_objects(filename)
        else:
            print_titles(filename, names)


main()
