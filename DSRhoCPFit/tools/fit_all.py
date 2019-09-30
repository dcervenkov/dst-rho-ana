#!/usr/bin/env python3

import argparse
import os
import shutil
import subprocess
import sys

DEFAULT_OPTIONS = ["--MC=1",
                   "--cpus=1",
                   "--log"]
COMMAND = "./DSRhoCPFit"
CHANNELS = ["Kpi", "Kpipi0", "K3pi"]
STREAMS = range(0, 6)
CONFIG_BASENAME = "config_stream"
RESULTS_DIR = "results"
TEMP_FILE = "/tmp/temp_commands"


def decode_arguments():
    """Decode CLI arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--channel", action="append",
                        help="channel to be processed")
    parser.add_argument("-d", "--td", action="store_true",
                        help="do time-dependent fits")
    parser.add_argument("-i", "--ti", action="store_true",
                        help="do time-independent fits")
    parser.add_argument("-m", "--component",
                        action="append", help="do CR-only fits")
    parser.add_argument("-f", "--force", action="store_true",
                        help="force rewrite of existing results")
    parser.add_argument(
        "-j", "--cpus", help="number of cores to use", default="1")
    args = parser.parse_args()
    return {"channels": args.channel,
            "components": args.component,
            "ti": args.ti,
            "td": args.td,
            "force": args.force,
            "cpus": args.cpus}


def are_valid(options):
    """Check if the options dict is valid"""
    return True


def create_fit_cmds(options):
    """Fit one channel with all the components, fit types, etc."""
    commands = []
    time_types = []
    if (options["ti"]):
        time_types.append("ti")
    if (options["td"]):
        time_types.append("td")

    dirs_to_create = []
    dirs_to_delete = []
    dirs_to_skip = []
    for time_type in time_types:
        for channel in options["channels"]:
            for component in options["components"]:
                output_dir = os.path.join(RESULTS_DIR, "_".join(
                    [channel.replace(",", "-"), time_type, component]))
                if os.path.exists(os.path.join("..", output_dir)):
                    if options["force"]:
                        dirs_to_delete.append(output_dir)
                    else:
                        dirs_to_skip.append(output_dir)
                        continue
                else:
                    dirs_to_create.append(output_dir)

                for stream in STREAMS:
                    cmd = []
                    cmd.append(COMMAND)
                    cmd += DEFAULT_OPTIONS
                    cmd.append("--components=" + component)
                    cmd.append("--config=" + CONFIG_BASENAME +
                               str(stream) + ".json")
                    excluded_channels = set(CHANNELS) - set(channel.split(','))
                    if excluded_channels:
                        cmd.append("--exclude-channels=" +
                                   ",".join(excluded_channels))
                    output_file = os.path.join(
                        output_dir, "stream" + str(stream))
                    cmd.append("--output=" + output_file)
                    if time_type == "ti":
                        cmd.append("--time-independent")

                    commands.append(cmd)

    return commands, dirs_to_create, dirs_to_delete, dirs_to_skip


def main():
    options = decode_arguments()
    if not are_valid(options):
        sys.exit(1)

    commands, dirs_to_create, dirs_to_delete, dirs_to_skip = create_fit_cmds(
        options)

    print("Will run the following commands:")
    for command in commands:
        print(" ".join(command))

    if (dirs_to_skip):
        print("\nWill skip the following dirs:")
        for dir_to_skip in dirs_to_skip:
            print(dir_to_skip)

    if (dirs_to_create):
        print("\nWill create the following dirs:")
        for dir_to_create in dirs_to_create:
            print(dir_to_create)

    if (dirs_to_delete):
        print("\nWill delete the following dirs:")
        for dir_to_delete in dirs_to_delete:
            print(dir_to_delete)

    print("\nWill run on " + options["cpus"] + " core(s)")

    user_input = input("\nContinue? [Y/n]: ")
    if (user_input == "" or user_input.lower() == "y"):
        with open(TEMP_FILE, 'w') as f:
            for command in commands:
                f.write(" ".join(command) + "\n")

        for dir in dirs_to_create:
            os.mkdir(os.path.join("..", dir))
        for dir in dirs_to_delete:
            shutil.rmtree(os.path.join("..", dir))

        subprocess.run(
            "parallel --nice 10 -j " + options["cpus"] + " < " + TEMP_FILE, cwd="..", shell=True)
        os.remove(TEMP_FILE)
    else:
        print("Aborting")
        sys.exit(2)


if __name__ == "__main__":
    main()
