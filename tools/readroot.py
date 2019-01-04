#!/usr/bin/env python
"""This script prints a TNamed's Title field.

It's intended use is to easily print strings stored in ROOT files in TNamed
objects. (The recommended way to store text in ROOT files is to create a TNamed
object and store the text in its title member.)
"""

import argparse
import rootpy
import os
from rootpy.io import root_open
import sys
from contextlib import contextmanager


def decode_arguments():
    """Decode CLI arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("root_file")
    parser.add_argument("-n", "--name", type=str,
                        help="name of the object whose title to print")
    args = parser.parse_args()

    if not args.name:
        args.name = 'pull_table'

    return args.root_file, args.name


def root_print(filename, name):
    with root_open(filename, 'r') as root_file:
        if name in root_file:
            print(root_file[name].GetTitle())
        else:
            sys.exit(1)


@contextmanager
def stdout_redirected(to=os.devnull):
    '''
    import os

    with stdout_redirected(to=filename):
        print("from Python")
        os.system("echo non-Python applications are also supported")
    '''
    fd = sys.stdout.fileno()

    ##### assert that Python and C stdio write using the same file descriptor
    ####assert libc.fileno(ctypes.c_void_p.in_dll(libc, "stdout")) == fd == 1

    def _redirect_stdout(to):
        sys.stdout.close() # + implicit flush()
        os.dup2(to.fileno(), fd) # fd writes to 'to' file
        sys.stdout = os.fdopen(fd, 'w') # Python writes to fd

    with os.fdopen(os.dup(fd), 'w') as old_stdout:
        with open(to, 'w') as file:
            _redirect_stdout(to=file)
        try:
            yield # allow code to be run with the redirected stdout
        finally:
            _redirect_stdout(to=old_stdout) # restore stdout.
                                            # buffering and flags such as
                                            # CLOEXEC may be different

def main():
    filename, name = decode_arguments()

    # This crazy thing is in here to prevent RooFit from printing its intro
    # message, which can't be disabled (?!).
    with stdout_redirected():
        # Couldn't figure anything simpler to try to init RooFit. without other
        # side-effects as in the case of rootpy.ROOT.RooRealVar, etc.
        root_print(filename, name)

    root_print(filename, name)

main()
