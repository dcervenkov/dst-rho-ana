#!/usr/bin/env python
"""This script prints a TNamed's Title field.

It's intended use is to easily print strings stored in ROOT files in TNamed
objects. (The recommended way to store text in ROOT files is to create a TNamed
object and store the text in its title member.)
"""

from __future__ import print_function
import argparse
import rootpy
import os
from rootpy.io import root_open
import sys
from contextlib import contextmanager


def decode_arguments():
    """Decode CLI arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("file", metavar='ROOT-FILE', nargs='+')
    parser.add_argument("-n", "--name", type=str, action='append', default=['pull_table'],
                        help="name of the object whose title to print (can be used multiple times)")
    args = parser.parse_args()

    return args.file, args.name


def root_print(filename, names):
    with root_open(filename, 'r') as root_file:
        for name in names:
            if name in root_file:
                value = root_file[name].GetTitle()
            else:
                print("ERROR: '" + name + "' not present in file '" + filename + "'", file=sys.stderr)
                print("The following objects are present in the file:\n" + "\n".join([obj.GetName() for obj in root_file]), file=sys.stderr)
                sys.exit(1)

            if len(names) > 1:
                if '\n' in value:
                    print(name + ':')
                else:
                    print(name + ':', end=' ')

            print(value)


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
    filenames, names = decode_arguments()

    # This crazy thing is in here to prevent RooFit from printing its intro
    # message, which can't be disabled (?!).
    with stdout_redirected():
        # Couldn't figure anything simpler to try to init RooFit. without other
        # side-effects as in the case of rootpy.ROOT.RooRealVar, etc.
        root_print(filenames[0], names)

    for filename in filenames:
        if len(filenames) > 1:
            print("=" * 40)
            print(filename)
        root_print(filename, names)

main()
