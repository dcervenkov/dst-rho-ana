#!/usr/bin/env python3

import root_pandas
import sys
import math


def main():
    df = root_pandas.read_root(sys.argv[1])
    df["thetat"] = math.pi - df["thetat"]
    df["phit"] = -df["phit"]
    df.to_root(sys.argv[2], key="h2000")


if __name__ == "__main__":
    main()
