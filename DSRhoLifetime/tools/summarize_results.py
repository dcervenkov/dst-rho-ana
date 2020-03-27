#!/usr/bin/env python3
import glob
from uncertainties import ufloat


def load_values(paths):
    """Load values from single-line CSV files

    Results are returned in the form of
    [
        [
            [val1_1, err1_1],
            [val1_2, err1_2]
        ],
        [
            [val2_1, err2_1],
            [val2_2, err2_2]
        ]
    ]
    """
    num_values = get_num_values(paths[0])
    results = [[] for _ in range(num_values)]
    for path in paths:
        with open(path, "r") as f:
            line = f.readline()
            value_row = [value.strip(' \n') for value in line.split(',')]
            # Make sure all files have the same number of values
            assert len(value_row) / 2 == num_values
            for val_no in range(num_values):
                results[val_no].append(
                    [float(value_row[2 * val_no]), float(value_row[2 * val_no + 1])])

    return results


def get_num_values(path):
    """Get number of (value, error) pairs in single-line CSV file"""
    with open(path, "r") as f:
        line = f.readline()
        value_row = [value.strip(' \n') for value in line.split(',')]
        # Make sure each number has (value, error)
        assert len(value_row) % 2 == 0
        return int(len(value_row) / 2)


def calc_averages(results):
    """Simple average calculation for a list of (val, err)

    Both are averaged separately, no error propagation occurs. This
    corresponds to completely correlated variables.
    """
    averages = [0, 0]
    for result in results:
        averages[0] += result[0] / len(results)
        averages[1] += result[1] / len(results)
    return averages


def calc_averages_w_errors(results):
    """Average calculation for a list of (val, err)

    Errors are propagated assuming uncorrelated variables.
    """
    uresults = []
    for result in results:
        uresults.append(ufloat(result[0], result[1]))

    return sum(uresults) / len(uresults)


def main():
    for type in ["lifetime", "mixing"]:
        print(type)
        averages = {
            "Kpi": {
                "CR": [],
                "CRSCF": [],
                "all": []
            },
            "Kpipi0": {
                "CR": [],
                "CRSCF": [],
                "all": []
            },
            "K3pi": {
                "CR": [],
                "CRSCF": [],
                "all": []
            },
            "together": {
                "CR": [],
                "CRSCF": [],
                "all": []
            }
        }

        components = ["CR", "CRSCF", "all"]
        for channel in ["Kpi", "Kpipi0", "K3pi", "together"]:
            print(channel)
            for component in components:
                print(f"{component:6}", end=" ")
                helper = "" if component == "all" else "0"
                for values_w_errors in load_values(
                        glob.glob(f"../results/{channel}/{component}_{type}_{helper}[0-5]")):
                    average = calc_averages(values_w_errors)
                    averages[channel][component].append(average)
                    print(f"{average[0]:.3f} +- {average[1]:.3f}", end="\t")
                print()
            print()

        print("Combined")
        for component in components:
            print(f"{component:6}", end=" ")
            for i in range(len(averages["Kpi"]["CR"])):
                average = calc_averages_w_errors(
                    [averages[channel][component][i] for channel in ["Kpi", "Kpipi0", "K3pi"]])
                print(f"{average.n:.3f} +- {average.s:.3f}", end="\t")
            print()
        print()


if __name__ == "__main__":
    main()
