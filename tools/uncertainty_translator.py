import cmath
from math import sqrt, pi
import numpy as np
import matplotlib.pyplot as plt


def translate(hpm, hpa, hmm, hma):
    hp = cmath.rect(hpm, hpa)
    hm = cmath.rect(hmm, hma)

    ap = (hp + hm)/sqrt(2)
    at = (hp - hm)/sqrt(2)

    apm = abs(ap)
    apa = cmath.phase(ap)
    atm = abs(at)
    ata = cmath.phase(at)

    if ata < 0:
        ata += 2*pi
    return apm, apa, atm, ata


def main():
    hpm = 0.107
    shpm = 0.033
    hpa = 1.42
    shpa = 0.27
    hmm = 0.322
    shmm = 0.030
    hma = 0.31
    shma = 0.13

    sample_size = 100000

    hpms = np.random.normal(hpm, shpm, sample_size)
    hpas = np.random.normal(hpa, shpa, sample_size)
    hmms = np.random.normal(hmm, shmm, sample_size)
    hmas = np.random.normal(hma, shma, sample_size)

    translate_vectorized = np.vectorize(translate)
    results = translate_vectorized(hpms, hpas, hmms, hmas)

    exact_results = translate(hpm, hpa, hmm, hma)

    bins = np.linspace(0, 4, 200)

    for i, var in enumerate(["apm", "apa", "atm", "ata"]):
        print(f"{var} = {exact_results[i]:.4f} +/- {np.std(results[i]):.4f}")
        plt.hist(results[i], bins, alpha=0.5, label=var)

    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
