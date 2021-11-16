'''
This script calculates descriptive statistics for the molecular fingerprint, proposed in the paper ***Topological indices in assessing molecular similarity***. Additionaly, it provides the same statistical parameters for the Morgan circular fingerprint as well.
'''

import json
import click
from statistics import mean, stdev
from progressbar import progressbar


def tanimoto(first_fp, second_fp):
    '''Calculates the Tanimoto similarity coeficient between the fingerprints
    of the **first** and **second** molecule.

    Parameters
    ----------
    first_fp : str
    second_fp : str

    Returns
    -------
    tanimoto_coeff : float
    '''
    a = first_fp.count('1')
    b = second_fp.count('1')
    c = 0
    for i, j in zip(first_fp, second_fp):
        if i == '1' and j == '1':
            c += 1

    return round(c / (a + b - c), 5)


@click.command()
@click.argument('path_to_file')
def main(path_to_file):
    with open(path_to_file) as infile:
        data = json.load(infile)

    fps = []
    mfps = []
    for molecule in data:
        fps.append(molecule['new_fp'])
        mfps.append(molecule['morgan_fp'])
    tanimotos = []
    tanimotos_morgan = []
    same_nfp = 0
    same_mfp = 0
    nbr = 0

    for i in progressbar(range(len(fps) - 1)):
        for j in range(i + 1, len(fps)):
            nbr += 1
            if mfps[i] == mfps[j]:
                same_mfp += 1
            if fps[i] == fps[j]:
                same_nfp += 1
            tanimotos.append(tanimoto(fps[i], fps[j]))
            tanimotos_morgan.append(tanimoto(mfps[i], mfps[j]))

    print(
        f"NEW FP:\nmean: {mean(tanimotos):.1%}\nmax: {max(tanimotos):.1%}\nmin: {min(tanimotos):.1%}\nsd: {stdev(tanimotos):.1%}\nIdentical: {same_nfp/nbr*10000:.2f} permyriads.\n"
    )

    print(
        f"MORGAN FP:\nmean: {mean(tanimotos_morgan):.1%}\nmax: {max(tanimotos_morgan):.1%}\nmin: {min(tanimotos_morgan):.1%}\nsd: {stdev(tanimotos_morgan):.1%}\nIdentical: {same_mfp/nbr*10000:.2f} permyriads.\n"
    )


if __name__ == '__main__':
    main()