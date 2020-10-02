import argparse
import sys

import numpy as np

def load_as_matrix(filename):
    _file = open(filename)
    header = _file.readline().replace("\n", "").split(",")
    names = header[2:]
    rows = [list() for i in range(len(names))]
    for line in _file:
        line_info = line.replace("\n", "").split(",")
        for i, element in enumerate(line_info[2:]):
            if element=="NA":
                rows[i].append(0)
            else:
                rows[i].append(int(element))
    return names, np.array(rows)

def replace_column_nans_by_mean(matrix):
    # Set the value of gaps/dashes in each column to be the average of the other values in the column.
    nan_indices = np.where(np.isnan(matrix))
    # Note: bn.nanmean() instead of np.nanmean() because it is a lot(!) faster.
    column_nanmeans = np.nanmean(matrix, axis=0)
    matrix[nan_indices] = np.take(column_nanmeans, nan_indices[1])
    return matrix

def calcGRM(SNP):
    N, M = SNP.shape
    NORM = (SNP-np.mean(SNP, 0))/(np.std(SNP, 0)+0.000000000000001)
    return np.dot(NORM, NORM.T)/M

def save_GRM(GRM, names):
    GRM = GRM.astype(str)
    sys.stdout.write("names/names,"+",".join(names))
    sys.stdout.write("\n")
    for i, name in enumerate(names):
        row = GRM[i,:].tolist()
        sys.stdout.write("{},".format(name)+",".join(row))
        sys.stdout.write("\n")

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        usage="python MakeGRM.py <genotypefile> <MAFlevel> > grm.csv"
    )

    parser.add_argument('genotypefile')

    args = parser.parse_args(sys.argv[1:])

    names, SNP = load_as_matrix(args.genotypefile)

    GRM = calcGRM(SNP)

    save_GRM(GRM, names)
