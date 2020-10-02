import argparse
import sys

def read_filter_file(filename):
    _file = open(filename)
    items = set(["chromosome", "position"])
    for line in _file:
        items.add(line.replace("\n",""))
    return items


if __name__=="__main__":
    parser = argparse.ArgumentParser(
        usage="python choose_genotypes.py <file_to_filter>.csv <file_to_filter_by>.txt > filtered_genotype_file"
    )

    parser.add_argument('genotypefile')
    parser.add_argument('samples_file')

    args = parser.parse_args(sys.argv[1:])

    genotypes_to_keep = read_filter_file(args.samples_file)

    #print(genotypes_to_keep)

    _file = open(args.genotypefile)

    header = _file.readline().replace("\n", "").split(",")

    indices_to_keep = set()
    for i, name in enumerate(header):
        if name in genotypes_to_keep:
            indices_to_keep.add(i)

    filtered_header = []
    for i, name in enumerate(header):
        if i in indices_to_keep:
            filtered_header.append(name.replace("_reseq","").replace("_rep1","").replace("_rep2", "").replace("_rep3",""))

    print(','.join(filtered_header))

    for line in _file:
        line_info = line.replace("\n", "").split(",")
        row = []
        for i, element in enumerate(line_info):
            if i in indices_to_keep:
                row.append(element)
        print(','.join(row))
