from VCF import VCF
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
        usage="python VCFfilterByGenotypeFile.py <vcf_file> <filter_by_genotype_file.csv> > filtered_vcf_file.csv"
    )

    parser.add_argument('vcffile')
    parser.add_argument('filterfile')

    args = parser.parse_args(sys.argv[1:])

    genotypes_to_keep = read_filter_file(args.filterfile)

    VCF_data = VCF(args.vcffile, colnames=True)
    for line in VCF_data.header_info.split("\n"):
        if "fileformat" in line:
            print("##fileformat=VCFv4.1")
            continue
        if "##"==line[:2]:
            print(line)
    #print(VCF_data.header_info, end="")

    genotypes_list = []
    for name  in VCF_data.individuals:
        if name in genotypes_to_keep:
            genotypes_list.append(name)

    columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    output_columns = columns+genotypes_list

    print("#"+"\t".join(output_columns))


    ## Initialize file reading
    VCFline = VCF_data.readline()

    ## Iteratively every line in the RNA file
    while not VCF_data._empty:
        elements, line = VCFline
        row = []
        for out_col in output_columns:
            row.append(elements[out_col])
        print("\t".join(row))
        VCFline = VCF_data.readline()
