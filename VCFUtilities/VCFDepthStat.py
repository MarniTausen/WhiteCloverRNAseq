from VCF import VCF
import argparse
import sys

def process_info(info):
    info_elements = {}
    for element in info.split(";"):
        name, item = element.split("=")
        info_elements[name] = item
    return info_elements

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        usage="python VCFfilterByGenotypeFile.py <vcf_file> <filter_by_genotype_file.csv> > filtered_vcf_file.csv"
    )

    parser.add_argument('vcffile')

    args = parser.parse_args(sys.argv[1:])

    VCF_data = VCF(args.vcffile)

    ## Initialize file reading
    VCFline = VCF_data.readline()

    Depth = []

    while not VCF_data._empty:
        elements, line = VCFline
        info_elements = process_info(elements["INFO"])
        if "DP" in info_elements:
            Depth.append(float(info_elements["DP"]))
        VCFline = VCF_data.readline()

    print(sum(Depth)/len(Depth))
