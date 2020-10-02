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

    individual_count = {}
    for individual in VCF_data.individuals:
        individual_count[individual] = {}

    while not VCF_data._empty:
        elements, _ = VCFline
        for individual in VCF_data.individuals:
            genotype = elements[individual].split(":")[0]
            if genotype=="./.": continue
            if genotype not in individual_count[individual]:
                individual_count[individual][genotype] = 0
            individual_count[individual][genotype] += 1
        #sys.exit()
        VCFline = VCF_data.readline()

    #heterozyogisty_table = {}
    print("individual,heterozygisty")
    for individual in VCF_data.individuals:
        het = individual_count[individual].get("0/1", 0)+individual_count[individual].get("1/0", 0)
        hom = individual_count[individual].get("0/0", 0)+individual_count[individual].get("1/1", 0)
        print("{},{}".format(individual,het/float(hom+het)))
