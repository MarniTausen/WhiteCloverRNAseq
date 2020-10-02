from VCF import VCF
from optparse import OptionParser
import sys

def get_genotype(genotype):
    return genotype.split(":")[0]

def convert_genotype(genotype):
    if genotype=="./.": return "NA"
    g1, g2 = genotype.split("/")
    return str(int(g1)+int(g2))

def collect_elements(elements, columns):
    output_row = []
    missing_counter = 0
    for name in columns:
        if name=="CHROM" or name=="POS":
            output_row.append(elements[name])
            continue
        genotype = convert_genotype(get_genotype(elements[name]))
        output_row.append(genotype)
        if genotype=="NA":
            missing_counter += 1
    return output_row, missing_counter

def process_info(info):
    info_elements = {}
    for element in info.split(";"):
        name, item = element.split("=")
        info_elements[name] = item
    return info_elements

def process_header_for_info(header):
    info_names = []
    for line in header.split("\n"):
        if "##INFO" in line:
            name = line.split("=<")[-1].split(",")[0].split("=")[-1]
            info_names.append(name)
    return info_names

if __name__=="__main__":
    usage = '''
    usage: python \033[4m%prog\033[24m \033[38;5;74m[options]\033[39m \033[32m<vcf file 1> <vcf file 2>\033[39m'''

    parser = OptionParser(usage)

    parser.add_option('-o', type="string", nargs=1, dest="output", default="output.csv",
                      help="Output file name. Default: output.csv")
    parser.add_option('--RNA', type="string", nargs=1, dest="RNA",
                      help="RNA vcf file")
    parser.add_option('--missingness', type="float", nargs=1, dest="Missingness", default="1",
                      help="Percentage allowed missingness")

    options, args = parser.parse_args()

    RNA = options.RNA
    output = options.output
    missingness = options.Missingness

    #if GBS==None: raise Exception("Please provide a GBS vcf file with --GBS filename")
    if RNA==None: raise Exception("Please provide a RNA vcf file with --RNA filename")

    RNA_file = VCF(RNA)

    all_info_names = process_header_for_info(RNA_file.header_info)


    ignore_list = set(['DS', 'END', 'HaplotypeScore', 'RAW_MQ'])
    info_names = []
    for name in all_info_names:
        if name in ignore_list: continue
        info_names.append(name)

    print(info_names)

    print(RNA_file.individuals)
    print(len(RNA_file.individuals))

    missingness_cutoff = int(len(RNA_file.individuals)*missingness)

    print(missingness_cutoff)

    ## Find common genotypes between files, and set order of VCF file.
    name_columns = RNA_file.individuals
    standard_columns = ["CHROM", "POS", "S9", "IS_SNP"]

    genotype_list = set(["S9"])
    info_list = set(info_names)

    columns = standard_columns+sorted(info_names)

    print(columns)
    out_file = open(output, "w")
    out_file.write(";".join(columns))
    out_file.write("\n")

    ## Initialize file reading
    RNAline = RNA_file.readline()

    #SNPs_removed = 0
    #SNPs_added = 0


    ## Mark INDELs and multi allelic sites
    ## Iteratively read RNA and GBS file in order of SNPs
    while not RNA_file._empty:
        elements, _ = RNAline
        info_elements = process_info(elements["INFO"])
        output_row = []
        for name in columns:
            if name=="IS_SNP":
                if len(elements["REF"])>1 or len(elements["ALT"])>1:
                    output_row.append("False")
                else:
                    output_row.append("True")
                continue
            if name in genotype_list:
                output_row.append(convert_genotype(get_genotype(elements[name])))
            elif name in info_list:
                output_row.append(info_elements.get(name, "NA"))
            else:
                output_row.append(elements[name])
        out_file.write(";".join(output_row))
        out_file.write("\n")
        #SNPs_added += 1
        RNAline = RNA_file.readline()
        if RNAline[0] is None: continue

    out_file.close()

    #print("{} SNPs removed with missingness more than: {}%".format(SNPs_removed, missingness*100))
    #print("{} SNPs writen to file")
