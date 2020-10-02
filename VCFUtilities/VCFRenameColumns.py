from VCF import VCF
from optparse import OptionParser

def IDmap_read(filename):
    if filename==None: raise Exception("Invalid filename. Please provide an IDmap file.")
    f = open(filename)
    header = f.readline().replace("\n", "").replace("\r", "").split(";")
    header[0] = header[0].replace("\ufeff", "")
    RNAtoVariety = {}
    GBStoVariety = {}
    for line in f:
        line = line.replace("\n", "").replace("\r", "")
        line = {name: value for name, value in zip(header, line.split(";"))}
        GBStoVariety[line['GBS']] = line['Variety']
        if line["RNA"]=='': continue
        if "," in line["RNA"]:
            samples = line["RNA"].split(",")
            for sample in samples:
                RNAtoVariety[sample] = line['Variety']
            continue
        RNAtoVariety[line["RNA"]] = line['Variety']
    return GBStoVariety, RNAtoVariety

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

if __name__=="__main__":
    usage = '''
    usage: python \033[4m%prog\033[24m \033[38;5;74m[options]\033[39m \033[32m<vcf file 1> <vcf file 2>\033[39m'''

    parser = OptionParser(usage)

    parser.add_option('-o', type="string", nargs=1, dest="output", default="output.csv",
                      help="Output file name. Default: output.csv")
    parser.add_option('--ID', type="string", nargs=1, dest="IDmap",
                      help="ID map file Default: IDmap.csv")
    parser.add_option('--VCF', type="string", nargs=1, dest="VCF",
                      help="Input vcf file")

    options, args = parser.parse_args()

    IDmap = options.IDmap
    VCF_filename = options.VCF
    output = options.output

    if IDmap==None: raise Exception("Please provide a IDmap file with --ID filename")
    GBStoVariety, RNAtoVariety = IDmap_read(IDmap)

    #if GBS==None: raise Exception("Please provide a GBS vcf file with --GBS filename")
    if VCF_filename==None: raise Exception("Please provide a vcf file with --VCF filename")

    VCF_file = VCF(VCF_filename)

    print(VCF_file.individuals)
    print(len(VCF_file.individuals))

    new_individual_names = []
    for rna_name, trans_name in zip(VCF_file.individuals, map(lambda x: RNAtoVariety.get(x, "SKIP"), VCF_file.individuals)):
        if trans_name=="SKIP": continue
        if trans_name in new_individual_names:
            continue
        new_individual_names.append(trans_name)
        i = VCF_file.header[rna_name]
        VCF_file.revheader[i] = trans_name
        del VCF_file.header[rna_name]
        VCF_file.header[trans_name] = i

    VCF_file.individuals = new_individual_names
    print(VCF_file.individuals)
    print(len(VCF_file.individuals))

    columns = [VCF_file.revheader[i] for i in range(0, len(VCF_file.revheader))]

    ## Find common genotypes between files, and set order of VCF file.
    #name_columns = VCF_file.individuals
    #standard_columns = ["CHROM", "POS"]
    #output_column_names = ["chromosome", "position"]

    #columns = standard_columns+sorted(name_columns)
    #out_columns = output_column_names+sorted(name_columns)

    print(columns)
    out_file = open(output, "w")
    out_file.write(VCF_file.header_info)
    out_file.write("#")
    out_file.write("\t".join(columns))
    out_file.write("\n")

    ## Initialize file reading
    VCFline = VCF_file.readline()

    SNPs_removed = 0
    SNPs_added = 0

    ## Iteratively read RNA and GBS file in order of SNPs
    while not VCF_file._empty:
        elements, line = VCFline
        out_file.write(line)
        VCFline = VCF_file.readline()
        if VCFline[0] is None: continue

    out_file.close()
