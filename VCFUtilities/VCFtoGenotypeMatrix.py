from VCF import VCF
from optparse import OptionParser

def IDmap_read(filename):
    if filename==None: raise Exception("Invalid filename. Please provide an IDmap file.")
    f = open(filename)
    header = f.readline().replace("\n", "").replace("\r", "").split(";")
    header[0] = header[0].replace("\ufeff", "")
    RNAtoVariety = {}
    GBStoVariety = {}
    VarietytoRNA = {}
    for line in f:
        line = line.replace("\n", "").replace("\r", "")
        line = {name: value for name, value in zip(header, line.split(";"))}
        GBStoVariety[line['GBS']] = line['Variety']
        if line["RNA"]=='': continue
        if "," in line["RNA"]:
            samples = line["RNA"].split(",")
            VarietytoRNA[line['Variety']] = samples
            for sample in samples:
                RNAtoVariety[sample] = line['Variety']
            continue
        RNAtoVariety[line["RNA"]] = line['Variety']
        VarietytoRNA[line['Variety']] = line["RNA"]
    return GBStoVariety, RNAtoVariety, VarietytoRNA

def cleanGBSname(name):
    return name.split("/")[-1].split(".")[0]

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

def contains_resequence(elements):
    for element in elements:
        if "A" in element:
            return True
    return False

def coerce_as_list(x):
    if type(x)==list:
        return x
    return [x]

if __name__=="__main__":
    usage = '''
    usage: python \033[4m%prog\033[24m \033[38;5;74m[options]\033[39m \033[32m<vcf file 1> <vcf file 2>\033[39m'''

    parser = OptionParser(usage)

    parser.add_option('-o', type="string", nargs=1, dest="output", default="output.csv",
                      help="Output file name. Default: output.csv")
    parser.add_option('--ID', type="string", nargs=1, dest="IDmap",
                      help="ID map file Default: IDmap.csv")
    parser.add_option('--RNA', type="string", nargs=1, dest="RNA",
                      help="RNA vcf file")
    parser.add_option('--missingness', type="float", nargs=1, dest="Missingness", default="1",
                      help="Percentage allowed missingness")

    options, args = parser.parse_args()

    IDmap = options.IDmap
    RNA = options.RNA
    output = options.output
    missingness = options.Missingness

    if IDmap==None: raise Exception("Please provide a IDmap file with --ID filename")
    GBStoVariety, RNAtoVariety, VarietytoRNA = IDmap_read(IDmap)

    #print(RNAtoVariety)
    #print(VarietytoRNA)

    #if GBS==None: raise Exception("Please provide a GBS vcf file with --GBS filename")
    if RNA==None: raise Exception("Please provide a RNA vcf file with --RNA filename")

    RNA_file = VCF(RNA)

    print(RNA_file.individuals)
    print(len(RNA_file.individuals))

    replicate_counter = {}

    new_individual_names = []
    for rna_name, trans_name in zip(RNA_file.individuals, map(lambda x: RNAtoVariety.get(x, "SKIP"), RNA_file.individuals)):
        if trans_name=="SKIP": continue
        matchingnames = coerce_as_list(VarietytoRNA.get(trans_name))
        #print(rna_name, trans_name,  matchingnames)
        if len(matchingnames)>1:
            if contains_resequence(matchingnames):
                if "A" in rna_name:
                    trans_name += "_reseq"
                #print(rna_name, trans_name)
                new_individual_names.append(trans_name)
                i = RNA_file.header[rna_name]
                RNA_file.revheader[i] = trans_name
                del RNA_file.header[rna_name]
                RNA_file.header[trans_name] = i
                continue
            else:
                if trans_name not in replicate_counter:
                    replicate_counter[trans_name] = 0
                replicate_counter[trans_name] += 1

                trans_name += "_rep{}".format(replicate_counter[trans_name])

                #print(rna_name, trans_name)
                new_individual_names.append(trans_name)
                i = RNA_file.header[rna_name]
                RNA_file.revheader[i] = trans_name
                del RNA_file.header[rna_name]
                RNA_file.header[trans_name] = i

                continue
        #if trans_name in new_individual_names:
        #    continue
        #print(rna_name, trans_name)
        new_individual_names.append(trans_name)
        i = RNA_file.header[rna_name]
        RNA_file.revheader[i] = trans_name
        del RNA_file.header[rna_name]
        RNA_file.header[trans_name] = i

    RNA_file.individuals = new_individual_names
    #print(RNA_file.individuals)
    #print(len(RNA_file.individuals))

    missingness_cutoff = int(len(RNA_file.individuals)*missingness)

    #print(missingness_cutoff)

    ## Find common genotypes between files, and set order of VCF file.
    name_columns = RNA_file.individuals
    standard_columns = ["CHROM", "POS"]
    output_column_names = ["chromosome", "position"]

    columns = standard_columns+sorted(name_columns)
    out_columns = output_column_names+sorted(name_columns)

    print(columns)
    out_file = open(output, "w")
    out_file.write(",".join(out_columns))
    out_file.write("\n")

    ## Initialize file reading
    RNAline = RNA_file.readline()

    SNPs_removed = 0
    SNPs_added = 0

    ## Iteratively read RNA and GBS file in order of SNPs
    while not RNA_file._empty:
        elements, _ = RNAline
        output_row, missing_counter = collect_elements(elements, columns)
        if missing_counter>missingness_cutoff:
            RNAline = RNA_file.readline()
            SNPs_removed += 1
            continue
        out_file.write(",".join(output_row))
        out_file.write("\n")
        SNPs_added += 1
        RNAline = RNA_file.readline()
        if RNAline[0] is None: continue

    out_file.close()

    print("{} SNPs removed with missingness more than: {}%".format(SNPs_removed, missingness*100))
    print("{} SNPs writen to file")
