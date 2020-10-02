from optparse import OptionParser
import sys
import os

class VCF:

    def __init__(self, filename):
        self._file = open(filename)
        self._empty = False
        self.read_header()
        self.find_pairs()

    def read_header(self):
        self.header_info = ""
        while True:
            self._last_line = self._file.tell()
            line = self._file.readline()
            if line[:2]=="##":
                self.header_info += line
                continue
            if line[0]=="#":
                self.header = {}
                self.revheader = {}
                self.header_info += line
                line = line.replace("\n","")
                for i, name in enumerate(line[1:].split("\t")):
                    self.header[name] = i
                    self.revheader[i] = name
                continue
            self._file.seek(self._last_line)
            break

    def find_pairs(self):
        exclude = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                   "FILTER", "INFO", "FORMAT"]
        self.individuals = []
        self.individual_pair = {}
        for name in self.header:
            if name in exclude: continue
            if ":" in name:
                true_name = name.split(":")[-1]
                self.individual_pair[true_name] = self.header[name]
                continue
            self.individuals.append(name)

    def readline(self):
        if self._empty: return None, ''
        line = self._file.readline()
        if line=='':
            self._empty = True
            return None, ''
        sline = line.split("\t")
        sline[-1] = sline[-1].replace("\n", "")
        elements = {}
        for i, column in enumerate(sline):
            elements[self.revheader[i]] = column
        return elements, line

    def __iter__(self):
        return self

    def __next__(self):
        return self.readline()

def clean(genotype):
    return (genotype[0], genotype[2])

if __name__=='__main__':
    usage = '''
    usage: python \033[4m%prog\033[24m \033[38;5;74m[options]\033[39m \033[32m<vcf file>\033[39m'''

    parser = OptionParser(usage)

    options, args = parser.parse_args()

    if args==[]:
        raise Exception("No vcf file was included. Please provide a vcf file as the first argument")
    VCF_file = VCF(args[0])
    sys.stdout.write(VCF_file.header_info)

    no_genotype = (".", ".")

    for info, line in VCF_file:
        if info is None: break
        genotype = info['S9'].split(":")[0]
        #if genotype!='0/0' and genotype!='./.': continue
        if genotype!='0/0': continue
        sys.stdout.write(line)
