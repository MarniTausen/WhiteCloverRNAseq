White Clover RNAseq supplementary scripts
================
28/09/2020 - 14:17:41

-   [Introduction](#introduction)
-   [RNA mapping and variant calling](#rna-mapping-and-variant-calling)
    -   [Add readgroups to BAM files](#add-readgroups-to-bam-files)
    -   [Mark duplicates and correct quality](#mark-duplicates-and-correct-quality)
    -   [SplitNTrim RNA reads](#splitntrim-rna-reads)
    -   [Variant calling (Unifiedgenotyper)](#variant-calling-unifiedgenotyper)
    -   [Variant calling (HaplotypeCaller)](#variant-calling-haplotypecaller)
    -   [Variant calling (HaplotypeCallerGVCF)](#variant-calling-haplotypecallergvcf)
    -   [FilterVariants](#filtervariants)

Introduction
============

This document lists the scripts and workflows/pipelines used for performing mapping and variant calling on white clover RNAseq data. Each step of the process has its own section and code in the correct order it was run.

RNA mapping and variant calling
===============================

Read mapping of the RNAseq data. The workflow script is a gwf pipeline (<http://gwf.readthedocs.io/en/latest/>), which keeps track of all job dependencies and submits the scripts in the correct order.

``` python
from gwf import Workflow
import os
from math import ceil
gwf = Workflow()
#################################################################### STAR
def star_index(reference):
    dir = "/".join(reference.split("/")[:-1])
    inputs = [reference]
    outputs = [dir+"/genomeParameters.txt"]
    options = {"cores":8, "memory":"64g", "account":"NChain", "walltime": "12:00:00"}
    spec = """
    source /com/extra/STAR/2.5.2b/load.sh
    STAR --runMode genomeGenerate --runThreadN 8 --genomeDir {dir} --genomeFastaFiles {ref}
    """.format(ref=reference, dir=dir)
    return inputs, outputs, options, spec
def star_mapping(read1, read2, output, reference):
    dir = "/".join(reference.split("/")[:-1])
    inputs = [dir+"/genomeParameters.txt", read1, read2]
    outputs = ["./"+output+"Aligned.sortedByCoord.out.bam"]
    options = {
        "cores": 4,
        "memory": "16g",
        "account": "NChain",
        "walltime": "12:00:00"}
    OFilePrefix = "./"+output
    spec = """
    source /com/extra/STAR/2.5.2b/load.sh
    STAR --runThreadN 8 --genomeDir {dir} --readFilesCommand zcat --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within KeepPairs --outFileNamePrefix {of} --readFilesIn {r1} {r2}
    """.format(r1 = read1, r2=read2, of=OFilePrefix, dir=dir)
    return inputs, outputs, options, spec
def add_readgroups(bam, rg_bam, rg):
     inputs = [bam]
     outputs = [rg_bam]
     options = {"cores": 1,
               "memory": "4g",
               "account": "NChain",
               "walltime": "12:00:00"}
     spec = """./picard.sh {bam} {rg_bam} {rg}
     """.format(bam=bam, rg_bam=rg_bam, rg=rg)
     return inputs, outputs, options, spec
def mark_duplicates(in_bam, out_bam):
    inputs = [in_bam]
    outputs = [out_bam]
    options = {"cores": 1,
              "memory": "8g",
              "account": "NChain",
              "walltime": "12:00:00"}
    spec = """./markduplicates.sh {bam} {out_bam}
    """.format(bam=in_bam, out_bam=out_bam)
    return inputs, outputs, options, spec
def samtools_index(bam_file):
     inputs = [bam_file]
     outputs = [bam_file+".bai"]
     options = {"cores": 1,
               "memory": "4g",
               "account": "NChain",
               "walltime": "12:00:00"}
     spec = """
     source /com/extra/samtools/1.6.0/load.sh
     samtools index {input}""".format(input=bam_file)
     return inputs, outputs, options, spec
def fix_contig_len(in_bam, out_bam):
    inputs = [in_bam]
    outputs = [out_bam]
    options = {"cores": 1,
               "memory": "4g",
               "account": "NChain",
               "walltime": "12:00:00"}
    spec = """
    source /com/extra/samtools/1.6.0/load.sh
    samtools view -h {bam} | python sam_reference_fix.py | samtools view -S -b -h > {out_bam}
    """.format(bam=in_bam, out_bam=out_bam)
    return inputs, outputs, options, spec
def SplitNTrim(in_bam, out_bam):
    inputs = [in_bam, in_bam+".bai"]
    outputs = [out_bam]
    options = {"cores": 1,
              "memory": "32g",
              "account": "NChain",
              "walltime": "12:00:00"}
    spec = """./SplitNTrim.sh {bam} {out_bam}
    """.format(bam=in_bam, out_bam=out_bam)
    return inputs, outputs, options, spec
def make_bam_list(bam_files, output):
    inputs = []
    for bam_file in bam_files:
        inputs.append(bam_file)
        inputs.append(bam_file+".bai")
    outputs = [output]
    options = {"cores": 1,
              "memory": "2g",
              "account": "NChain",
              "walltime": "1:00:00"}
    spec = ""
    for bam_file in bam_files:
        spec += "echo $PWD/{} >> {}\n".format(bam_file, output)
    return inputs, outputs, options, spec
def unifiedgenotyper(bam_list, output):
     inputs = [bam_list]
     outputs = [output]
     options = {"cores": 12,
               "memory": "32g",
               "account": "NChain",
               "walltime": "72:00:00"}
     spec = """
     ./Unifiedgenotyper.sh {} {}
     """.format(bam_list, output)
     return inputs, outputs, options, spec
def HaplotypeCaller(bam_list, output):
    inputs = [bam_list]
    outputs = [output]
    options = {"cores": 16,
               "memory": "128g",
               "account": "NChain",
               "walltime": "360:00:00"}
    spec = """
    ./HaplotypeCaller.sh {} {}
    """.format(bam_list, output)
    return inputs, outputs, options, spec
def HaplotypeCallerGVCF(bam_file, output, cores=6, time="24:00:00", memory="32g"):
    if ".list" in bam_file:
        inputs = [bam_file]
    else:
        inputs = [bam_file, bam_file+".bai"]
    outputs = [output]
    options = {"cores": cores,
               "memory": memory,
               "account": "NChain",
               "walltime": time}
    spec = """
    ./HaplotypeCallerGVCF.sh {} {} {}
    """.format(bam_file, output, cores)
    return inputs, outputs, options, spec
def make_vcf_list(vcf_files, output):
    inputs = []
    for vcf_file in vcf_files:
        inputs.append(vcf_file)
    outputs = [output]
    options = {"cores": 1,
              "memory": "4g",
              "account": "NChain",
              "walltime": "12:00:00"}
    spec = ""
    for vcf_file in vcf_files:
        spec += "echo $PWD/{} >> {}\n".format(vcf_file, output)
    return inputs, outputs, options, spec
def CombineGVCFs(vcf_files, output):
    inputs = vcf_files
    outputs = [output]
    options = {"cores": 1,
               "memory": "46g",
               "account": "NChain",
               "walltime": "36:00:00"}
    spec = """
source /com/extra/java/8/load.sh
source /com/extra/GATK/3.6/load.sh
java -Xmx64g -jar GATK/version3.8/GenomeAnalysisTK.jar \
-T CombineGVCFs \
-o {} \
-R /home/marnit/NChain/faststorage/WHITE_CLOVER/SofieWhiteClover/references/TrR/TrR.v5.fasta \
""".format(output)
    for vcf_file in vcf_files:
        spec += "--variant {} ".format(vcf_file)
    return inputs, outputs, options, spec
def GenotypeGVCFs(vcf_file, output):
    inputs = [vcf_file]
    outputs = [output]
    options = {"cores": 1,
               "memory": "64g",
               "account": "NChain",
               "walltime": "24:00:00"}
    spec = """
source /com/extra/java/8/load.sh
source /com/extra/GATK/3.6/load.sh
java -Xmx64g -jar GATK/version3.8/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
--variant {} \
-o {} \
-R /home/marnit/NChain/faststorage/WHITE_CLOVER/SofieWhiteClover/references/TrR/TrR.v5.fasta
""".format(vcf_file, output)
    return inputs, outputs, options, spec
def FilterVariants(in_vcf, out_vcf):
    inputs = [in_vcf]
    outputs = [out_vcf]
    options = {"cores": 1,
               "memory": "32g",
               "account": "NChain",
               "walltime": "12:00:00"}
    spec = """
    ./Selectvariantsall.sh {in_vcf} {out_vcf}
    """.format(in_vcf=in_vcf, out_vcf=out_vcf)
    return inputs, outputs, options, spec
############################################################################################################################################################
Reference = "../REFERENCE/TrR.v5.fasta"
gwf.target_from_template("star_index",
                         star_index(Reference))
Read_directories = ["../RAWDATA/5samples/raw_data",
                    "../RAWDATA/140samples/",
                    "../RAWDATA/new41samples/X201SC20030426-Z01-F001_01",
                    "../RAWDATA/new41samples/X201SC20030426-Z01-F001_02",
                    "../RAWDATA/new41samples/X201SC20030426-Z01-F001_03",
                    "../RAWDATA/new41samples/X201SC20030426-Z01-F001"]
read_files = {}
## Collect all of the read pairs in the directory
for rootdir in Read_directories:
    for root, subdirs, files in os.walk(rootdir):
        for filename in files:
            if filename.endswith(".fq.gz"):
                if filename.startswith("WC"):
                    sample_name = filename[:6]
                    if sample_name not in read_files:
                        read_files[sample_name] = list()
                    read_files[sample_name].append(os.path.join(root,filename))
                if filename.startswith("A"):
                    sample_name = filename.split("_")[0]
                    if sample_name not in read_files:
                        read_files[sample_name] = list()
                    read_files[sample_name].append(os.path.join(root,filename))
individuals = list(read_files.keys())
read_files["S9"] = {}
## ADD all S9 reads:
for root, subdirs, files in os.walk("/faststorage/project/NChain/WHITE_CLOVER/20181211_RNASEQ/new_rnaseq/raw_data"):
    for filename in files:
        if filename.endswith(".fq.gz") and filename.startswith("S9"):
            sample_name = "_".join(filename.split("_")[:3])
            if sample_name not in read_files["S9"]:
                read_files["S9"][sample_name] = list()
            read_files["S9"][sample_name].append(os.path.join(root, filename))
#print(read_files)
#print(individuals)
bam_files = []
## Mapping Individuals using STAR
for indv in individuals:
    gwf.target_from_template("STARmapping{indv}".format(indv=indv),
                             star_mapping(read_files[indv][0],
                                          read_files[indv][1],
                                          "MAPPED_READS/{}".format(indv),
                                          Reference))
    bam_files.append("MAPPED_READS/{}".format(indv)+"Aligned.sortedByCoord.out.bam")
S9_tissues = list(read_files["S9"].keys())
#print(S9_tissues)
S9_bam_files = []
## Readmapping for all tissues of the S9 samples
for indv in S9_tissues:
    gwf.target_from_template("STARmapping{indv}".format(indv=indv),
                             star_mapping(read_files["S9"][indv][0],
                                          read_files["S9"][indv][1],
                                          "MAPPED_READS/{}".format(indv),
                                          Reference))
    S9_bam_files.append("MAPPED_READS/{}".format(indv)+"Aligned.sortedByCoord.out.bam")
#print(S9_bam_files)
fixed_bam_files = []
for indv, bam_file in zip(individuals, bam_files):
    gwf.target_from_template("Add_RG{}".format(indv),
                             add_readgroups(bam_file,
                                            "RG_READS/"+indv+".bam",
                                            indv))
    fixed_bam_files.append("RG_READS/"+indv+".bam")
for indv, bam_file in zip(S9_tissues, S9_bam_files):
    gwf.target_from_template("Add_RG{}".format(indv),
                             add_readgroups(bam_file,
                                            "RG_READS/"+indv+".bam",
                                            "S9"))
    fixed_bam_files.append("RG_READS/"+indv+".bam")
#print(fixed_bam_files)
md_bam_files = []
for bam_file in fixed_bam_files:
    gwf.target_from_template("Mark_Dups{}".format(bam_file.split("/")[-1].split(".")[0]),
                             mark_duplicates(bam_file,
                                            "MD_READS/"+bam_file.split("/")[-1].split(".")[0]+".bam"))
    md_bam_files.append("MD_READS/"+bam_file.split("/")[-1].split(".")[0]+".bam")
#print(md_bam_files)
lenfix_bam_files = []
for bam_file in md_bam_files:
    gwf.target_from_template("Fix_Contig{}".format(bam_file.split("/")[-1].split(".")[0]),
                             fix_contig_len(bam_file,
                                            "FIXED_READS/"+bam_file.split("/")[-1].split(".")[0]+".bam"))
    lenfix_bam_files.append("FIXED_READS/"+bam_file.split("/")[-1].split(".")[0]+".bam")
    gwf.target_from_template("Index{}".format(lenfix_bam_files[-1].split("/")[-1].split(".")[0]),
                                                   samtools_index(lenfix_bam_files[-1]))
#print(lenfix_bam_files)
final_bam_files = []
for bam_file in lenfix_bam_files:
    gwf.target_from_template("SplitNTrim{}".format(bam_file.split("/")[-1].split(".")[0]),
                             SplitNTrim(bam_file,
                                            "FINAL_READS/"+bam_file.split("/")[-1].split(".")[0]+".bam"))
    final_bam_files.append("FINAL_READS/"+bam_file.split("/")[-1].split(".")[0]+".bam")
for bam_file in final_bam_files:
    gwf.target_from_template("LastIndex{}".format(bam_file.split("/")[-1].split(".")[0]),
                             samtools_index(bam_file))
gvcf_files = []
## Variant call S9 samples
S9_samples = list(filter(lambda x: "S9" in x, final_bam_files))
gwf.target_from_template("MakeS9BamList",
                         make_bam_list(S9_samples, "S9bam.list"))
gwf.target_from_template("S9variant_calling",
                         HaplotypeCallerGVCF("S9bam.list",
                                             "VCF_FILES/S9.g.vcf",
                                             cores=8,
                                             memory="64g",
                                             time="120:00:00"))
gvcf_files.append("VCF_FILES/S9.g.vcf")
remaining_bam_files = list(filter(lambda x: not "S9" in x, final_bam_files))
## Variant call for each individual sample
for bam_file in remaining_bam_files:
    gwf.target_from_template("{}Variant_calling".format(bam_file.split("/")[-1].split(".")[0]),
                             HaplotypeCallerGVCF(bam_file,
                                                 "VCF_FILES/{}.g.vcf".format(bam_file.split("/")[-1].split(".")[0])))
    gvcf_files.append("VCF_FILES/{}.g.vcf".format(bam_file.split("/")[-1].split(".")[0]))
#gwf.target_from_template("makeVcfList",
#                         make_vcf_list(gvcf_files,
#                                       "vcf_files.txt"))
batch_vcf_files = []
batch_size = 20
for i in range(ceil(len(gvcf_files)/batch_size)):
    current_batch = gvcf_files[(i*batch_size):((i+1)*batch_size)]
    gwf.target_from_template("CombineVCFsBatch{}".format(i),
                             CombineGVCFs(current_batch,
                                          "JOINED_VCFS/batch{}.g.vcf".format(i)))
    batch_vcf_files.append("JOINED_VCFS/batch{}.g.vcf".format(i))
gwf.target_from_template("FinalCombineVCF",
                         CombineGVCFs(batch_vcf_files,
                                      "joined.g.vcf"))
gwf.target_from_template("GenotypeGVCFs",
                         GenotypeGVCFs("joined.g.vcf",
                                       "RNAseq_unfiltered_variants.vcf"))
gwf.target_from_template("GATKfiltervariants",
                         FilterVariants("RNAseq_unfiltered_variants.vcf",
                                        "RNAseq_variants.vcf"))
#gwf.target_from_template("MakeBamList",
#                         make_bam_list(final_bam_files, "bam.list"))
# #gwf.target_from_template("MakeBamListNoTrimming",
# #                         make_bam_list(lenfix_bam_files, "bam.no.trim.list"))
#
#
#gwf.target_from_template("GATKvariantcalling",
#                         HaplotypeCaller("bam.list", "RNAseq_unfiltered_variants.vcf"))
# #gwf.target_from_template("GATKvariantcallingNoTrim",
# #                         HaplotypeCaller("bam.no.trim.list", "white_clover_unfiltered_variants.no.trim.vcf"))
#
#gwf.target_from_template("GATKfiltervariants",
#                         FilterVariants("RNAseq_unfiltered_variants.vcf",
#                                        "RNAseq_variants.vcf"))
```

Add readgroups to BAM files
---------------------------

Using picard to add read groups. Usage: ./picard.sh <in_bam_file> <out_bam_file> <read group>. Used to append read groups to the bam files after read mapping with STAR is done. This step is required for GATK to process the reads.

``` bash
#!/bin/bash

source /com/extra/picard/2.0.1/load.sh
source /com/extra/java/8/load.sh

picard AddOrReplaceReadGroups INPUT=$1 OUTPUT=$2 SO=coordinate RGID=$3 RGLB=library1 RGPL=illumina RGPU=H0164ALXX140820.2 RGSM=$3 
```

Mark duplicates and correct quality
-----------------------------------

Usage: ./markduplicates.sh <in_bam_file> <out_bam_file>. Mark any duplicate reads, however main usage is for the process to correct mapping quality used by STAR and standardizing it for usage with GATK.

``` bash
#!/bin/bash

source /com/extra/picard/2.0.1/load.sh
source /com/extra/java/8/load.sh

picard MarkDuplicates I=$1 O=$2 CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metric 
```

SplitNTrim RNA reads
--------------------

Usage: ./SplitNTrim.sh <in_bam_file> <out_bam_file>. Fixes RNA reads and prepares them for input to GATK.

``` bash
#!/bin/bash

source /com/extra/java/8/load.sh
source /com/extra/GATK/3.6/load.sh

java -Xmx64g -Djava.io.tmpdir=./gtak_tmp -jar /com/extra/GATK/3.6/jar-bin/GenomeAnalysisTK.jar \
     -T SplitNCigarReads \
     -I $1 \
     -o $2 \
     -rf ReassignOneMappingQuality \
     -RMQF 255 \
     -RMQT 60 \
     -R /home/marnit/NChain/faststorage/WHITE_CLOVER/SofieWhiteClover/references/TrR/TrR.v5.fasta \
     -U ALLOW_N_CIGAR_READS \
     2>&1
```

Variant calling (Unifiedgenotyper)
----------------------------------

Usage: ./Unifiedgenotyper.sh <bamlist>
<output vcf file>
. Uses the Unifiedgenotyper from GATK to do variant calling.

``` bash
#!/bin/bash

source /com/extra/java/8/load.sh
source /com/extra/GATK/3.6/load.sh

java -Xmx64g -jar /com/extra/GATK/3.6/jar-bin/GenomeAnalysisTK.jar \
     -T UnifiedGenotyper \
     -I $1 \
     -o $2 \
     -stand_call_conf 50 \
     -stand_emit_conf 20.0 \
     --sample_ploidy 2 \
     -nct 12 --genotyping_mode DISCOVERY \
     -R /home/marnit/NChain/faststorage/WHITE_CLOVER/SofieWhiteClover/references/TrR/TrR.v5.fasta \
     --output_mode EMIT_ALL_CONFIDENT_SITES \
     -rf BadCigar 2>&1 > log
```

Variant calling (HaplotypeCaller)
---------------------------------

Usage: ./HaplotypeCaller.sh <bamlist>
<output vcf file>
. Uses the HaplotypeCaller from GATK to do variant calling (this is recommended caller when using RNAseq data).

``` bash
#!/bin/bash

source /com/extra/java/8/load.sh
source /com/extra/GATK/3.6/load.sh

java -Xmx64g -jar GATK/version3.8/GenomeAnalysisTK.jar \
     -T HaplotypeCaller \
     -I $1 \
     -o $2 \
     -stand_call_conf 20.0 \
     -dontUseSoftClippedBases \
     -nct 16 \
     -R /home/marnit/NChain/faststorage/WHITE_CLOVER/SofieWhiteClover/references/TrR/TrR.v5.fasta \
     2>&1 > log
```

Variant calling (HaplotypeCallerGVCF)
-------------------------------------

Usage: ./HaplotypeCallerGVCF.sh <bamlist>
<output vcf file>
. Uses the HaplotypeCaller from GATK to do variant calling, but produces a gvcf file instead which outputs all callable sites. This is useful when performing variant calling on each individual separately when parallelizing the workflow. In the case of more than 80+ individuals, this is required due to the resource requirements when doing variant calling on all of the individuals at the same time.

``` bash
#!/bin/bash

source /com/extra/java/8/load.sh
source /com/extra/GATK/3.6/load.sh

java -Xmx64g -jar GATK/version3.8/GenomeAnalysisTK.jar \
     -T HaplotypeCaller \
     -I $1 \
     -o $2 \
     -stand_call_conf 20.0 \
     -dontUseSoftClippedBases \
     --emitRefConfidence GVCF \
     -nct $3 \
     -R /home/marnit/NChain/faststorage/WHITE_CLOVER/SofieWhiteClover/references/TrR/TrR.v5.fasta \
     2>&1 > log
```

FilterVariants
--------------

Usage: ./Selectvariantsall.sh <bamlist>
<output vcf file>
. Filter variants using the Selectvariants command in GATK. Settings are limited to basic quality checks, and a depth filter based on the amount of individuals in the variant calling.

``` bash
#!/bin/bash

source /com/extra/java/8/load.sh

java -Xmx64g -jar GATK/version3.8/GenomeAnalysisTK.jar \
     -R /home/marnit/NChain/faststorage/WHITE_CLOVER/SofieWhiteClover/references/TrR/TrR.v5.fasta \
     -T SelectVariants \
     -o $2 \
     --variant $1 \
     --restrictAllelesTo ALL -select "MQ>30.00 && DP>160 && QUAL>20.00"
```
