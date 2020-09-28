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
