#!/bin/bash

source /com/extra/java/8/load.sh

java -Xmx64g -jar GATK/version3.8/GenomeAnalysisTK.jar \
     -R /home/marnit/NChain/faststorage/WHITE_CLOVER/SofieWhiteClover/references/TrR/TrR.v5.fasta \
     -T SelectVariants \
     -o $2 \
     --variant $1 \
     --restrictAllelesTo ALL -select "MQ>30.00 && DP>160 && QUAL>20.00"
