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
