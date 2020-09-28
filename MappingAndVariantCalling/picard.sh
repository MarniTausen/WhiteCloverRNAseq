#!/bin/bash

source /com/extra/picard/2.0.1/load.sh
source /com/extra/java/8/load.sh

picard AddOrReplaceReadGroups INPUT=$1 OUTPUT=$2 SO=coordinate RGID=$3 RGLB=library1 RGPL=illumina RGPU=H0164ALXX140820.2 RGSM=$3 
