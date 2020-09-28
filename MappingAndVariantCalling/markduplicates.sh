#!/bin/bash

source /com/extra/picard/2.0.1/load.sh
source /com/extra/java/8/load.sh

picard MarkDuplicates I=$1 O=$2 CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metric 
