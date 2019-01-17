#!/bin/bash

GTF=$1
BED=$1.bed.gz

zcat -f $GTF | sed 's/[\"\;]//g' | awk 'BEGIN{OFS="\t"} {print $1,$4,$5,$10,$6,$7}' | gzip -nc > $BED
