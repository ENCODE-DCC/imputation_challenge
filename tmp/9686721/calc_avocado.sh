#!/bin/bash

GIT_REPO=~/code/imputation_challenge
GENE_ANNOT=$GIT_REPO/annot/hg38/gencode.v29.genes.gtf.bed.gz
ENH_ANNOT=$GIT_REPO/annot/hg38/F5.hg38.enhancers.bed.gz

DIR_BW=/mnt/imputation-challenge/data/evaluation_data/avocado
DIR_BW_TRUTH=/mnt/imputation-challenge/data/validation_data
DIR_OUT=/mnt/imputation-challenge/output/score/avocado
DIR_LOG=/mnt/imputation-challenge/output/score/avocado/logs

ALL_CHRS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX"
mkdir -p $DIR_OUT $DIR_LOG

BIGWIGS=(C02M22
C03M02
C04M16
C09M20
C10M17
C12M16
C12M32
C13M20
C16M17
C17M04
C17M19
C17M29
C17M32
C18M21
C18M25
C20M22
C23M03
C23M07
C23M26
C23M34
C24M17
C24M25
C25M21
C25M26
C27M03
C27M13
C27M24
C27M26
C29M29
C31M25
C32M08
C32M12
C32M20
C34M02
C34M32
C36M18
C37M29
C45M22
C46M10
C46M18
C46M21
C46M35
C47M18
C48M16
C50M02)

IDS=(p1
p2
p3
p4
p5
p6
p7
p8
p9
p10)

N=30
(
for B in "${BIGWIGS[@]}"; do  
  BW_TRUTH=$DIR_BW_TRUTH/$B.bigwig
  for IDX in "${IDS[@]}"; do  
    ((i=i%N)); ((i++==0)) && wait    
    BW=$DIR_BW/$B.$IDX.bigwig
    OUT=$DIR_OUT/$B.$IDX.score.txt
    STDOUT=$DIR_LOG/$B.$IDX.o
    STDERR=$DIR_LOG/$B.$IDX.e
    echo "Scoring $BW, O=$STDOUT, E=$STDERR"
    python $GIT_REPO/score.py $BW $BW_TRUTH \
      --gene-annotations $GENE_ANNOT \
      --enh-annotations $ENH_ANNOT \
      --chrom $ALL_CHRS \
      --window-size 25 --prom-loc 80 \
      --out $OUT 1> $STDOUT 2> $STDERR &
  done
done
)



