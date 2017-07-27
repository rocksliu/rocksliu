#!/bin/bash
bamfile=$1
bedfile=$2
sample=$3
outpath=$4

samtools view -H $bamfile > $outpath/$sample.Intervals
awk 'OFS="\t"{$4="+";$5=$1":"$2"-"$3;print $1,$2,$3,$4,$5}' $bedfile >> $outpath/$sample.Intervals
picard CalculateHsMetrics I=$bamfile o=$outpath/$sample.CalculateHSmetrics \
 R=/path/to/ucsc.hg19.fasta TI=$outpath/$sample.Intervals \
 BI=$outpath/$sample.Intervals VALIDATION_STRINGENCY=SILENT PER_TARGET_COVERAGE=$outpath/$sample.coverage.bed
