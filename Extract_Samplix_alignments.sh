#!/bin/bash

#
# Copyright 2020 Simone Maestri. All rights reserved.
# Simone Maestri <simone.maestri@univr.it>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

READS=$1
REFERENCE=$2
BED=$3

PCR_CLIP_READS="java -jar /path/to/jvarkit/dist/pcrclipreads.jar" #to be installed
SAM_EXTRACT_CLIP="java -jar /path/to/jvarkit/dist/samextractclip.jar" #to be installed
EXTRACT_AMPLICONS="/path/to/Extract_amplicons.sh"
SEQTK=seqtk
MINIMAP2=minimap2
SAMTOOLS=samtools
BEDTOOLS=bedtools

CHR=$(cat $BED | cut -f1)
FW_START=$(cat $BED | cut -f2)
FW_END=$(echo $FW_START + 20 | bc)
RV_END=$(cat $BED | cut -f3)
RV_START=$(echo $RV_END - 20 | bc)
FW_PRIMER=$($SAMTOOLS faidx $REFERENCE $CHR":"$FW_START"-"$FW_END | grep -v ">")
RV_PRIMER=$($SAMTOOLS faidx $REFERENCE $CHR":"$RV_START"-"$RV_END | grep -v ">")

SAMPLE_NAME=$(echo $(basename $READS) | sed 's/\.fast.*//')
BED_NAME=$(echo $(basename $BED) | sed 's/\.bed//')
WORKING_DIR=$(realpath $(dirname $READS))
BAM=$WORKING_DIR"/"$SAMPLE_NAME".bam"
BAM_target=$WORKING_DIR"/"$SAMPLE_NAME"_"$BED_NAME".bam"
FASTQ_forward=$WORKING_DIR"/"$SAMPLE_NAME"_"$BED_NAME"_fw.fastq"
BAM_target_forward=$WORKING_DIR"/"$SAMPLE_NAME"_"$BED_NAME"_fw.bam"
SAM_target_clipped=$WORKING_DIR"/"$SAMPLE_NAME"_"$BED_NAME"_fw_clipped.sam"
FASTQ_target_clipped=$WORKING_DIR"/"$SAMPLE_NAME"_"$BED_NAME"_fw_clipped.fastq"
FASTQ_target_clipped_ids=$WORKING_DIR"/"$SAMPLE_NAME"_"$BED_NAME"_fw_clipped_only_reads_ids.txt"
FASTQ_extracted=$WORKING_DIR"/"$SAMPLE_NAME"_"$BED_NAME"_extracted.fastq"

#map reads to reference
$MINIMAP2 -ax map-ont --MD $REFERENCE $READS | $SAMTOOLS view -h | $SAMTOOLS sort -o $BAM -T reads.tmp
$SAMTOOLS index $BAM
#extract reads overlapping to region of interest
$BEDTOOLS intersect -a $BAM -b $BED -F 0.9 > $BAM_target
#extract fastq sequence oriented in forward direction
$SAMTOOLS view $BAM_target -F260 | awk 'BEGIN {FS="\t"} {print "@" $1"_"$2"_"$3"_"$4"_"$12 "\n" $10 "\n+\n" $11}' > $FASTQ_forward
#map forward-oriented reads to reference
$MINIMAP2 -ax map-ont --MD $REFERENCE $FASTQ_forward | $SAMTOOLS view -hSb -F2068 | $BEDTOOLS intersect -a stdin -b $BED -F 0.9 | $SAMTOOLS sort -o $BAM_target_forward -T reads.tmp
#clip portions of reads outside of bed
$PCR_CLIP_READS -B $BED $BAM_target_forward | $SAMTOOLS view -q 1 -hF 4 > $SAM_target_clipped
#extract clipped and not clipped portions of reads
$SAM_EXTRACT_CLIP -c $SAM_target_clipped -o $FASTQ_target_clipped
#extract the id of the portions of interest of reads
cat $FASTQ_target_clipped | grep "clipped" | sed 's/^@//g' > $FASTQ_target_clipped_ids
#extract the portions of interest of reads
$SEQTK subseq $FASTQ_target_clipped $FASTQ_target_clipped_ids |  $SEQTK rename - $SAMPLE_NAME"_" > $FASTQ_extracted
#extract complete reads
$EXTRACT_AMPLICONS $FASTQ_extracted $FW_PRIMER $RV_PRIMER 0.7
mkdir tmp
mv $BAM_target $FASTQ_forward $BAM_target_forward $SAM_target_clipped $FASTQ_target_clipped $FASTQ_target_clipped_ids $FASTQ_extracted "*in_silico_pcr_*" tmp
