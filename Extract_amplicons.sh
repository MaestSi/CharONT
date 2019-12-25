#!/bin/bash

#
# Copyright 2019 Simone Maestri. All rights reserved.
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

FASTQ_READS=$1
PRIMER_SEQ_ONE=$2
PRIMER_SEQ_TWO=$3

MSA=msa.sh
CUTPRIMERS=cutprimers.sh
SEQTK=seqtk

working_dir=$(dirname $(realpath $FASTQ_READS))
reads_full=$(realpath $FASTQ_READS)
sample_name=$(echo $(basename $(realpath $FASTQ_READS)) | sed 's/\.fa.*//' | sed 's/\.fq.*//')
sam_file_one=$working_dir"/"$sample_name"_in_silico_pcr_one.sam"
sam_file_two=$working_dir"/"$sample_name"_in_silico_pcr_two.sam"
trimmed_reads_fq=$working_dir"/BC01.fastq"
trimmed_reads_fa=$working_dir"/BC01.fasta"

$MSA in=$reads_full out=$sam_file_one literal=$PRIMER_SEQ_ONE qin=33 cutoff=0.8
$MSA in=$reads_full out=$sam_file_two literal=$PRIMER_SEQ_TWO qin=33 cutoff=0.8
$CUTPRIMERS in=$reads_full out=$trimmed_reads_fq sam1=$sam_file_one sam2=$sam_file_two qin=33 fake=f include=t fixjunk
$SEQTK seq -A $trimmed_reads_fq > $trimmed_reads_fa
