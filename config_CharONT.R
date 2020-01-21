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

####################################################################################################
##Note: rows starting with '#' are notes for the user, and are ignored by the software
#BC_int <- c("BC01", "BC02", "BC03", "BC04", "BC05", "BC06", "BC07", "BC08", "BC09", "BC10", "BC11", "BC12")
BC_int <- c("BC01", "BC02", "BC03", "BC04", "BC05", "BC06", "BC07", "BC08", "BC09", "BC10", "BC11", "BC12")
#if do_subsampling_flag <- 1, subsampling of num_fast5_files fast5 files is performed; otherwise set do_subsampling_flag <- 0
do_subsampling_flag <- 0
#num_fast5_files is the number of fast5 files to be subsampled/analysed (if do_subsampling_flag <- 1)
num_fast5_files <- 25
#kit (1D/1D^2 reads/rapid 16S)
kit <- "SQK-LSK109"
#flowcell chemistry (R9.4/R9.5 chemistry)
flowcell <- "FLO-MIN106"
#fast_basecalling_flag <- 1 if you want to use the fast basecalling algorithm; otherwise set fast_basecalling_flag <- 0 if you want to use the accurate but slow one (FLO-MIN106 only)
fast_basecalling_flag <- 1
#pair_strands_flag <- 1 if, in case a 1d2 kit and FLO-MIN107 flow-cell have been used, you want to perform 1d2 basecalling; otherwise set pair_strands_flag <- 0
pair_strands_flag <- 0
#set the maximum number of threads to be used
num_threads <- 8
#set length [bp] of PCR primers to be trimmed at both sides of consensus sequences
primers_length <- 25
#do_in_silico_pcr <- 1 in case you are not expecting PCR-like amplicons and want to extract on-target trimmed reads; otherwise, set do_insilico_pcr <- 0
do_in_silico_pcr <- 0
#if do_in_silico_pcr <- 1, extract portion of reads between sequences pcr_silico_primer_one and pcr_silico_primer_two, pcr_silico_primer sequences included
pcr_silico_primer_one <- "sequence_of_interest"
pcr_silico_primer_two <- "sequence_of_interest"
#if skip_demultiplexing_flag <- 1 demultiplexing is skipped; otherwise set skip_demultiplexing_flag <- 0
skip_demultiplexing_flag <- 0
#require_two_barcodes_flag <- 1 if you want to keep only reads with a barcode (tag) at both ends of the read; otherwise set require_two_barcodes_flag <- 0
require_two_barcodes_flag <- 0
#min read quality value
min_qual <- 7
#minimum minor allele frequency; if less than min_maf*100% of reads are assigned to Allele #2, the sample is assumed homozygous
min_maf <- 0.2
#set haploid_flag <- 1 if you are studying haploid chromsomes (e.g. sexual chrosomomes in man); otherwise set haploid_flag <- 0
haploid_flag <- 0
#label as outliers reads with score > 3rd_QR + IQR_outliers_coef*IQR or score < 1st_QR - IQR_outliers_coef*IQR
IQR_outliers_coef <- 2
#set fast_alignment_flag <- 1 if you want to perform fast multiple sequence alignment; otherwise set fast_alignment_flag <- 0
fast_alignment_flag <- 1
########################################################################################################
PIPELINE_DIR <- "/path/to/CharONT"
#MINICONDA DIR
MINICONDA_DIR <- "/path/to/miniconda3"
#basecaller_dir
BASECALLER_DIR <- "/path/to/ont-guppy-cpu/bin/"
########### End of user editable region ################################################################
#load BioStrings package
suppressMessages(library(Biostrings))
#load stats package
suppressMessages(library(stats))
#path to CharONT.R
CharONT <- paste0(PIPELINE_DIR, "/CharONT.R")
#path to DecONT.sh
DECONT <- paste0(PIPELINE_DIR, "/decONT.sh")
#path to subsample fast5
subsample_fast5 <- paste0(PIPELINE_DIR, "/subsample_fast5.sh")
#MAFFT
MAFFT <- paste0(MINICONDA_DIR, "/envs/CharONT_env/bin/mafft")
#VSEARCH
VSEARCH <- paste0(MINICONDA_DIR, "/envs/CharONT_env/bin/vsearch")
#NANOPOLISH
NANOPOLISH <- paste0(MINICONDA_DIR, "/envs/CharONT_env/bin/nanopolish")
#EMBOSS cons
CONS <- paste0(MINICONDA_DIR, "/envs/CharONT_env/bin/cons")
#SEQTK
SEQTK <- paste0(MINICONDA_DIR, "/envs/CharONT_env/bin/seqtk")
#MINIMAP2
MINIMAP2 <- paste0(MINICONDA_DIR, "/envs/CharONT_env/bin/minimap2")
#SAMTOOLS
SAMTOOLS <- paste0(MINICONDA_DIR, "/envs/CharONT_env/bin/samtools")
#TRF
TRF <- paste0(MINICONDA_DIR, "/envs/CharONT_env/bin/trf")
#PYCOQC
PYCOQC <- paste0(MINICONDA_DIR, "/envs/CharONT_env/bin/pycoQC")
#NANOFILT
NANOFILT <- paste0(MINICONDA_DIR, "/envs/CharONT_env/bin/NanoFilt")
#MSA
MSA <- paste0(MINICONDA_DIR, "/envs/CharONT_env/bin/msa.sh")
#CUTPRIMERS
CUTPRIMERS <- paste0(MINICONDA_DIR, "/envs/CharONT_env/bin/cutprimers.sh")