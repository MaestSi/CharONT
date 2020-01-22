# CharONT

**CharONT** is a consensus calling pipeline meant for characterizing long genomic regions from diploid organisms. Starting from ONT reads including a shared flanking sequence, it provides consensus sequences for the two alleles and tandem repeats annotations. In case you used an enrichment method different to PCR, amplicons can be extracted _in-silico_ based on known flanking sequences. Moreover, a preprocessing pipeline is provided, so to make the whole bioinformatic analysis from raw fast5 files to consensus sequences straightforward and simple.

<p align="center">
  <img src="Figures/CharONT_logo.png" alt="drawing" width=450" title="CharONT_logo">
</p>

## Getting started

**Prerequisites**

* Miniconda3.
Tested with conda 4.8.1.
```which conda``` should return the path to the executable.
If you don't have Miniconda3 installed, you could download and install it with:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod 755 Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
```

Then, after completing _CharONT_ installation, set the _MINICONDA_DIR_ variable in **config_CharONT.R** to the full path to miniconda3 directory.

* Guppy, the software for basecalling and demultiplexing provided by ONT. Tested with Guppy v3.4.4.
If you don't have [Guppy](https://community.nanoporetech.com/downloads) installed, choose an appropriate version and install it.
For example, you could download and unpack the archive with:
```
wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_version_of_interest.tar.gz
tar -xf ont-guppy-cpu_version_of_interest.tar.gz
```
A directory _ont-guppy-cpu_ should have been created in your current directory.
Then, after completing _CharONT_ installation, set the _BASECALLER_DIR_ variable in **config_CharONT.R** to the full path to _ont-guppy-cpu/bin_ directory.

* gcc ver. >= 4.8

* cmake ver. >= 3.2


**Installation**

```
git clone https://github.com/MaestSi/CharONT.git
cd CharONT
chmod 755 *
./install.sh
```

A conda environment named _CharONT_env_ is created, where emboss, vsearch, seqtk, mafft, minimap2, samtools, racon, medaka, NanoFilt, Tandem Repeat Finder, BBMap, pycoQC and R with package Biostrings are installed.
Then, you can open the **config_CharONT.R** file with a text editor and set the variables _PIPELINE_DIR_ and _MINICONDA_DIR_ to the value suggested by the installation step.

## Overview

<p align="center">
  <img src="Figures/CharONT_pipeline_flowchart.png" alt="drawing" width="700" title="CharONT_pipeline_flowchart">
</p>

## Usage

The CharONT pipeline can be applied either starting from raw fast5 files, or from already basecalled and demultiplexed sequences. 
In both cases, the first step of the pipeline requires you to open the **config_CharONT.R** file with a text editor and to modify it according to the features of your sequencing experiment and your preferences.
If you have already basecalled and demultiplexed your sequences, you can run the pipeline using the **CharONT.R** script.
Otherwise, you can run the pipeline using the **Launch_CharONT.sh** script. If you have basecalled sequences that need filtering and trimming to look like PCR amplicons, the script **Extract_amplicons.sh** can be applied prior to **CharONT.R**.

**CharONT.R**

Usage: Rscript CharONT.R \<analysis_dir\>

Note: Activate the virtual environment with ```source activate CharONT_env``` before running. The script is run by **CharONT_preprocessing.R**, but can be also run as a main script if you have already basecalled and demultiplexed your sequences.

Inputs:
* \<analysis_dir\>: directory containing fastq files for each sample named BC\<numbers\>.fastq

Outputs (saved in <analysis_dir>):
* \<"sample_name"\_first_allele.fasta\>: consensus sequence for first allele in fasta format
* \<"sample_name"\_first_allele.fasta."trf scores".html\>: Tandem Repeat Finder report for first allele sequence
* \<"sample_name"\_second_allele.fasta\>: consensus sequence for second allele in fasta format
* \<"sample_name"\_second_allele.fasta."trf scores".html\>: Tandem Repeat Finder report for second allele sequence
* \<"sample_name"\>: directory including intermediate files

**Launch_CharONT.sh**

Usage:
Launch_CharONT.sh \<fast5_dir\>

Note: modify **config_CharONT.R** before running; the script runs the full pipeline from raw fast5 files to consensus sequences.

Input
* \<fast5_dir\>: directory containing raw fast5 files

Outputs (saved in \<fast5_dir\>_analysis/analysis):

* \<"sample_name"\_first_allele.fasta\>: consensus sequence for first allele in fasta format
* \<"sample_name"\_first_allele.fasta."trf scores".html\>: Tandem Repeat Finder report for first allele sequence
* \<"sample_name"\_second_allele.fasta\>: consensus sequence for second allele in fasta format
* \<"sample_name"\_second_allele.fasta."trf scores".html\>: Tandem Repeat Finder report for second allele sequence
* \<"sample_name"\>: directory including intermediate files

Outputs (saved in \<fast5_dir\>_analysis/qc):
* Read length distributions and pycoQC report

Outputs (saved in \<fast5_dir\>_analysis/basecalling):
* Temporary files for basecalling

Outputs (saved in \<fast5_dir\>_analysis/preprocessing):
* Temporary files for demultiplexing, filtering based on read length and adapters trimming

**Extract_amplicons.sh**

Usage: Extract_amplicons.sh \<fastq_reads\> \<primer_sequence_one\> \<primer_sequence_two\>

Note: Activate the virtual environment with ```source activate CharONT_env``` before running.

Inputs:
* \<fastq_reads\>: fastq file containing reads that need filtering and trimming to look like PCR amplicons
* \<primer_sequence_one\>: sequence of first _in-silico_ PCR primer to look for, flanking the region of interest
* \<primer_sequence_two\>: sequence of second _in-silico_ PCR primer to look for, flanking the region of interest

Outputs:
* BC01.fast*: fastq and fasta files containing amplicon-like sequences extracted from \<fastq_reads\> based on \<primer_sequence_one\> and \<primer_sequence_two\> sequences
* \<"sample_name"\_in_silico_pcr_one.sam\>: sam file containing alignments between \<primer_sequence_one\> and \<fastq_reads\>
* \<"sample_name"\_in_silico_pcr_two.sam\>: sam file containing alignments between \<primer_sequence_two\> and \<fastq_reads\>
 

## Auxiliary scripts

In the following, auxiliary scripts run either by **CharONT.R** or by **Launch_CharONT.sh** are listed. These scripts should not be called directly.

**CharONT_preprocessing.R**

Note: script run by _Launch_CharONT.sh_.

**config_CharONT.R**

Note: configuration script, must be modified before running _Launch_CharONT.sh_ or _CharONT.R_.

**subsample_fast5.sh**

Note: script run by _CharONT_preprocessing.R_ if _do_subsampling_flag_ variable is set to 1 in _config_CharONT.R_.
