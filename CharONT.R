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

args = commandArgs(trailingOnly=TRUE)

if (args[1] == "-h" | args[1] == "--help") {
  cat("", sep = "\n")
  cat(paste0("Usage: Rscript CharONT.R <home_dir> <fast5_dir> <sequencing_summary.txt>"), sep = "\n")
  cat(paste0("Note that config_CharONT.R must be in the same directory of CharONT.R"), sep = "\n")
  cat(paste0("<home_dir>: directory containing fastq and fasta files for each sample"), sep = "\n")
  cat(paste0("<fast5_dir>: directory containing raw fast5 files for nanopolish polishing, optional"), sep = "\n")
  cat(paste0("<sequencing_summary.txt>: sequencing summary file generated during base-calling, used to speed-up polishing, optional"), sep = "\n")
  stop(simpleError(sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "))))
}

if (length(args) == 1) {
  home_dir <- args[1]
  if (!dir.exists(home_dir)) {
    stop(paste0(home_dir, " directory does not exist!"))
  }
} else if (length(args) == 2) {
  home_dir <- args[1]
  fast5_dir <- args[2]
  if (!dir.exists(home_dir)) {
    stop(paste0(home_dir, " directory does not exist!"))
  } else if (!dir.exists(fast5_dir)) {
    stop(paste0(fast5_dir, " directory does not exist!"))
  }
} else if (length(args) == 3) {
  home_dir <- args[1]
  fast5_dir <- args[2]
  sequencing_summary <- args[3]
  if (!dir.exists(home_dir)) {
    stop(paste0(home_dir, " directory does not exist!"))
  } else if (!dir.exists(fast5_dir)) {
    stop(paste0(fast5_dir, " directory does not exist!"))
  }
  
} else {
  stop("At least one input argument must be provided")
}

PIPELINE_DIR <- dirname(strsplit(commandArgs(trailingOnly = FALSE)[4],"=")[[1]][2])

CONFIG_FILE <- paste0(PIPELINE_DIR, "/config_CharONT.R")
source(CONFIG_FILE)

if (!exists("sequencing_summary")) {
  seq_sum_flag <- 0
} else {
  seq_sum_flag <- 1
}

if (!exists("num_threads")) {
  num_threads <- 8
}

if (!exists("haploid_flag")) {
  haploid_flag <- 0
}

if (!exists("pair_strands_flag")) {
  pair_strands_flag <- 0
}

if (!exists("primers_length")) {
  primers_length <- 25
}

if (!exists("min_maf")) {
  min_maf <- 0.2
}

if (!exists("IQR_outliers_coef")) {
  IQR_outliers_coef <- 2
}

if (!exists("fast_alignment_flag")) {
  fast_alignment_flag <- 1
}


#target reads for creating consensus
TRC <- 200
#target reads for polishing
TRP <- 200

logfile <- paste0(home_dir, "/logfile.txt")

fasta_files <- list.files(path = home_dir, pattern = "BC\\d+\\.fasta", full.names = TRUE)
fastq_files <- paste0(home_dir, "/", gsub(pattern = "\\.fasta$", replacement = "\\.fastq", x = basename(fasta_files)))

if (length(fasta_files) > 0) {
  cat(text = paste0("Processing fasta files ", paste0(basename(fasta_files), collapse = ", ")), sep = "\n")
} else {
  stop(paste0("No fasta files in directory ", home_dir))
}

target_reads_contig <- TRC
target_reads_polishing <- TRP
THR <- 0.85
plurality_value <- 0.15*target_reads_contig

#cycle over fasta files
for (i in 1:length(fasta_files)) {
  skip_first_allele_flag <- 0
  skip_second_allele_flag <- 0
  sample_dir <- gsub(pattern = "\\.fasta", replacement = "", x = fasta_files[i])
  sample_name <- basename(sample_dir)
  dir.create(sample_dir)
  #perform clustering and create a preliminary version for Allele #1
  decont_fa_first_preliminary <- paste0(sample_dir, "/", sample_name, "_decont.fasta")
  decont_fq_first_preliminary <- paste0(sample_dir, "/", sample_name, "_decont.fastq")
  system(paste0(DECONT, " ", fasta_files[i], " ", VSEARCH, " ", SEQTK, " ", THR))
  system(paste0("mv ", home_dir, "/decontam_tmp_", sample_name, " ", sample_dir))
  system(paste0("mv ", home_dir, "/", sample_name, "_decont.fasta ", sample_dir))
  system(paste0("mv ", home_dir, "/", sample_name, "_decont.fastq ", sample_dir))
  num_reads_sample <- as.double(system(paste0("cat ", fasta_files[i], " | grep \"^>\" | wc -l"), intern=TRUE))
  num_reads_mac_first_preliminary <- as.double(system(paste0("cat ", decont_fa_first_preliminary, " | grep \"^>\" | wc -l"), intern=TRUE))
  target_reads_contig <- TRC
  target_reads_polishing <- TRP
  first_allele_preliminary <- paste0(home_dir, "/", sample_name, "/", sample_name, "_preliminary.first.contig.fasta")
  if (num_reads_mac_first_preliminary < 3) {
    system(paste0("head -n2 ", home_dir, "/", sample_name, "/decontam_tmp_", sample_name, "/consensus_", sample_name, ".fasta > ", first_allele_preliminary))
  } else {
    if (num_reads_mac_first_preliminary < target_reads_contig) {
      target_reads_contig <- num_reads_mac_first_preliminary
      target_reads_polishing <- num_reads_mac_first_preliminary
    }
    plurality_value <- 0.15*target_reads_contig
    first_allele_preliminary_tmp1 <- paste0(sample_dir, "/", sample_name, "_preliminary.first.contig_tmp1.fasta")
    first_allele_preliminary_tmp2 <- paste0(sample_dir, "/", sample_name, "_preliminary.first.contig_tmp2.fasta")
    first_allele_preliminary_first <- paste0(sample_dir, "/", sample_name, "_preliminary.first.contig.fasta")
    sequences <- readDNAStringSet(fasta_files[i], "fasta")
    ws <- width(sequences)
    amplicon_length <- ceiling(mean(ws))
    draft_reads_fq_first_preliminary <- paste0(sample_dir, "/", sample_name, "_draft_", target_reads_contig, "_reads_first_preliminary.fastq")
    draft_reads_fa_first_preliminary <- paste0(sample_dir, "/", sample_name, "_draft_", target_reads_contig, "_reads_first_preliminary.fasta")
    seed <- 1
    system(paste0(SEQTK, " sample -s ", seed , " ", decont_fq_first_preliminary, " ",  target_reads_contig, " > ", draft_reads_fq_first_preliminary))
    system(paste0(SEQTK, " seq -A ", draft_reads_fq_first_preliminary, " > ", draft_reads_fa_first_preliminary))
    mfa_file_first_preliminary <- gsub(pattern = "\\.fasta$", replacement = ".mfa", x = draft_reads_fa_first_preliminary)
    system(paste0(MAFFT ," --auto --thread ", num_threads, " --adjustdirectionaccurately ", draft_reads_fa_first_preliminary, " > ", mfa_file_first_preliminary))
    system(paste0(CONS, " -sequence ", mfa_file_first_preliminary, " -plurality ", plurality_value, " -outseq ", first_allele_preliminary_tmp1))
    system(paste0("sed 's/[nN]//g' ", first_allele_preliminary_tmp1, " > ", first_allele_preliminary_tmp2))
    DNAStringSet_obj <- readDNAStringSet(first_allele_preliminary_tmp2, "fasta")
    DNAStringSet_obj_renamed <- DNAStringSet_obj
    original_headers <- names(DNAStringSet_obj)
    sequences <- seq(DNAStringSet_obj)
    names(DNAStringSet_obj_renamed) <- "Allele_number_1_preliminary"
    writeXStringSet(x = DNAStringSet_obj_renamed, filepath = first_allele_preliminary, format = "fasta", width = 20000)
  }
  #map all reads to the preliminary Allele #1 in order to identify reads coming from the two alleles
  sam_file_reads_to_first_allele <- paste0(sample_dir, "/", sample_name, "_reads_to_first_allele_preliminary.sam")
  system(paste0(MINIMAP2, " -ax map-ont ", first_allele_preliminary, " ", fastq_files[i], " | ", SAMTOOLS, " view -h -F 2308 -o " , sam_file_reads_to_first_allele))
  sam_reads_to_first_allele_tmp <- read.table(file = sam_file_reads_to_first_allele, fill = TRUE, comment.char = "@", stringsAsFactors = FALSE, sep = "\t", quote = "")
  ind_rm <- grep(pattern = "NM", x =  sam_reads_to_first_allele_tmp[, 1])
  if (length(ind_rm) > 0) {
    sam_reads_to_first_allele <- sam_reads_to_first_allele_tmp[-ind_rm, ]
  } else {
    sam_reads_to_first_allele <- sam_reads_to_first_allele_tmp
  }
  #extract reads names and cigar strings
  reads_names <- sam_reads_to_first_allele[, 1]
  cigar_strings <- sam_reads_to_first_allele[, 6]
  score <- c()
  del_block_lengths <- c()
  max_del_block_lengths <- c()
  ins_block_lengths <- c()
  max_ins_block_lengths <- c()
  clipping_block_lengths <- c()
  max_clipping_block_lengths <- c()
  min_clipped_len <- 50
  
  #cycle over reads and extract for each read the longest portion non matching the reference
  for (k in 1:length(cigar_strings)) {
    cigar_string <- cigar_strings[k]
    nonref_symbols_coords <- gregexpr(pattern= "[DIS]", text = cigar_string)[[1]]
    ref_symbols_coords <- gregexpr(pattern= "M", text = cigar_string)[[1]]
    num_starting_coords <- gregexpr(pattern = "[0-9]+", text = cigar_string)[[1]]
    nonref_starting_coords <- gregexpr(pattern = "[0-9]+[DIS]+", text = cigar_string)[[1]]
    del_starting_coords <- gregexpr(pattern= "[0-9]+D", text = cigar_string)[[1]]
    ins_starting_coords <- gregexpr(pattern= "[0-9]+I", text = cigar_string)[[1]]
    clipping_starting_coords <- gregexpr(pattern = "[0-9]+[S]+", text = cigar_string)[[1]]
    del_block_lengths_curr_read <- c()
    ins_block_lengths_curr_read <- c()
    clipping_block_lengths_curr_read <- c()
    if (nonref_starting_coords[1] != -1) {
      if (del_starting_coords[1] != -1) {
        for (j in 1:length(del_starting_coords)) {
          del_block_curr <- as.numeric(substr(x = cigar_string, start = (del_starting_coords[j]), stop = (del_starting_coords[j] + attr(del_starting_coords, "match.length")[j] - 2)))
          del_block_lengths_curr_read <- c(del_block_lengths_curr_read, del_block_curr)
        }
        del_block_lengths <- c(del_block_lengths, sum(del_block_lengths_curr_read))
        max_del_block_lengths <- c(max_del_block_lengths, max(del_block_lengths_curr_read))
      } else {
        del_block_lengths <- c(del_block_lengths, 0)
        max_del_block_lengths <- c(max_del_block_lengths, 0)
      }
      if (ins_starting_coords[1] != -1) {
        for (m in 1:length(ins_starting_coords)) {
          ins_block_curr <- as.numeric(substr(x = cigar_string, start = (ins_starting_coords[m]), stop = (ins_starting_coords[m] + attr(ins_starting_coords, "match.length")[m] - 2)))
          ins_block_lengths_curr_read <- c(ins_block_lengths_curr_read, ins_block_curr)
        }
        ins_block_lengths <- c(ins_block_lengths, sum(ins_block_lengths_curr_read))
        max_ins_block_lengths <- c(max_ins_block_lengths, max(ins_block_lengths_curr_read))
      } else {
        ins_block_lengths <- c(ins_block_lengths, 0)
        max_ins_block_lengths <- c(max_ins_block_lengths, 0)
      }
      if (clipping_starting_coords[1] != -1) {
        for (l in 1:length(clipping_starting_coords)) {
          clipping_block_curr <- as.numeric(substr(x = cigar_string, start = (clipping_starting_coords[l]), stop = (clipping_starting_coords[l] + attr(clipping_starting_coords, "match.length")[l] - 2)))
          if (clipping_block_curr < min_clipped_len) {
            clipping_block_lengths_curr_read <- c(clipping_block_lengths_curr_read, 0)
          } else {
            clipping_block_lengths_curr_read <- c(clipping_block_lengths_curr_read, clipping_block_curr)
          }
        }
        clipping_block_lengths <- c(clipping_block_lengths, sum(clipping_block_lengths_curr_read))
        max_clipping_block_lengths <- c(max_clipping_block_lengths, max(clipping_block_lengths_curr_read))
      } else {
        clipping_block_lengths <- c(clipping_block_lengths, 0)
        max_clipping_block_lengths <- c(max_clipping_block_lengths, 0)
      }
    } else {
      del_block_lengths_curr_read <- 0
      max_del_block_lengths <- c(max_del_block_lengths, 0)
      ins_block_lengths_curr_read <- 0
      max_ins_block_lengths <- c(max_ins_block_lengths, 0)
      clipping_block_lengths_curr_read <- 0
      max_clipping_block_lengths <- c(max_clipping_block_lengths, 0)
    }
  }
  #linear model with double weigth for longest indel
  #sc_thr <- 0 to consider reads carrying a big soft-clipped portion different to reads with big indel (sc_flanking_correction <- 1)
  #otherwise, if sc_thr <- 0.3 soft-clipped portions' lengths are corrected to account for flanking regions length
  sc_thr <- 0
  for (k in 1:length(cigar_strings)) {
    #if allele #1 is the one with the shortest repeat, some reads should carry a big insertion, and soft-clipped portions are probably associated with insertions
    if (mean(max_ins_block_lengths) > mean(max_del_block_lengths)) {
      #soft-clipped regions also include flanking regions -> apply correction
      if (max(max_clipping_block_lengths) > 0) {
        #if there is at least one read carrying an insertion as big as sc_thr*100% the longest soft-clipped portion, evaluate sc_flanking_correction to correct for flanking region length
        if (max(max_ins_block_lengths)/max(max_clipping_block_lengths) > sc_thr) {
          sc_flanking_correction <- max(max_ins_block_lengths)/mean(max_clipping_block_lengths[which(max_clipping_block_lengths != 0)])
        } else {
          sc_flanking_correction <- 1
        }
      } else {
        sc_flanking_correction <- 1
      }
      #if read k carries bigger insertion than deletion it probably comes from allele #2 (score >> 0)
      if (max(max_ins_block_lengths[k], sc_flanking_correction*max_clipping_block_lengths[k]) > max_del_block_lengths[k]) {
        score[k] <- sc_flanking_correction*clipping_block_lengths[k] + ins_block_lengths[k] - del_block_lengths[k] + max(sc_flanking_correction*max_clipping_block_lengths[k], max_ins_block_lengths[k])
        #if read k carries bigger deletion than insertion it probably comes from allele #1 (score ~ 0)
      } else {
        score[k] <- sc_flanking_correction*clipping_block_lengths[k] + ins_block_lengths[k] - del_block_lengths[k] - max(sc_flanking_correction*max_clipping_block_lengths[k], max_del_block_lengths[k])
      }
      #if allele #1 is the one with the longest repeat, some reads should carry a big deletion, and soft-clipped portions are probably associated with deletions
    } else {
      #soft-clipped regions also include flanking regions -> apply correction
      if (max(max_clipping_block_lengths) > 0) {
        #if there is at least one read carrying a deletion as big as sc_thr*100% the longest soft-clipped portion, evaluate sc_flanking_correction to correct for flanking region length
        if (max(max_del_block_lengths)/max(max_clipping_block_lengths) > sc_thr) {
          sc_flanking_correction <- max(max_del_block_lengths)/mean(max_clipping_block_lengths[which(max_clipping_block_lengths != 0)])
        } else {
          sc_flanking_correction <- 1
        }
      } else {
        sc_flanking_correction <- 1
      }
      #if read k carries bigger deletion than insertion, it probably comes from allele #2 (score << 0)
      if (max(max_del_block_lengths[k], max_clipping_block_lengths[k]) > max_ins_block_lengths[k]) {
        score[k] <- ins_block_lengths[k] - sc_flanking_correction*clipping_block_lengths[k] - del_block_lengths[k] - max(sc_flanking_correction*max_clipping_block_lengths[k], max_del_block_lengths[k])
        #if read k carries bigger insertion than deletion it probably comes from allele #1 (score ~ 0)
      } else {
        score[k] <- ins_block_lengths[k] - sc_flanking_correction*clipping_block_lengths[k] - del_block_lengths[k] + max(sc_flanking_correction*max_clipping_block_lengths[k], max_del_block_lengths[k])
      }
    }
  }
  #skip clustering and assign all reads to one allele if studying haploid chromosome or if there are not 2 different maximum non-matching lengths across all reads
  if (haploid_flag == 1 || length(score) < 2) {
    first_allele_reads_fq <- paste0(sample_dir, "/", sample_name, "_reads_first_allele.fastq")
    first_allele_reads_fa <- paste0(sample_dir, "/", sample_name, "_reads_first_allele.fasta")
    first_allele_reads_names <- paste0(sample_dir, "/", sample_name, "_reads_names_first_allele.txt")
    first_allele_reads_fq <- paste0(sample_dir, "/", sample_name, "_reads_first_allele.fastq")
    first_allele_reads_fa <- paste0(sample_dir, "/", sample_name, "_reads_first_allele.fasta")
    #remove outliers which may be associated with somatic mutations
    outliers_score <- boxplot.stats(score, coef = IQR_outliers_coef)$out
    ind_outliers <- which(sort(score) %in% outliers_score)
    score_no_outliers <- score[!score %in% outliers_score]
    num_outliers <- length(score) - length(score_no_outliers)
    if (num_outliers > 0) {
      cat(text = paste0("Sample ", sample_name, ": ", sprintf("%d", num_outliers), " reads possibly associated with somatic mutations have been discarded"), sep = "\n")
      cat(text = paste0("Sample ", sample_name, ": ", sprintf("%d", num_outliers), " reads possibly associated with somatic mutations have been discarded"),  file = logfile, sep = "\n", append = TRUE)
      outliers_readnames <- reads_names[ind_outliers]
      outliers_reads_names <- paste0(sample_dir, "/", sample_name, "_reads_names_outliers.txt")
      write.table(x=outliers_readnames, quote = FALSE, file = outliers_reads_names, row.names = FALSE, col.names = FALSE)
      outliers_reads_fq <- paste0(sample_dir, "/", sample_name, "_reads_outliers.fastq")
      system(paste0(SEQTK, " subseq ", fastq_files[i], " ", outliers_reads_names, " > ", outliers_reads_fq))
      reads_names_no_outliers <- reads_names[-ind_outliers]
    } else {
      reads_names_no_outliers <- reads_names
    }
    num_outliers_reference <- num_outliers
    num_outliers_alternative <- 0
    ind_outliers_reference <- ind_outliers
    ind_outliers_alternative <- c()
    outliers_reference_score <- outliers_score
    outliers_alternative_score <- c()
    write.table(x=reads_names_no_outliers, quote = FALSE, file = first_allele_reads_names, row.names = FALSE, col.names = FALSE)
    system(paste0(SEQTK, " subseq ", fastq_files[i], " ", first_allele_reads_names, " > ", first_allele_reads_fq))
    system(paste0(SEQTK, " seq -A ", first_allele_reads_fq, " > ", first_allele_reads_fa))
    skip_second_allele_flag <- 1
    num_reads_second_allele <- 0
    allelic_ratio_perc_second <- 0
    #cluster reads into 2 groups (ref/alt allele) based on length of the longest portion non matching the reference (first allele)
  } else {
    clusters <- kmeans(score, 2, iter.max = 1000, nstart = 5)
    cluster_reference_id <- which(clusters$centers == min(clusters$centers))[1]
    cluster_alternative_id <- which(clusters$centers == max(clusters$centers))[1]
    cluster_reference_index <- which(clusters$cluster == cluster_reference_id)
    cluster_alternative_index <- which(clusters$cluster == cluster_alternative_id)
    cluster_reference_score <- score[cluster_reference_index]
    cluster_alternative_score <- score[cluster_alternative_index]
    #remove outliers which may be associated with somatic mutations
    outliers_reference_score <- boxplot.stats(cluster_reference_score, coef = IQR_outliers_coef)$out
    ind_outliers_reference <- which(sort(cluster_reference_score) %in% outliers_reference_score)
    score_reference_no_outliers <- cluster_reference_score[!cluster_reference_score %in% outliers_reference_score]
    num_outliers_reference <- length(cluster_reference_score) - length(score_reference_no_outliers)
    if (num_outliers_reference > 0) {
      cluster_reference_index_no_outliers <- cluster_reference_index[-ind_outliers_reference]
      cluster_reference_readnames <- (reads_names[cluster_reference_index])[-ind_outliers_reference]
    } else {
      cluster_reference_index_no_outliers <- cluster_reference_index
      cluster_reference_readnames <- reads_names[cluster_reference_index]
    }
    outliers_alternative_score <- boxplot.stats(cluster_alternative_score, coef = 2)$out
    ind_outliers_alternative <- which(sort(cluster_alternative_score) %in% outliers_alternative_score)
    score_alternative_no_outliers <- cluster_alternative_score[!cluster_alternative_score %in% outliers_alternative_score]
    num_outliers_alternative <- length(cluster_alternative_score) - length(score_alternative_no_outliers)
    if (num_outliers_alternative > 0) {
      cluster_alternative_index_no_outliers <- cluster_alternative_index[-ind_outliers_alternative]
      cluster_alternative_readnames <- (reads_names[cluster_alternative_index])[-ind_outliers_alternative]
    } else {
      cluster_alternative_index_no_outliers <- cluster_alternative_index
      cluster_alternative_readnames <- reads_names[cluster_alternative_index]
    }
    median_cluster_reference_nonreflen <- clusters$centers[cluster_reference_id]
    median_cluster_alternative_nonreflen <- clusters$centers[cluster_alternative_id]
    clusters$median <- c(median_cluster_reference_nonreflen, median_cluster_alternative_nonreflen)
    first_allele_reads_names <- paste0(sample_dir, "/", sample_name, "_reads_names_first_allele.txt")
    first_allele_reads_fq <- paste0(sample_dir, "/", sample_name, "_reads_first_allele.fastq")
    first_allele_reads_fa <- paste0(sample_dir, "/", sample_name, "_reads_first_allele.fasta")
    write.table(x=cluster_reference_readnames, quote = FALSE, file = first_allele_reads_names, row.names = FALSE, col.names = FALSE)
    system(paste0(SEQTK, " subseq ", fastq_files[i], " ", first_allele_reads_names, " > ", first_allele_reads_fq))
    system(paste0(SEQTK, " seq -A ", first_allele_reads_fq, " > ", first_allele_reads_fa))
    second_allele_reads_names <- paste0(sample_dir, "/", sample_name, "_reads_names_second_allele.txt")
    second_allele_reads_fq <- paste0(sample_dir, "/", sample_name, "_reads_second_allele.fastq")
    second_allele_reads_fa <- paste0(sample_dir, "/", sample_name, "_reads_second_allele.fasta")
    write.table(x=cluster_alternative_readnames, quote = FALSE, file = second_allele_reads_names, row.names = FALSE, col.names = FALSE)
    system(paste0(SEQTK, " subseq ", fastq_files[i], " ", second_allele_reads_names, " > ", second_allele_reads_fq))
    system(paste0(SEQTK, " seq -A ", second_allele_reads_fq, " > ", second_allele_reads_fa))
    outliers_readnames <- reads_names[c(ind_outliers_reference, ind_outliers_alternative)]
    outliers_reads_names <- paste0(sample_dir, "/", sample_name, "_reads_names_outliers.txt")
    write.table(x=outliers_readnames, quote = FALSE, file = outliers_reads_names, row.names = FALSE, col.names = FALSE)
    outliers_reads_fq <- paste0(sample_dir, "/", sample_name, "_reads_outliers.fastq")
    system(paste0(SEQTK, " subseq ", fastq_files[i], " ", outliers_reads_names, " > ", outliers_reads_fq))
  }
  png(paste0(sample_dir, "/", sample_name, "_per_read_score.png"))
  plot(1:length(score), sort(score), xlab = "Reads", ylab = "Score (bp)", main = "Per-read score")
  if (num_outliers_reference > 0) {
    lines(ind_outliers_reference, sort(outliers_reference_score), col = "red", type = "p")
  }
  if (num_outliers_alternative > 0) {
    lines((length(cluster_reference_score) + ind_outliers_alternative), sort(outliers_alternative_score), col = "red", type = "p")
  }
  dev.off()
  num_outliers <- num_outliers_reference + num_outliers_alternative
  if (num_outliers > 0) {
    cat(text = paste0("Sample ", sample_name, ": ", sprintf("%d", num_outliers), " reads possibly associated with somatic mutations have been discarded"), sep = "\n")
    cat(text = paste0("Sample ", sample_name, ": ", sprintf("%d", num_outliers), " reads possibly associated with somatic mutations have been discarded"),  file = logfile, sep = "\n", append = TRUE)
    ind_outliers <- which(score %in% c(outliers_reference_score, outliers_alternative_score))
    outliers_readnames <- reads_names[ind_outliers]
    outliers_reads_names <- paste0(sample_dir, "/", sample_name, "_reads_names_outliers.txt")
    write.table(x=outliers_readnames, quote = FALSE, file = outliers_reads_names, row.names = FALSE, col.names = FALSE)
    outliers_reads_fq <- paste0(sample_dir, "/", sample_name, "_reads_outliers.fastq")
    system(paste0(SEQTK, " subseq ", fastq_files[i], " ", outliers_reads_names, " > ", outliers_reads_fq))
  }
  if (num_outliers > 0) {
    reads_names_no_outliers <- reads_names[-ind_outliers]
  } else {
    reads_names_no_outliers <- reads_names
  }
  if (skip_first_allele_flag != 1) {
    num_reads_first_allele <- as.double(system(paste0("cat ", first_allele_reads_fa, " | grep \"^>\" | wc -l"), intern=TRUE))
  }
  if (skip_second_allele_flag != 1) {
    num_reads_second_allele <- as.double(system(paste0("cat ", second_allele_reads_fa, " | grep \"^>\" | wc -l"), intern=TRUE))
  }
  #create consensus sequence for Allele #1
  allelic_ratio_first <- num_reads_first_allele/num_reads_sample
  allelic_ratio_perc_first <- allelic_ratio_first*100
  allelic_ratio_second <- num_reads_second_allele/num_reads_sample
  allelic_ratio_perc_second <- allelic_ratio_second*100
  unassigned_reads_perc <- 100 - allelic_ratio_perc_first - allelic_ratio_perc_second
  cat(text = paste0("Sample ", sample_name, ": ", sprintf("%.2f", allelic_ratio_perc_first), "% reads assigned to Allele #1; ", sprintf("%.2f", allelic_ratio_perc_second), "% reads assigned to Allele #2; ", sprintf("%.2f", unassigned_reads_perc), "% reads unassigned"), sep = "\n")
  cat(text = paste0("Sample ", sample_name, ": ", sprintf("%.2f", allelic_ratio_perc_first), "% reads assigned to Allele #1; ", sprintf("%.2f", allelic_ratio_perc_second), "% reads assigned to Allele #2; ", sprintf("%.2f", unassigned_reads_perc), "% reads unassigned"),  file = logfile, sep = "\n", append = TRUE)
  #if not at least 3 reads are assigned to one of the two alleles, the sample is skipped
  if (max(num_reads_first_allele, num_reads_second_allele) < 3) {
    max_reads_allele <- max(num_reads_first_allele, num_reads_second_allele)
    num_allele <- which(c(num_reads_first_allele, num_reads_second_allele) == max_reads_allele)[1]                        
    cat(text = paste0("WARNING: Only ", max_reads_allele, " reads available for sample ", sample_name, " for Allele #", num_allele, "; skipping sample"), sep = "\n")
    cat(text = paste0("WARNING: Only ", max_reads_allele, " reads available for sample ", sample_name, " for Allele #", num_allele, "; skipping sample"),  file = logfile, sep = "\n", append = TRUE)
    next
  }
  if (num_reads_first_allele < 3 || allelic_ratio_first < min_maf) {
    skip_first_allele_flag <- 1
  } else if (num_reads_second_allele < 3 || allelic_ratio_second < min_maf) {
    skip_second_allele_flag <- 1
  } 
  
  target_reads_contig <- TRC
  target_reads_polishing <- TRP
  
  if (skip_first_allele_flag == 1) {
    cat(text = paste0("WARNING: Only ", num_reads_first_allele, " reads (", sprintf("%.2f", allelic_ratio_perc_first), "%) from sample ", sample_name, " have been assigned to Allele #1; skipping"), sep = "\n")
    cat(text = paste0("WARNING: Only ", num_reads_first_allele, " reads (", sprintf("%.2f", allelic_ratio_perc_first), "%) from sample ", sample_name, " have been assigned to Allele #1; skipping"),  file = logfile, sep = "\n", append = TRUE)
  } else {
    if (num_reads_first_allele < target_reads_contig) {
      target_reads_contig <- num_reads_first_allele
      target_reads_polishing <- num_reads_first_allele
      cat(text = paste0("WARNING: Only ", num_reads_first_allele, " reads available for sample ", sample_name, " for Allele #1"), sep = "\n")
      cat(text = paste0("WARNING: Only ", num_reads_first_allele, " reads available for sample ", sample_name, " for Allele #1"),  file = logfile, sep = "\n", append = TRUE)
    } 
    plurality_value <- 0.15*target_reads_contig
    draft_contig_first_tmp1 <- paste0(sample_dir, "/", sample_name, "_non_polished.contig_first_tmp1.fasta")
    draft_contig_first_tmp2 <- paste0(sample_dir, "/", sample_name, "_non_polished.contig_first_tmp2.fasta")
    draft_contig_first <- paste0(sample_dir, "/", sample_name, "_non_polished.contig_first.fasta")
    first_allele_untrimmed <- paste0(home_dir, "/", sample_name, "/", sample_name, "_first_allele_untrimmed.fasta")
    first_allele <- paste0(home_dir, "/", sample_name, "_first_allele.fasta")
    sequences <- readDNAStringSet(fasta_files[i], "fasta")
    ws <- width(sequences)
    amplicon_length <- ceiling(mean(ws))
    draft_reads_fq_first <- paste0(sample_dir, "/", sample_name, "_draft_", target_reads_contig, "_reads_first.fastq")
    draft_reads_fa_first <- paste0(sample_dir, "/", sample_name, "_draft_", target_reads_contig, "_reads_first.fasta")
    seed <- 1
    system(paste0(SEQTK, " sample -s ", seed , " ", first_allele_reads_fq, " ",  target_reads_contig, " > ", draft_reads_fq_first))
    system(paste0(SEQTK, " seq -A ", draft_reads_fq_first, " > ", draft_reads_fa_first))
    mfa_file_first <- gsub(pattern = "\\.fasta$", replacement = ".mfa", x = draft_reads_fa_first)
    #system(paste0(MAFFT, " --auto --thread ", num_threads, " --threadit 0 --adjustdirectionaccurately ", draft_reads_fa_first, " > ", mfa_file_first)) #threadit 0 if you need reproducible results
    if (fast_alignment_flag == 1) {
      system(paste0(MAFFT, " --auto --thread ", num_threads, " --adjustdirectionaccurately ", draft_reads_fa_first, " > ", mfa_file_first))
    } else {
      system(paste0(MAFFT , "-linsi --thread ", num_threads, " --adjustdirectionaccurately ", draft_reads_fa_first, " > ", mfa_file_first))
    }
    system(paste0(CONS, " -sequence ", mfa_file_first, " -plurality ", plurality_value, " -outseq ", draft_contig_first_tmp1))
    system(paste0("sed 's/[nN]//g' ", draft_contig_first_tmp1, " > ", draft_contig_first_tmp2))
    DNAStringSet_obj <- readDNAStringSet(draft_contig_first_tmp2, "fasta")
    DNAStringSet_obj_renamed <- DNAStringSet_obj
    original_headers <- names(DNAStringSet_obj)
    sequences <- seq(DNAStringSet_obj)
    names(DNAStringSet_obj_renamed) <- "Allele_number_1"
    writeXStringSet(x = DNAStringSet_obj_renamed, filepath = draft_contig_first, format = "fasta", width = 20000)
    if (exists("fast5_dir") && pair_strands_flag != 1) {
      qual_filter <- 0
      reads_polishing_fq_first <- paste0(sample_dir, "/", sample_name, "_polishing_", target_reads_polishing, "_reads_first.fastq")
      reads_polishing_fa_first <- paste0(sample_dir, "/", sample_name, "_polishing_", target_reads_polishing, "_reads_first.fasta")  
      seed_polishing <- 2 
      system(paste0(SEQTK, " sample -s ", seed_polishing, " ", first_allele_reads_fq, " ", target_reads_polishing, " > ", reads_polishing_fq_first))
      system(paste0(SEQTK, " seq -A ", reads_polishing_fq_first, " > ", reads_polishing_fa_first))
      bam_file_first <- paste0(sample_dir, "/", sample_name, "_first.bam")
      system(paste0(MINIMAP2, " -ax map-ont ", draft_contig_first, " ", reads_polishing_fa_first, " | ", SAMTOOLS, " view -h -q 55 -F 2048 | " , SAMTOOLS, " sort -o ", bam_file_first, " -T reads.tmp"))
      system(paste0(SAMTOOLS," index ", bam_file_first))
      cat(text = paste0("Running nanopolish for sample ", sample_name, " - Allele #1"), sep = "\n")
      cat(text = "Indexing reads", sep = "\n")
      if (seq_sum_flag == 1) {
        system(paste0(NANOPOLISH, " index -d ", fast5_dir, " -s ", sequencing_summary, " ", reads_polishing_fa_first))
      } else {
        system(paste0(NANOPOLISH, " index -d ", fast5_dir, " ", reads_polishing_fa_first))
      }
      cat(text = paste0("Running nanopolish for consensus polishing of sample ", sample_name, " - Allele #1"), sep = "\n")
      output_vcf_first <- paste0(sample_dir, "/", "nanopolish_output_first.vcf")
      output_vcf_filtered_first <- paste0(sample_dir, "/", "nanopolish_output_filtered_first.vcf")
      system(paste0(NANOPOLISH, " variants --consensus --reads ", reads_polishing_fa_first, " --bam ", bam_file_first, " --genome ", draft_contig_first, " -p 1 --threads ", num_threads, " -o ", output_vcf_first))
      system(paste0("cat ", output_vcf_first, " | grep \"^#\" > ", output_vcf_filtered_first))
      system(paste0("cat ", output_vcf_first, " | grep -v \"^#\" | awk '$6 > ", qual_filter, " {print}' >> ", output_vcf_filtered_first))
      system(paste0(NANOPOLISH, " vcf2fasta -g ", draft_contig_first, " ", output_vcf_filtered_first, " > ", first_allele_untrimmed))
      system(paste0(SEQTK, " trimfq ", first_allele_untrimmed, " -b ", primers_length, " -e ", primers_length, " > ", first_allele))
    } else {
      system(paste0(SEQTK, " trimfq ", draft_contig_first, " -b ", primers_length, " -e ", primers_length, " > ", first_allele))
    }
    setwd(home_dir)
    system(paste0(TRF, " ", first_allele, " 2 7 7 80 10 50 500"))
    system(paste0(TRF, " ", first_allele, " 2 3 5 80 10 14 500"))
    system(paste0(TRF, " ", first_allele, " 2 500 500 80 10 50 500"))
    default_trf_par_files <-  list.files(path = home_dir, pattern = "\\.2\\.7\\.7\\.80\\.10\\.50\\.500", full.names = TRUE)
    stringent_trf_par_files <- list.files(path = home_dir, pattern = "\\.2\\.500\\.500\\.80\\.10\\.50\\.500", full.names = TRUE)
    lenient_trf_par_files <- list.files(path = home_dir, pattern = "\\.2\\.3\\.5\\.80\\.10\\.14\\.500", full.names = TRUE)
    system(paste0("mv ", paste0(lenient_trf_par_files, collapse = " "), " ", sample_dir))
    system(paste0("mv ", paste0(stringent_trf_par_files, collapse = " "), " ", sample_dir))
  }
  #create consensus sequence for Allele #2
  target_reads_contig <- TRC
  target_reads_polishing <- TRP
  
  if (skip_second_allele_flag == 1) {
    cat(text = paste0("WARNING: Only ", num_reads_second_allele, " reads (", sprintf("%.2f", allelic_ratio_perc_second), "%) from sample ", sample_name, " have been assigned to Allele #2; skipping"), sep = "\n")
    cat(text = paste0("WARNING: Only ", num_reads_second_allele, " reads (", sprintf("%.2f", allelic_ratio_perc_second), "%) from sample ", sample_name, " have been assigned to Allele #2; skipping"),  file = logfile, sep = "\n", append = TRUE)
  } else {
    if (num_reads_second_allele < target_reads_contig) {
      target_reads_contig <- num_reads_second_allele
      target_reads_polishing <- num_reads_second_allele
      cat(text = paste0("WARNING: Only ", num_reads_second_allele, " reads available for sample ", sample_name, " for Allele #2"), sep = "\n")
      cat(text = paste0("WARNING: Only ", num_reads_second_allele, " reads available for sample ", sample_name, " for Allele #2"),  file = logfile, sep = "\n", append = TRUE)
    } 
    plurality_value <- 0.15*target_reads_contig
    draft_contig_second_tmp1 <- paste0(sample_dir, "/", sample_name, "_non_polished.contig_second_tmp1.fasta")
    draft_contig_second_tmp2 <- paste0(sample_dir, "/", sample_name, "_non_polished.contig_second_tmp2.fasta")
    draft_contig_second <- paste0(sample_dir, "/", sample_name, "_non_polished.contig_second.fasta")
    second_allele_untrimmed <- paste0(home_dir, "/", sample_name, "/", sample_name, "_second_allele_untrimmed.fasta")
    second_allele <- paste0(home_dir, "/", sample_name, "_second_allele.fasta")
    sequences <- readDNAStringSet(fasta_files[i], "fasta")
    ws <- width(sequences)
    amplicon_length <- ceiling(mean(ws))
    draft_reads_fq_second <- paste0(sample_dir, "/", sample_name, "_draft_", target_reads_contig, "_reads_second.fastq")
    draft_reads_fa_second <- paste0(sample_dir, "/", sample_name, "_draft_", target_reads_contig, "_reads_second.fasta")
    seed <- 1
    system(paste0(SEQTK, " sample -s ", seed , " ", second_allele_reads_fq, " ",  target_reads_contig, " > ", draft_reads_fq_second))
    system(paste0(SEQTK, " seq -A ", draft_reads_fq_second, " > ", draft_reads_fa_second))
    mfa_file_second <- gsub(pattern = "\\.fasta$", replacement = ".mfa", x = draft_reads_fa_second)
    #system(paste0(MAFFT, " --auto --thread ", num_threads, " --threadit 0 --adjustdirectionaccurately ", draft_reads_fa_second, " > ", mfa_file_second)) #threadit 0 if you need reproducible results
    if (fast_alignment_flag == 1) {
      system(paste0(MAFFT, " --auto --thread ", num_threads, " --adjustdirectionaccurately ", draft_reads_fa_second, " > ", mfa_file_second))
    } else {
      system(paste0(MAFFT, "-linsi --thread ", num_threads, " --adjustdirectionaccurately ", draft_reads_fa_second, " > ", mfa_file_second))
    }
    system(paste0(CONS, " -sequence ", mfa_file_second, " -plurality ", plurality_value, " -outseq ", draft_contig_second_tmp1))
    system(paste0("sed 's/[nN]//g' ", draft_contig_second_tmp1, " > ", draft_contig_second_tmp2))
    DNAStringSet_obj <- readDNAStringSet(draft_contig_second_tmp2, "fasta")
    DNAStringSet_obj_renamed <- DNAStringSet_obj
    original_headers <- names(DNAStringSet_obj)
    sequences <- seq(DNAStringSet_obj)
    names(DNAStringSet_obj_renamed) <- "Allele_number_2"
    writeXStringSet(x = DNAStringSet_obj_renamed, filepath = draft_contig_second, format = "fasta", width = 20000)
    if (exists("fast5_dir") && pair_strands_flag != 1) {
      qual_filter <- 0
      reads_polishing_fq_second <- paste0(sample_dir, "/", sample_name, "_polishing_", target_reads_polishing, "_reads_second.fastq")
      reads_polishing_fa_second <- paste0(sample_dir, "/", sample_name, "_polishing_", target_reads_polishing, "_reads_second.fasta")  
      seed_polishing <- 2 
      system(paste0(SEQTK, " sample -s ", seed_polishing, " ", draft_reads_fq_second, " ", target_reads_polishing, " > ", reads_polishing_fq_second))
      system(paste0(SEQTK, " seq -A ", reads_polishing_fq_second, " > ", reads_polishing_fa_second))
      bam_file_second <- paste0(sample_dir, "/", sample_name, "_second.bam")
      system(paste0(MINIMAP2, " -ax map-ont ", draft_contig_second, " ", reads_polishing_fa_second, " | ", SAMTOOLS, " view -h -q 55 -F 2048 | " , SAMTOOLS, " sort -o ", bam_file_second, " -T reads.tmp"))
      system(paste0(SAMTOOLS," index ", bam_file_second))
      cat(text = paste0("Running nanopolish for sample ", sample_name, " - Allele #2"), sep = "\n")
      cat(text = "Indexing reads", sep = "\n")
      if (seq_sum_flag == 1) {
        system(paste0(NANOPOLISH, " index -d ", fast5_dir, " -s ", sequencing_summary, " ", reads_polishing_fa_second))
      } else {
        system(paste0(NANOPOLISH, " index -d ", fast5_dir, " ", reads_polishing_fa_second))
      }
      cat(text = paste0("Running nanopolish for consensus polishing of sample ", sample_name, " - Allele #2"), sep = "\n")
      output_vcf_second <- paste0(sample_dir, "/", "nanopolish_output_second.vcf")
      output_vcf_filtered_second <- paste0(sample_dir, "/", "nanopolish_output_filtered_second.vcf")
      system(paste0(NANOPOLISH, " variants --consensus --reads ", reads_polishing_fa_second, " --bam ", bam_file_second, " --genome ", draft_contig_second, " -p 1 --threads ", num_threads, " -o ", output_vcf_second))
      system(paste0("cat ", output_vcf_second, " | grep \"^#\" > ", output_vcf_filtered_second))
      system(paste0("cat ", output_vcf_second, " | grep -v \"^#\" | awk '$6 > ", qual_filter, " {print}' >> ", output_vcf_filtered_second))
      system(paste0(NANOPOLISH, " vcf2fasta -g ", draft_contig_second, " ", output_vcf_filtered_second, " > ", second_allele_untrimmed))
      system(paste0(SEQTK, " trimfq ", second_allele_untrimmed, " -b ", primers_length, " -e ", primers_length, " > ", second_allele))
    } else {
      system(paste0(SEQTK, " trimfq ", draft_contig_second, " -b ", primers_length, " -e ", primers_length, " > ", second_allele))
    }
    setwd(home_dir)
    system(paste0(TRF, " ", second_allele, " 2 7 7 80 10 50 500"))
    system(paste0(TRF, " ", second_allele, " 2 3 5 80 10 14 500"))
    system(paste0(TRF, " ", second_allele, " 2 500 500 80 10 50 500"))
    default_trf_par_files <-  list.files(path = home_dir, pattern = "\\.2\\.7\\.7\\.80\\.10\\.50\\.500", full.names = TRUE)
    stringent_trf_par_files <- list.files(path = home_dir, pattern = "\\.2\\.500\\.500\\.80\\.10\\.50\\.500", full.names = TRUE)
    lenient_trf_par_files <- list.files(path = home_dir, pattern = "\\.2\\.3\\.5\\.80\\.10\\.14\\.500", full.names = TRUE)
    system(paste0("mv ", paste0(lenient_trf_par_files, collapse = " "), " ", sample_dir))
    system(paste0("mv ", paste0(stringent_trf_par_files, collapse = " "), " ", sample_dir))
  }
}
