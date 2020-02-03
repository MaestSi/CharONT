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
  cat(paste0("Usage: Rscript CharONT.R <analysis_dir>"), sep = "\n")
  cat(paste0("Note that config_CharONT.R must be in the same directory of CharONT.R"), sep = "\n")
  cat(paste0("<analysis_dir>: directory containing fastq and fasta files for each sample"), sep = "\n")
  stop(simpleError(sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "))))
}

if (length(args) == 1) {
  analysis_dir <- args[1]
  if (!dir.exists(analysis_dir)) {
    stop(paste0(analysis_dir, " directory does not exist!"))
  }
} else {
  stop("Home directory has to be provided")
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

logfile <- paste0(analysis_dir, "/logfile.txt")

fastq_files <- list.files(path = analysis_dir, pattern = "BC\\d+\\.fastq", full.names = TRUE)
fasta_files <- list.files(path = analysis_dir, pattern = "BC\\d+\\.fasta", full.names = TRUE)

if (length(fastq_files) == 0) {
  stop(paste0("No fastq files BC<numbers>.fastq in ", analysis_dir, " directory"))
}

if (length(fasta_files) > 0) {
  cat(text = paste0("Processing fasta files ", paste0(basename(fasta_files), collapse = ", ")), sep = "\n")
} else {
  for (i in 1:length(fastq_files)) {
    fasta_file_curr <- paste0(analysis_dir, "/", gsub(pattern = "\\.fastq$", replacement = "\\.fasta", x = basename(fastq_files[i])))
    fasta_files <- c(fasta_files, fasta_file_curr)
    system(paste0(SEQTK, " seq -A ", fastq_files[i], " > ", fasta_file_curr))
  }
  cat(text = paste0("Processing fasta files ", paste0(basename(fasta_files), collapse = ", ")), sep = "\n")
}

available_medaka_models <- system(paste0(MEDAKA, " tools list_models"), intern = TRUE)[1]

if (length(grep(pattern = medaka_model, x = available_medaka_models)) > 0) {
  cat(text = paste0("Medaka model: ", medaka_model), sep = "\n")
  cat(text = paste0("Medaka model: ", medaka_model), file = logfile, sep = "\n", append = TRUE)
  cat(text = "\n", file = logfile, append = TRUE)
} else {
  cat(text = paste0("Medaka model: ", medaka_model, " is not available; default model r941_min_high_g344 was selected"), sep = "\n")
  cat(text = paste0("Medaka model: ", medaka_model, " is not available; default model r941_min_high_g344 was selected"), file = logfile, sep = "\n", append = TRUE)
  cat(text = "\n", file = logfile, append = TRUE)  
  medaka_model <- "r941_min_high_g344"
}

target_reads_consensus <- TRC
target_reads_polishing <- TRP
THR <- 0.85
plurality_value <- 0.15*target_reads_consensus
max_num_reads_clustering <- 500
seed <- 1

#cycle over fasta files
for (i in 1:length(fasta_files)) {
  skip_first_allele_flag <- 0
  skip_second_allele_flag <- 0
  sample_dir <- gsub(pattern = "\\.fasta", replacement = "", x = fasta_files[i])
  sample_name <- basename(sample_dir)
  subset_reads_fa <- paste0(sample_dir, "/", sample_name, "_subset_reads.fasta")
  dir.create(sample_dir)
  vsearch_clustering_dir <- paste0(sample_dir, "/clustering_vsearch")
  dir.create(vsearch_clustering_dir)
  #perform vsearch clustering and create a preliminary version for Allele #1
  system(paste0(SEQTK, " sample -s ", seed , " ", fasta_files[i], " ",  max_num_reads_clustering, " > ", subset_reads_fa))
  ids_mac_first_preliminary <- paste0(sample_dir, "/", sample_name, "_reads_ids_mac.txt")
  mac_fa_first_preliminary <- paste0(sample_dir, "/", sample_name, "_reads_mac.fasta")
  mac_fq_first_preliminary <- paste0(sample_dir, "/", sample_name, "_reads_mac.fastq")
  system(paste0(VSEARCH, " --cluster_smallmem ", subset_reads_fa, " --usersort --id ", THR, " --iddef 2 --clusterout_sort --fasta_width 0 --strand both --sizeout --consout ", vsearch_clustering_dir, "/", sample_name, "_consensus.fasta --clusters ", vsearch_clustering_dir, "/", sample_name, "_cluster"))
  centroid_mac <- system(paste0("head -n1 ", vsearch_clustering_dir, "/", sample_name, "_consensus.fasta"), intern = TRUE)
  id_centroid <-  system(paste0("echo ", "\"", centroid_mac, "\"", " | sed 's/centroid=//g' | sed 's/;seqs.*$//g'"), intern = TRUE)
  clusters_vsearch <- list.files(path = vsearch_clustering_dir, pattern = paste0(sample_name, "_cluster"), full.names = TRUE)
  for (j in 1:length(clusters_vsearch)) {
    cluster_match <- system(paste0("cat ", clusters_vsearch[j], " | grep \"", id_centroid, "\""), intern = TRUE)
    if (length(cluster_match) > 0) {
      mac_file <- clusters_vsearch[j]
      system(paste0("cat ", mac_file, " | grep ", "\"",  "^>", "\"", "  | sed 's/^>//' | sed 's/;.*$//' > ", ids_mac_first_preliminary))
      system(paste0(SEQTK, " subseq ", fasta_files[i], " ", ids_mac_first_preliminary, " > ", mac_fa_first_preliminary))
      system(paste0(SEQTK, " subseq ", fastq_files[i], " ", ids_mac_first_preliminary, " > ", mac_fq_first_preliminary))
      break
    }
  }
  num_reads_sample <- as.double(system(paste0("cat ", fasta_files[i], " | grep \"^>\" | wc -l"), intern=TRUE))
  num_reads_mac_first_preliminary <- as.double(system(paste0("cat ", mac_fa_first_preliminary, " | grep \"^>\" | wc -l"), intern=TRUE))
  target_reads_consensus <- TRC
  first_allele_preliminary_tmp1 <- paste0(sample_dir, "/", sample_name, "_preliminary_first_allele_tmp1.fasta")
  first_allele_preliminary_tmp2 <- paste0(sample_dir, "/", sample_name, "_preliminary_first_allele_tmp2.fasta")
  first_allele_preliminary <- paste0(analysis_dir, "/", sample_name, "/", sample_name, "_preliminary_first_allele.fasta")
  if (num_reads_mac_first_preliminary < 3) {
    system(paste0("head -n2 ", vsearch_clustering_dir, "/", sample_name, "_consensus.fasta > ", first_allele_preliminary))
  } else {
    if (num_reads_mac_first_preliminary < target_reads_consensus) {
      target_reads_consensus <- num_reads_mac_first_preliminary
    }
    plurality_value <- 0.15*target_reads_consensus
    sequences <- readDNAStringSet(fasta_files[i], "fasta")
    ws <- width(sequences)
    amplicon_length <- ceiling(mean(ws))
    draft_reads_fq_first_preliminary <- paste0(sample_dir, "/", sample_name, "_draft_", target_reads_consensus, "_reads_preliminary_first_allele.fastq")
    draft_reads_fa_first_preliminary <- paste0(sample_dir, "/", sample_name, "_draft_", target_reads_consensus, "_reads_preliminary_first_allele.fasta")
    system(paste0(SEQTK, " sample -s ", seed , " ", mac_fq_first_preliminary, " ",  target_reads_consensus, " > ", draft_reads_fq_first_preliminary))
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
  min_clipped_len <- 500
  
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
    #if there are non matching bases
    if (nonref_starting_coords[1] != -1) {
      #if there are deletions
      if (del_starting_coords[1] != -1) {
        #extract lengths of deletions
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
      #if there are insertions
      if (ins_starting_coords[1] != -1) {
        #extract lengths of insertions
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
      #if there are soft-clipped bases
      if (clipping_starting_coords[1] != -1) {
        #extract length of soft-clipped bases if they are longer than min_clipped_len
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
  #longest indel model
  #cycle over reads
  for (k in 1:length(cigar_strings)) {
    #if allele #1 is the one with the shortest repeat, some reads should carry a big insertion, and soft-clipped portions are probably associated with insertions
    if (mean(max_ins_block_lengths) > mean(max_del_block_lengths)) {
      #if read k carries bigger insertion than deletion it probably comes from allele #2 (score >> 0)
      if (max(max_ins_block_lengths[k], max_clipping_block_lengths[k]) > max_del_block_lengths[k]) {
        score[k] <- max(max_clipping_block_lengths[k], max_ins_block_lengths[k])
      #if read k carries bigger deletion than insertion it probably comes from allele #1 (score ~ 0)
      } else {
        score[k] <- - max_del_block_lengths[k]
      }
    #if allele #1 is the one with the longest repeat, some reads should carry a big deletion, and soft-clipped portions are probably associated with deletions
    } else {
        #if read k carries bigger deletion than insertion, it probably comes from allele #2 (score << 0)
        if (max(max_del_block_lengths[k], max_clipping_block_lengths[k]) > max_ins_block_lengths[k]) {
          score[k] <- - max(max_clipping_block_lengths[k], max_del_block_lengths[k])
          #if read k carries bigger insertion than deletion it probably comes from allele #1 (score ~ 0)
        } else {
          score[k] <- max_ins_block_lengths[k]
        }
    }
    if (length(which(is.na(score) == TRUE)) > 0) {
      ind_NA_score <- which(is.na(score) == TRUE)
      score[ind_NA_score] <- 0
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
    ind_outliers <- which(score %in% outliers_score)
    score_no_outliers <- score[!score %in% outliers_score]
    num_outliers <- length(score) - length(score_no_outliers)
    if (num_outliers > 0) {
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
    score_thr <- min(score)
    #cluster reads into 2 groups (ref/alt allele) based on length of the longest portion non matching the reference (first allele)
  } else {
    #remove outliers that might negatively affect clustering
    preclustering_outliers_score <- boxplot.stats(score, coef = IQR_outliers_coef)$out
    preclustering_score_no_outliers <- score[!score %in% preclustering_outliers_score]
    #do clustering
    clusters <- kmeans(preclustering_score_no_outliers, 2, iter.max = 1000, nstart = 5)
    #check if preclustering outliers should be assigned to reference or alternative cluster
    outliers_centers_sorted <- sort(c(preclustering_outliers_score, clusters$centers))
    ind_center_ref <- which(abs(outliers_centers_sorted) == min(abs(clusters$centers)))
    ind_center_alt <- which(abs(outliers_centers_sorted) == max(abs(clusters$centers)))
    #determine which center is the center of reference cluster, and assign outliers to one of the clusters
    if (ind_center_ref < ind_center_alt) {
      preclustering_outliers_reference_score <- intersect(preclustering_outliers_score, outliers_centers_sorted[1:(ind_center_ref - 1)])
      preclustering_outliers_alternative_score <- intersect(preclustering_outliers_score, outliers_centers_sorted[(ind_center_ref + 2):length(outliers_centers_sorted)])
    } else {
      preclustering_outliers_alternative_score <- intersect(preclustering_outliers_score, outliers_centers_sorted[1:(ind_center_alt - 1)])
      preclustering_outliers_reference_score <- intersect(preclustering_outliers_score, outliers_centers_sorted[(ind_center_alt + 2):length(outliers_centers_sorted)])
    }
    #find index of preclustering outliers
    ind_outliers_reference_preclustering <- which(score %in% preclustering_outliers_reference_score)
    ind_outliers_alternative_preclustering <- which(score %in% preclustering_outliers_alternative_score)
    #find index and score of reference reads
    cluster_reference_id <- which(abs(clusters$centers) == min(abs(clusters$centers)))[1]
    cluster_reference_index_tmp <- which(clusters$cluster == cluster_reference_id)
    cluster_reference_score <- c(preclustering_score_no_outliers[cluster_reference_index_tmp], preclustering_outliers_reference_score)
    cluster_reference_index <- which(score %in% cluster_reference_score)
    #find index and score of alternative reads
    cluster_alternative_id <- which(abs(clusters$centers) == max(abs(clusters$centers)))[1]
    cluster_alternative_index_tmp <- which(clusters$cluster == cluster_alternative_id)
    cluster_alternative_score <- c(preclustering_score_no_outliers[cluster_alternative_index_tmp], preclustering_outliers_alternative_score)
    cluster_alternative_index <- which(score %in% cluster_alternative_score)
    #calculate score threshold to split reads into two clusters 
    if (ind_center_ref < ind_center_alt) {
      score_thr <- (max(cluster_reference_score) + min(cluster_alternative_score))/2
    } else {
      score_thr <- (max(cluster_alternative_score) +  min(cluster_reference_score))/2
    }
    #remove outliers from reference cluster which may be associated with somatic mutations
    outliers_reference_score <- boxplot.stats(cluster_reference_score, coef = IQR_outliers_coef)$out
    ind_outliers_reference <- which(cluster_reference_score %in% outliers_reference_score)
    score_reference_no_outliers <- cluster_reference_score[!cluster_reference_score %in% outliers_reference_score]
    num_outliers_reference <- length(cluster_reference_score) - length(score_reference_no_outliers)
    if (num_outliers_reference > 0) {
      cluster_reference_index_no_outliers <- cluster_reference_index[-ind_outliers_reference]
      cluster_reference_readnames <- reads_names[cluster_reference_index_no_outliers]
    } else {
      cluster_reference_index_no_outliers <- cluster_reference_index
      cluster_reference_readnames <- reads_names[cluster_reference_index]
    }
    #remove outliers from alternative cluster which may be associated with somatic mutations
    outliers_alternative_score <- boxplot.stats(cluster_alternative_score, coef = IQR_outliers_coef)$out
    ind_outliers_alternative <- which(sort(cluster_alternative_score) %in% outliers_alternative_score)
    score_alternative_no_outliers <- cluster_alternative_score[!cluster_alternative_score %in% outliers_alternative_score]
    num_outliers_alternative <- length(cluster_alternative_score) - length(score_alternative_no_outliers)
    if (num_outliers_alternative > 0) {
      cluster_alternative_index_no_outliers <- cluster_alternative_index[-ind_outliers_alternative]
      cluster_alternative_readnames <- reads_names[cluster_alternative_index_no_outliers]
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
  num_outliers <- num_outliers_reference + num_outliers_alternative
  png(paste0(sample_dir, "/", sample_name, "_reads_scores.png"))
  plot(1:length(score), sort(score), xlab = "Reads", ylab = "Score (bp)", main = "Reads scores")
  abline(h = score_thr, col = "blue", lty = 2, lwd = 3)
  if(num_outliers > 0) {
    outliers_score <- sort(c(outliers_reference_score, outliers_alternative_score))
    ind_outliers <- which(sort(score) %in% outliers_score)
    lines(ind_outliers, sort(score)[ind_outliers], col = "red", type = "p")
  }
  dev.off()
  if (num_outliers > 0) {
    png(paste0(sample_dir, "/", sample_name, "_reads_scores_no_outliers.png"))
    score_no_outliers <- sort(score)[-ind_outliers]
    plot(1:length(score_no_outliers), score_no_outliers, xlab = "Reads", ylab = "Score (bp)", main = "Reads scores")
    abline(h = score_thr, col = "blue", lty = 2, lwd = 3)
    dev.off()
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
  
  target_reads_consensus <- TRC
  
  if (skip_first_allele_flag == 1) {
    cat(text = paste0("WARNING: Only ", num_reads_first_allele, " reads (", sprintf("%.2f", allelic_ratio_perc_first), "%) from sample ", sample_name, " have been assigned to Allele #1; skipping"), sep = "\n")
    cat(text = paste0("WARNING: Only ", num_reads_first_allele, " reads (", sprintf("%.2f", allelic_ratio_perc_first), "%) from sample ", sample_name, " have been assigned to Allele #1; skipping"),  file = logfile, sep = "\n", append = TRUE)
  } else {
    if (num_reads_first_allele < target_reads_consensus) {
      target_reads_consensus <- num_reads_first_allele
      target_reads_polishing <- num_reads_first_allele
      cat(text = paste0("WARNING: Only ", num_reads_first_allele, " reads available for sample ", sample_name, " for Allele #1"), sep = "\n")
      cat(text = paste0("WARNING: Only ", num_reads_first_allele, " reads available for sample ", sample_name, " for Allele #1"),  file = logfile, sep = "\n", append = TRUE)
    } 
    plurality_value <- 0.15*target_reads_consensus
    draft_consensus_first_tmp1 <- paste0(sample_dir, "/", sample_name, "_draft_first_allele_tmp1.fasta")
    draft_consensus_first_tmp2 <- paste0(sample_dir, "/", sample_name, "_draft_first_allele_tmp2.fasta")
    draft_consensus_first <- paste0(sample_dir, "/", sample_name, "_draft_first_allele.fasta")
    first_allele_untrimmed <- paste0(analysis_dir, "/", sample_name, "/", sample_name, "_draft_first_allele_untrimmed.fasta")
    first_allele <- paste0(analysis_dir, "/", sample_name, "_first_allele.fasta")
    sequences <- readDNAStringSet(fasta_files[i], "fasta")
    ws <- width(sequences)
    amplicon_length <- ceiling(mean(ws))
    draft_reads_fq_first <- paste0(sample_dir, "/", sample_name, "_draft_", target_reads_consensus, "_reads_first_allele.fastq")
    draft_reads_fa_first <- paste0(sample_dir, "/", sample_name, "_draft_", target_reads_consensus, "_reads_first_allele.fasta")
    seed <- 1
    system(paste0(SEQTK, " sample -s ", seed , " ", first_allele_reads_fq, " ",  target_reads_consensus, " > ", draft_reads_fq_first))
    system(paste0(SEQTK, " seq -A ", draft_reads_fq_first, " > ", draft_reads_fa_first))
    mfa_file_first <- gsub(pattern = "\\.fasta$", replacement = ".mfa", x = draft_reads_fa_first)
    #system(paste0(MAFFT, " --auto --thread ", num_threads, " --threadit 0 --adjustdirectionaccurately ", draft_reads_fa_first, " > ", mfa_file_first)) #threadit 0 if you need reproducible results
    if (fast_alignment_flag == 1) {
      system(paste0(MAFFT, " --auto --thread ", num_threads, " --adjustdirectionaccurately ", draft_reads_fa_first, " > ", mfa_file_first))
    } else {
      system(paste0(MAFFT , "-linsi --thread ", num_threads, " --adjustdirectionaccurately ", draft_reads_fa_first, " > ", mfa_file_first))
    }
    system(paste0(CONS, " -sequence ", mfa_file_first, " -plurality ", plurality_value, " -outseq ", draft_consensus_first_tmp1))
    system(paste0("sed 's/[nN]//g' ", draft_consensus_first_tmp1, " > ", draft_consensus_first_tmp2))
    DNAStringSet_obj <- readDNAStringSet(draft_consensus_first_tmp2, "fasta")
    DNAStringSet_obj_renamed <- DNAStringSet_obj
    original_headers <- names(DNAStringSet_obj)
    sequences <- seq(DNAStringSet_obj)
    names(DNAStringSet_obj_renamed) <- "Allele_number_1"
    writeXStringSet(x = DNAStringSet_obj_renamed, filepath = draft_consensus_first, format = "fasta", width = 20000)
    if (pair_strands_flag != 1) {
      num_threads_medaka <- min(num_threads, 8)
      paf_file_first <- gsub(pattern = "\\.fastq$", replacement = ".paf", x = first_allele_reads_fq)
      racon_consensus_first <- paste0(sample_dir, "/", sample_name, "_racon_first_allele.fasta")
      polishing_reads_fq_first <- paste0(sample_dir, "/", sample_name, "_polishing_", target_reads_polishing, "_reads_first_allele.fastq")
      seed <- 2
      system(paste0(SEQTK, " sample -s ", seed , " ", first_allele_reads_fq, " ",  target_reads_polishing, " > ", polishing_reads_fq_first))
      cat(text = paste0("Running Racon for consensus polishing of sample ", sample_name, " - Allele #1"), sep = "\n")
      system(paste0(MINIMAP2, " -x ava-ont ", draft_consensus_first, " ", polishing_reads_fq_first, " > ", paf_file_first))
      system(paste0(RACON, " -t ", num_threads, " -m 8 -x -6 -g -8 -w 500 --no-trimming ", polishing_reads_fq_first, " ", paf_file_first, " ", draft_consensus_first, " > ", racon_consensus_first))
      cat(text = paste0("Running Medaka for consensus polishing of sample ", sample_name, " - Allele #1"), sep = "\n")
      system(paste0(MEDAKA, "_consensus -i ", polishing_reads_fq_first, " -d ", racon_consensus_first, " -m ", medaka_model, " -t ", num_threads_medaka, " -o ", sample_dir, "/medaka_first_allele"))
      system(paste0("cp ", sample_dir, "/medaka_first_allele/consensus.fasta ", first_allele_untrimmed))
      system(paste0(SEQTK, " trimfq ", first_allele_untrimmed, " -b ", primers_length, " -e ", primers_length, " > ", first_allele))
    } else {
      system(paste0(SEQTK, " trimfq ", draft_consensus_first, " -b ", primers_length, " -e ", primers_length, " > ", first_allele))
    }
    setwd(analysis_dir)
    system(paste0(TRF, " ", first_allele, " 2 7 7 80 10 50 500"))
    system(paste0(TRF, " ", first_allele, " 2 3 5 80 10 14 500"))
    system(paste0(TRF, " ", first_allele, " 2 500 500 80 10 50 500"))
    default_trf_par_files <-  list.files(path = analysis_dir, pattern = "\\.2\\.7\\.7\\.80\\.10\\.50\\.500", full.names = TRUE)
    stringent_trf_par_files <- list.files(path = analysis_dir, pattern = "\\.2\\.500\\.500\\.80\\.10\\.50\\.500", full.names = TRUE)
    lenient_trf_par_files <- list.files(path = analysis_dir, pattern = "\\.2\\.3\\.5\\.80\\.10\\.14\\.500", full.names = TRUE)
    system(paste0("mv ", paste0(lenient_trf_par_files, collapse = " "), " ", sample_dir))
    system(paste0("mv ", paste0(stringent_trf_par_files, collapse = " "), " ", sample_dir))
  }
  #create consensus sequence for Allele #2
  target_reads_consensus <- TRC
  target_reads_polishing <- TRP  

  if (skip_second_allele_flag == 1) {
    cat(text = paste0("WARNING: Only ", num_reads_second_allele, " reads (", sprintf("%.2f", allelic_ratio_perc_second), "%) from sample ", sample_name, " have been assigned to Allele #2; skipping"), sep = "\n")
    cat(text = paste0("WARNING: Only ", num_reads_second_allele, " reads (", sprintf("%.2f", allelic_ratio_perc_second), "%) from sample ", sample_name, " have been assigned to Allele #2; skipping"),  file = logfile, sep = "\n", append = TRUE)
  } else {
    if (num_reads_second_allele < target_reads_consensus) {
      target_reads_consensus <- num_reads_second_allele
      target_reads_polishing <- num_reads_second_allele
      cat(text = paste0("WARNING: Only ", num_reads_second_allele, " reads available for sample ", sample_name, " for Allele #2"), sep = "\n")
      cat(text = paste0("WARNING: Only ", num_reads_second_allele, " reads available for sample ", sample_name, " for Allele #2"),  file = logfile, sep = "\n", append = TRUE)
    } 
    plurality_value <- 0.15*target_reads_consensus
    draft_consensus_second_tmp1 <- paste0(sample_dir, "/", sample_name, "_draft_second_allele_tmp1.fasta")
    draft_consensus_second_tmp2 <- paste0(sample_dir, "/", sample_name, "_draft_second_allele_tmp2.fasta")
    draft_consensus_second <- paste0(sample_dir, "/", sample_name, "_draft_second_allele.fasta")
    second_allele_untrimmed <- paste0(analysis_dir, "/", sample_name, "/", sample_name, "_draft_second_allele_untrimmed.fasta")
    second_allele <- paste0(analysis_dir, "/", sample_name, "_second_allele.fasta")
    sequences <- readDNAStringSet(fasta_files[i], "fasta")
    ws <- width(sequences)
    amplicon_length <- ceiling(mean(ws))
    draft_reads_fq_second <- paste0(sample_dir, "/", sample_name, "_draft_", target_reads_consensus, "_reads_second_allele.fastq")
    draft_reads_fa_second <- paste0(sample_dir, "/", sample_name, "_draft_", target_reads_consensus, "_reads_second_allele.fasta")
    seed <- 1
    system(paste0(SEQTK, " sample -s ", seed , " ", second_allele_reads_fq, " ",  target_reads_consensus, " > ", draft_reads_fq_second))
    system(paste0(SEQTK, " seq -A ", draft_reads_fq_second, " > ", draft_reads_fa_second))
    mfa_file_second <- gsub(pattern = "\\.fasta$", replacement = ".mfa", x = draft_reads_fa_second)
    #system(paste0(MAFFT, " --auto --thread ", num_threads, " --threadit 0 --adjustdirectionaccurately ", draft_reads_fa_second, " > ", mfa_file_second)) #threadit 0 if you need reproducible results
    if (fast_alignment_flag == 1) {
      system(paste0(MAFFT, " --auto --thread ", num_threads, " --adjustdirectionaccurately ", draft_reads_fa_second, " > ", mfa_file_second))
    } else {
      system(paste0(MAFFT, "-linsi --thread ", num_threads, " --adjustdirectionaccurately ", draft_reads_fa_second, " > ", mfa_file_second))
    }
    system(paste0(CONS, " -sequence ", mfa_file_second, " -plurality ", plurality_value, " -outseq ", draft_consensus_second_tmp1))
    system(paste0("sed 's/[nN]//g' ", draft_consensus_second_tmp1, " > ", draft_consensus_second_tmp2))
    DNAStringSet_obj <- readDNAStringSet(draft_consensus_second_tmp2, "fasta")
    DNAStringSet_obj_renamed <- DNAStringSet_obj
    original_headers <- names(DNAStringSet_obj)
    sequences <- seq(DNAStringSet_obj)
    names(DNAStringSet_obj_renamed) <- "Allele_number_2"
    writeXStringSet(x = DNAStringSet_obj_renamed, filepath = draft_consensus_second, format = "fasta", width = 20000)
    if (pair_strands_flag != 1) {
      num_threads_medaka <- min(num_threads, 8)
      paf_file_second <- gsub(pattern = "\\.fastq$", replacement = ".paf", x = second_allele_reads_fq)
      racon_consensus_second <- paste0(sample_dir, "/", sample_name, "_racon_second_allele.fasta")
      polishing_reads_fq_second <- paste0(sample_dir, "/", sample_name, "_polishing_", target_reads_polishing, "_reads_second_allele.fastq")
      seed <- 2
      system(paste0(SEQTK, " sample -s ", seed , " ", second_allele_reads_fq, " ",  target_reads_polishing, " > ", polishing_reads_fq_second))
      cat(text = paste0("Running Racon for consensus polishing of sample ", sample_name, " - Allele #2"), sep = "\n")
      system(paste0(MINIMAP2, " -x ava-ont ", draft_consensus_second, " ", polishing_reads_fq_second, " > ", paf_file_second))
      system(paste0(RACON, " -t ", num_threads, " -m 8 -x -6 -g -8 -w 500 --no-trimming ", polishing_reads_fq_second, " ", paf_file_second, " ", draft_consensus_second, " > ", racon_consensus_second))
      cat(text = paste0("Running Medaka for consensus polishing of sample ", sample_name, " - Allele #2"), sep = "\n")
      system(paste0(MEDAKA, "_consensus -i ", polishing_reads_fq_second, " -d ", racon_consensus_second, " -m ", medaka_model, " -t ", num_threads_medaka, " -o ", sample_dir, "/medaka_second_allele"))
      system(paste0("cp ", sample_dir, "/medaka_second_allele/consensus.fasta ", second_allele_untrimmed))
      system(paste0(SEQTK, " trimfq ", second_allele_untrimmed, " -b ", primers_length, " -e ", primers_length, " > ", second_allele))
    } else {
      system(paste0(SEQTK, " trimfq ", draft_consensus_second, " -b ", primers_length, " -e ", primers_length, " > ", second_allele))
    }
    setwd(analysis_dir)
    system(paste0(TRF, " ", second_allele, " 2 7 7 80 10 50 500"))
    system(paste0(TRF, " ", second_allele, " 2 3 5 80 10 14 500"))
    system(paste0(TRF, " ", second_allele, " 2 500 500 80 10 50 500"))
    default_trf_par_files <-  list.files(path = analysis_dir, pattern = "\\.2\\.7\\.7\\.80\\.10\\.50\\.500", full.names = TRUE)
    stringent_trf_par_files <- list.files(path = analysis_dir, pattern = "\\.2\\.500\\.500\\.80\\.10\\.50\\.500", full.names = TRUE)
    lenient_trf_par_files <- list.files(path = analysis_dir, pattern = "\\.2\\.3\\.5\\.80\\.10\\.14\\.500", full.names = TRUE)
    system(paste0("mv ", paste0(lenient_trf_par_files, collapse = " "), " ", sample_dir))
    system(paste0("mv ", paste0(stringent_trf_par_files, collapse = " "), " ", sample_dir))
  }
}
