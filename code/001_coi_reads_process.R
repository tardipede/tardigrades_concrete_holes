#!/usr/bin/env Rscript

## Script for processing paied reads (Illumina - AVITI) into ASV tables

# Load argparse library to parse arguments
suppressPackageStartupMessages({
  library(argparse)
  })
## Main program starts here
###############################################################################

# Arguments parsing
parser <- ArgumentParser()
parser$add_argument("-i", "--input", help = "input directory path where raw reads are stored")
parser$add_argument("-m", "--intermediates", help = "intermediate files directory path")
parser$add_argument("-o", "--output", help = "output directory path")
parser$add_argument("--run_name", help = "run name (prefix for output files)")
parser$add_argument("--forward_reads_regex", default = ".*R1(_001)?.fastq.gz", help = "regex to find forward reads files")
parser$add_argument("--reverse_reads_regex", default = ".*R2(_001)?.fastq.gz", help = "regex to find reverse reads files")
parser$add_argument("--run_regex", help = "regex to extract run code from raw reads filenames")
parser$add_argument("--sample_regex", help = "regex to extract sample code from raw reads filenames")
parser$add_argument("--database_file", help = "database for classification")
parser$add_argument("--forward_primer", help = "Forward primer")
parser$add_argument("--reverse_primer", help = "Reverse primer")
parser$add_argument("--min_amplicon_size", help = "Minimum amplicon size to keep")
parser$add_argument("--max_amplicon_size", help = "Maximum amplicon size to keep")
parser$add_argument("--processors", default = 12, help = "Number of processors for ASV classification and alignment")
parser$add_argument("--clustering_treshold", default = 0.03, help = "Number of processors for ASV classification and alignment")

arguments <- parser$parse_args()

INPUT_FOLDER = arguments$input
INTERMEDIATE_FOLDER = arguments$intermediates
OUTPUT_FOLDER = arguments$output
RUN_NAME = arguments$run_name
FORWARD_READS_REGEX = arguments$forward_reads_regex
REVERSE_READS_REGEX = arguments$reverse_reads_regex
RUN_REGEX = arguments$run_regex
SAMPLE_REGEX = arguments$sample_regex
DATABASE_PATH = arguments$database_file
FWD = arguments$forward_primer
REV = arguments$reverse_primer
min_amplicon_size = arguments$min_amplicon_size
max_amplicon_size = arguments$max_amplicon_size
PROCESSORS = arguments$processors
CLUSTERING_TRESHOLD = arguments$clustering_treshold

# Load libraryes and custom scripts
suppressPackageStartupMessages({
  options(tidyverse.quiet = TRUE)
  library(tidyverse, quietly = TRUE)
  library(magrittr, quietly = TRUE)
  library(dada2, quietly = TRUE)
  library(ShortRead, quietly = TRUE)
  library(Biostrings, quietly = TRUE)
  library(DECIPHER, quietly = TRUE)
  source("./code/999_utils.R")
  })

#### Get the samples organized ####

# define paths
raw_path = file.path(INPUT_FOLDER)
trim_path = file.path(INTERMEDIATE_FOLDER, "01_trim")
filt_path = file.path(INTERMEDIATE_FOLDER, "02_filt")
backup_path = file.path(INTERMEDIATE_FOLDER, "backups")
results_path = file.path(OUTPUT_FOLDER)

# create paths if missing
if (!dir.exists(filt_path)) dir.create(filt_path, recursive = TRUE)
if (!dir.exists(trim_path)) dir.create(trim_path, recursive = TRUE)
if (!dir.exists(backup_path)) dir.create(backup_path, recursive = TRUE)
if (!dir.exists(results_path)) dir.create(results_path, recursive = TRUE)


# find files
sample_table <- tibble::tibble(
  fastq_R1 = sort(list.files(raw_path, FORWARD_READS_REGEX, 
                             full.names = TRUE)),
  fastq_R2 = sort(list.files(raw_path, REVERSE_READS_REGEX, 
                             full.names = TRUE))
) %>% dplyr::mutate(
  runcode = stringr::str_extract(fastq_R1, RUN_REGEX), # Extract run code
  sample = stringr::str_extract(fastq_R1, SAMPLE_REGEX) # Extract sample code
)

# generate filenames for trimmed and filtered reads
sample_table <- sample_table %>% dplyr::mutate(
  trim_R1 = file.path(trim_path, paste(runcode, 
                                       sample, 
                                       "R1_trim.fastq.gz", 
                                       sep = "_")),
  trim_R2 = file.path(trim_path, paste(runcode, 
                                       sample, 
                                       "R2_trim.fastq.gz", 
                                       sep = "_")),
  filt_R1 = file.path(filt_path, paste(runcode, 
                                       sample, 
                                       "R1_filt.fastq.gz", 
                                       sep = "_")),
  filt_R2 = file.path(filt_path, paste(runcode, 
                                       sample, 
                                       "R2_filt.fastq.gz", 
                                       sep = "_"))
)

# Check that there are no duplicates
assertthat::assert_that(
  !any(duplicated(sample_table$fastq_R1)),
  !any(duplicated(sample_table$fastq_R2)),
  !any(duplicated(sample_table$trim_R1)),
  !any(duplicated(sample_table$trim_R2)),
  !any(duplicated(sample_table$filt_R1)),
  !any(duplicated(sample_table$filt_R2))
)

write.table(sample_table, 
            file = paste0(backup_path,"/filenames_table.txt"), 
            sep = "\t", 
            quote = FALSE, 
            col.names = TRUE)

cat("Found", nrow(sample_table), "samples.\n")


# Flags for Cutadapt and trimming primers, filtering reads without primers
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

CUTADAPT_FLAGS <- paste0("cutadapt --cores ",PROCESSORS," -a ^", 
                         FWD, 
                         "...", 
                         REV.RC, 
                         " -A ^", 
                         REV, 
                         "...", 
                         FWD.RC, 
                         " --discard-untrimmed -o")

for (i in seq_along(sample_table$fastq_R1)) {
  system(paste(CUTADAPT_FLAGS, 
               sample_table$trim_R1[i],
               "-p", 
               sample_table$trim_R2[i], 
               sample_table$fastq_R1[i], 
               sample_table$fastq_R2[i]))
}

# Quality-filtering of reads
out <- filterAndTrim(sample_table$trim_R1, 
                     sample_table$filt_R1,
                     sample_table$trim_R2, 
                     sample_table$filt_R2,
  maxN = 0, truncQ = 2, rm.phix = TRUE,
  compress = TRUE, multithread = PROCESSORS, minLen = 190
)

# Error-rate estimation, merging and chimera filtering done separ. for each run
error_results <- list()
backups <- list()
for (i in 1:length(unique(sample_table$runcode))) {
  sample_table_temp <- subset(sample_table, 
                              runcode == unique(sample_table$runcode)[i])

  ## calculating errors - on a subset of data to save time 
  ## see(https://benjjneb.github.io/dada2/bigdata.html)
  set.seed(100)
  errF <- learnErrors(sample_table_temp$filt_R1, 
                      multithread = PROCESSORS, nbases = 1e+8, randomize = TRUE)
  errR <- learnErrors(sample_table_temp$filt_R2, 
                      multithread = PROCESSORS, nbases = 1e+8, randomize = TRUE)
  # Sample inference
  dadaFs <- dada(sample_table_temp$filt_R1, err = errF, multithread = PROCESSORS)
  dadaRs <- dada(sample_table_temp$filt_R2, err = errR, multithread = PROCESSORS)
  # Merging paired reads
  mergers <- mergePairs(dadaFs, 
                        sample_table_temp$filt_R1, 
                        dadaRs, 
                        sample_table_temp$filt_R2, 
                        minOverlap = 8)
  # Creating sequence table
  seqtab <- makeSequenceTable(mergers)
  # Removing chimeras
  seqtab.nochim <- removeBimeraDenovo(seqtab, 
                                      method = "consensus", 
                                      multithread = PROCESSORS, 
                                      verbose = TRUE)
  error_results[[i]] <- seqtab.nochim
  backups[[i]] <- list(errF, errR, dadaFs, dadaRs, mergers, seqtab)

  rm(sample_table_temp)
}

# Merge the results from different runs together
if (length(error_results) == 1) {
  seqtab.nochim <- error_results[[1]]
}
if (length(error_results) != 1) {
  seqtab.nochim <- dada2::mergeSequenceTables(tables = error_results)}


# Save the sequences table
saveRDS(seqtab.nochim, file.path(backup_path, "savepoint1_seqtable.rds"))
# Read the sequence table (only a part of the script can be run from here)
seqtab.nochim <- readRDS(file.path(backup_path, "savepoint1_seqtable.rds"))


# Abundance filtering, length filtering
# This step removes lots of erronneous sequences and saves a lot of RAM to the next steps
seqtab.nochim.filtered_length <- 
  seqtab.nochim[, nchar(colnames(seqtab.nochim)) >= min_amplicon_size & nchar(colnames(seqtab.nochim)) <= max_amplicon_size]
seqtab.nochim.filtered_length.filtered_min_abundance <- 
  seqtab.nochim.filtered_length[, apply(seqtab.nochim.filtered_length, 2, sum) >= 10] # filtering reads with counts < 10
# sorting by abundance
seqtab.nochim.filtered_length.filtered_min_abundance <- 
  seqtab.nochim.filtered_length.filtered_min_abundance[, order(-colSums(seqtab.nochim.filtered_length.filtered_min_abundance))]


# Clustering
seqtab_clustered = cluster_seqs(seqtab.nochim.filtered_length.filtered_min_abundance, 
                                  nproc = PROCESSORS, 
                                  cutoff = CLUSTERING_TRESHOLD, 
                                  align = TRUE)










# ASV classification
trainingSet <- readRDS(DATABASE_PATH) # database loading
ids <- DECIPHER::IdTaxa(dna, trainingSet, strand = "both", processors = 12)

# exporting classification results
output <- sapply(
  ids,
  function(id) {
    paste(id$taxon,
      " (",
      round(id$confidence, digits = 1),
      "%)",
      sep = "",
      collapse = "; "
    )
  }
)

otu_seqtab <- curated_result

colnames(otu_seqtab) <- gsub("-", 
                             ".", 
                             sample_table$sample[match(colnames(otu_seqtab), 
                                                       basename(sample_table$filt_R1))], 
                             fixed = TRUE)

otu_seqtab <- cbind(data.frame(ASV = rownames(otu_seqtab)), 
                    otu_seqtab)
output <- cbind(data.frame(ASV = names(output)), 
                data.frame(output = output))
otutable_classified <- merge(output, 
                             otu_seqtab, 
                             by = "ASV")

otutable_classified <- otutable_classified %>% 
  select(output, everything())

write.table(otutable_classified, 
            file = paste0(OUTPUT_FOLDER,"/", RUN_NAME, "_classified_otus.csv"), 
            quote = FALSE, 
            sep = "\t")

            # Align OTUs and save alignment
aln <- DECIPHER::AlignSeqs(dna, processors = PROCESSORS)
Biostrings::writeXStringSet(aln, 
                            paste0(OUTPUT_FOLDER,"/", 
                                   RUN_NAME, "_aligned_otus.fas"))


# Make OTUs phylogenetic tree with Fastree
system(paste0("FastTree -gtr -nt ",OUTPUT_FOLDER,"/", 
              RUN_NAME, 
              "_aligned_otus.fas > ",OUTPUT_FOLDER, "/", 
              RUN_NAME, 
              "_otus_phylotree.tre"))


