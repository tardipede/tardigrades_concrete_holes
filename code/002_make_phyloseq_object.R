#!/usr/bin/env Rscript

## Script for creating phyloseq objects

# Load argparse library to parse arguments
suppressPackageStartupMessages({
    library(argparse)
    })

## Main program starts here
###############################################################################

# Arguments parsing
parser <- ArgumentParser()
parser$add_argument("--otu_table", help = "path to OTU or ASV table with taxonomy")
parser$add_argument("--sequences", help = "path to aligned sequences")
parser$add_argument("--samples_data", help = "path to samples metadata")
parser$add_argument("--run_name", help = "run name")
parser$add_argument("--output", help = "path to output folder")

arguments <- parser$parse_args()

OTU_TABLE = arguments$otu_table
SEQS = arguments$sequences
SAMPLES_DATA = arguments$samples_data
RUN_NAME = arguments$run_name
OUTPUT_FOLDER = arguments$output

# Load libraries
suppressPackageStartupMessages({
    library(ape)
    library(Biostrings)
    library(magrittr)
    library(phyloseq) 
    })

# Load files
otu_tab <- read.table(OTU_TABLE, 
                        sep = "\t", 
                        header = TRUE)
rownames(otu_tab) <- otu_tab$otu

sample_dat <- read.table(SAMPLES_DATA, 
                           sep = "\t", 
                           header = TRUE)
rownames(sample_dat) <- sample_dat$sample

seqs <- Biostrings::readDNAStringSet(SEQS)


# Prepare the components to build a phyloseq object

## OTU table
phyloseq_otu_tab <- phyloseq::otu_table(otu_tab[, colnames(otu_tab) %in% sample_dat$sample], 
                                        taxa_are_rows = TRUE)

## Taxonomy table
### Extract taxonomy from otutable
taxonomy_table <- otu_tab[,c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species")]
taxonomy_table <- as.matrix(taxonomy_table)

phyloseq_tax_tab <- phyloseq::tax_table(taxonomy_table)

## Sample data
phyloseq_sample_data <- phyloseq::sample_data(sample_dat)

## Sequences
phyloseq_seqs <- phyloseq::refseq(seqs)

# Assembly phyoseq object
dataset <- phyloseq(phyloseq_otu_tab, 
                    phyloseq_sample_data, 
                    phyloseq_tax_tab, 
                    phyloseq_seqs)

# Create rds filename
rds_filename = paste0(RUN_NAME,"_phyloseq.rds")
saveRDS(dataset, file.path(OUTPUT_FOLDER,rds_filename))