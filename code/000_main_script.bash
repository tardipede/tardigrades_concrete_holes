#!/usr/bin/bash

# Process raw reads into ASV tables

Rscript code/001_coi_reads_process.R --input "./data"  \
                                           --intermediates "./intermediates" \
                                           --output "./results" \
                                           --run_name "COI" \
                                           --run_regex "(ID[0-9]{4})" \
                                           --sample_regex "((BL|IT)-[0-9]{3})" \
                                           --database_file "./databases/db_all.fas" \
                                           --forward_primer "GCNCCNGAYATRKSNTTYCC" \
                                           --reverse_primer "TCDGGRTGNCCRAARAAYCA" \
                                           --min_amplicon_size 415 \
                                           --max_amplicon_size 427 \
                                           --amplicon_size_step 3


INPUT_FOLDER = "./data"
INTERMEDIATE_FOLDER = "./intermediates"
OUTPUT_FOLDER = "./results"
RUN_NAME = "COI"
FORWARD_READS_REGEX = ".*R1(_001)?.fastq.gz"
REVERSE_READS_REGEX = ".*R2(_001)?.fastq.gz"
RUN_REGEX = "(ID[0-9]{4})" 
SAMPLE_REGEX = "((BL|IT)-[0-9]{3})"
DATABASE_PATH = "./databases/db_all.fas"
FWD = "GCNCCNGAYATRKSNTTYCC"
REV = "TCDGGRTGNCCRAARAAYCA"
min_amplicon_size = 415
max_amplicon_size = 427
amplicon_size_step = 3
PROCESSORS = 12