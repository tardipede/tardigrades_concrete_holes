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
                                           --amplicon_size_step 3 \
                                           --translation_frame 2

# Create phyloseq object
Rscript code/002_make_phyloseq_object.R --otu_table "./results/COI_classified_otutab.tsv"  \
                                        --sequences "./results/COI_OTU.fasta"  \
                                        --samples_data "./samples_data.txt"  \
                                        --run_name "COI" \
                                        --output "./results"

# Run analysis and visualization
Rscript code/003_analysis_and_visualization.R