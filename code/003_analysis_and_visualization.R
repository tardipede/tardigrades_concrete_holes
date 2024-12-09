#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages({
    library(phyloseq)
    library(magrittr)
    options(tidyverse.quiet = TRUE)
    library(tidyverse)
    library(patchwork)
    source("./code/999_utils.R")
    })


# Load dataset, remove blanks and otus in samples with < 0.1% relative abundance
dataset = readRDS("results/COI_phyloseq.rds")

dataset_clean = dataset %>% 
  subtract_blank("BL.001") %>%
  prune_samples(sample_sums(.) != 0, .) %>%
  saveRDS_pipe("results/COI_phyloseq_noBlank.rds")

dataset_temp = transform_sample_counts(dataset_clean, function(x) x/sum(x))
otu_table(dataset_clean)[otu_table(dataset_temp) <= 1/1000] = 0
rm(dataset_temp)


### Tardigrada proportion plot
dataset_tardigrades = dataset_clean
tax_table(dataset_tardigrades)[tax_table(dataset_tardigrades)[,2] != "Tardigrada",2] = "Others"
tax_table(dataset_tardigrades) = tax_table(dataset_tardigrades)[,2:ncol(tax_table(dataset_tardigrades))]
dataset_tardigrades = tax_glom(dataset_tardigrades, taxrank = "Phylum")
dataset_tardigrades = transform_sample_counts(dataset_tardigrades, function(x) x/sum(x))

t_density = data.frame(sample_data(dataset_tardigrades))[,c("sample","T_density")]
t_density$T_density = log10(t_density$T_density + 1)
colnames(t_density)[2] = "Tardigrada"
t_density$Tardigrada = scales::rescale(t_density$Tardigrada, from = c(0,4), to = c(0,1))
t_density$Others = 1- t_density$Tardigrada
t_density = pivot_longer(t_density, 2:ncol(t_density))
t_density$type = rep("density", nrow(t_density))

dataset_tardigrades = cbind(data.frame(tax = tax_table(dataset_tardigrades)[,1]),data.frame(otu_table(dataset_tardigrades)))
dataset_tardigrades = pivot_longer(dataset_tardigrades, 2:ncol(dataset_tardigrades))
dataset_tardigrades$type = rep("a_tax", nrow(dataset_tardigrades))

dataset_tardigrades = dataset_tardigrades[,c(2,1,3,4)]
colnames(dataset_tardigrades) = colnames(t_density)
dataset_tardigrades = rbind(dataset_tardigrades, t_density)

plot_t1 = ggplot(dataset_tardigrades)+
  theme_bw()+
  geom_bar(stat = "identity", aes(x = as.numeric(interaction(type, sample)), y = value, fill = name), position = "stack", color = "black", alpha = 0.75)+
  xlab("")+ylab("Reads proportion")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous( name = "Reads proportion", sec.axis = sec_axis( trans=~.*4, name="Tardigrades/gr dry sedimens"))+
  scale_fill_viridis_d()


#### Tardigrada OTUs plot
dataset_tardiotu = subset_taxa(dataset_clean, Phylum == "Tardigrada")
dataset_tardiotu = tax_glom(dataset_tardiotu, "Species", NArm= FALSE)

  
otu_labels = data.frame(tax_table(dataset_tardiotu))
otu_labels[otu_labels == ""] = NA
otu_labels = apply(otu_labels, 1, function(x) tail(na.omit(x), 1))
otu_labels = paste0(names(otu_labels), "_", otu_labels)

dataset_tardiotu = cbind(data.frame(otu = otu_labels,data.frame(otu_table(dataset_tardiotu))))
dataset_tardiotu = dataset_tardiotu[rowSums(dataset_tardiotu[,2:ncol(dataset_tardiotu)]) !=0,]
dataset_tardiotu = pivot_longer(dataset_tardiotu, 2:ncol(dataset_tardiotu))

sumtable = aggregate(dataset_tardiotu$value, by = list(dataset_tardiotu$otu), FUN = sum)
sumtable = sumtable[order(sumtable$x, decreasing = T),]

dataset_tardiotu$otu = factor(dataset_tardiotu$otu, levels = sumtable$Group.1, ordered = TRUE)

plot_t2 = ggplot(dataset_tardiotu)+
  theme_bw()+
  geom_bar(stat = "identity", aes(x = name, y = value, fill = otu), position = "fill", color = "black", alpha = 0.75)+
  xlab("")+ylab("Reads proportion")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_viridis_d(option = "turbo")


# Merge and save plots
plot_merged = plot_t1/plot_t2
ggsave("./results/barplots_concrete.pdf", plot = plot_merged, height = 10, width = 20)

# Save intermediate tables just in case
write.table(dataset_tardigrades, "./results/tardigrade_reads_proportion.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(dataset_tardiotu,"./results/tardigrade_otus_proportion.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Get number of reads per sample
reads_all = colSums(otu_table(dataset))
write.table(data.frame(reads_all),"./results/all_reads_count.tsv", sep = "\t", quote = FALSE)

reads_clean = colSums(otu_table(dataset_clean))
write.table(data.frame(reads_clean),"./results/clean_reads_count.tsv", sep = "\t", quote = FALSE)