#################################################################
# Plots to visualise taxonomic 
# diversity of chemical impacted communities
#
# T. Smith 2025
#################################################################

rm(list = ls())

# load some packages
library(qiime2R)
library(tidyverse)

library(ggplot2)
library(vegan)
library(dplyr)
library(lme4)
library(lmerTest)
#library(ggbeeswarm)

setwd("~/Documents/chems-vs-comms/code/")

### --- Load the data --- ###

# load in the sequencing data (ASV table)
SVs <- read_qza("../data/Earlham/rarefied_table.qza")$data

# load in the metadata
metadata <- read.csv("../data/Earlham/metadata.tsv", sep = "\t") %>%
  rename(SampleID = sample.id)

metadata[metadata$chem.code == "",]$chem.code <- "Control"

# load in the taxonomy
taxonomy <- read_qza("../data/Earlham/taxonomy.qza")$data %>% parse_taxonomy()

# load in the phylogeny
tree <- read_qza("../data/Earlham/rooted_tree.qza")$data

### --- Taxonomy barplots --- ###

# summarise the taxa to the genus level
taxasums <- summarize_taxa(SVs, taxonomy)$Genus

# plot it
#options(repr.plot.width = 20, repr.plot.height = 6) # seetting plot size in jupyter
svg("../results/taxonomy-barplots.svg", width = 20, height = 4)
taxa_barplot(taxasums, metadata, "community.number") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 14),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 6),
        legend.position = "bottom")
dev.off()

# also do one for oxytetracycline, which we can cleverly piece together into our figure later

svg("../results/taxonomy-barplots-oxytet.svg", width = 20, height = 4)
taxa_barplot(taxasums, metadata, "Oxytetracycline") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 14),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 6),
        legend.position = "bottom")
dev.off()


### ------ Diversity Metrics ------ ###

# transform matrix for vegan community metrics
rarefied_OTUs_t <- t(SVs)

diversity_shannon <- diversity(rarefied_OTUs_t)
diversity_simpson <- diversity(rarefied_OTUs_t, "simpson")
diversity_invsimp <- diversity(rarefied_OTUs_t, "inv")
diversity_alpha <- fisher.alpha(rarefied_OTUs_t)

species_richness <- specnumber(rarefied_OTUs_t)
evenness <- diversity_shannon/log(species_richness)

diversity_metrics <- data.frame(diversity_shannon, diversity_simpson, diversity_invsimp, diversity_alpha, species_richness,
                                evenness)
diversity_metrics$SampleID <- rownames(diversity_metrics)

# add the metadata back
diversity_metrics <- left_join(diversity_metrics, metadata)

# do some plotting
svg("../results/diversity-richness.svg", width = 4, height = 4)
ggplot(diversity_metrics, aes(x = complexity, y = species_richness, fill = as.factor(Oxytetracycline))) +
  geom_beeswarm(shape = 21) +
  scale_fill_manual(values = c("black", "orange")) +
  geom_smooth(method = lm, col = "black") +
  labs(x = "Number of chemicals", y = "Species Richness") +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.7, 0.8),
        legend.title = element_blank(),
        legend.text = element_blank(),
        aspect.ratio = 1)
dev.off()

summary(lmer(scale(species_richness) ~ complexity + Oxytetracycline + (1|community.name), data = diversity_metrics))

svg("../results/diversity-shannon.svg", width = 4, height = 4)
ggplot(diversity_metrics, aes(x = complexity, y = diversity_shannon, fill = as.factor(Oxytetracycline))) +
  geom_beeswarm(shape = 21) +
  scale_fill_manual(values = c("black", "orange")) +
  geom_smooth(method = lm, col = "black") +
  labs(x = "Number of chemicals", y = "Shannon Diversity") +
  theme_bw() +
  theme(legend.position = "none",
        aspect.ratio = 1)
dev.off()

summary(lmer(scale(diversity_shannon) ~ complexity + Oxytetracycline + (1|community.name), data = diversity_metrics))

### --- Functional diversity (enzymatic pathways)--- ###

### Enzymes first ###

# read in the picrust data
library(data.table)
# abundance of ECs in metagenome
EC_abundance <- as.data.frame(fread("../data/Earlham/picrust/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz"))

# transform these data into the sort of matrix that vegan likes
EC_matrix <- data.matrix(EC_abundance[,3:376])
dimnames(EC_matrix)[[1]] <- EC_abundance$description
EC_matrix_t <- t(EC_matrix)

# diversity metrics
EC_richness <- specnumber(EC_matrix_t)
EC_shannon <- diversity(EC_matrix_t)
EC_simpson <- diversity(EC_matrix_t, "simpson")
EC_evenness <- EC_shannon/log(EC_richness)

EC_diversity_metrics <- data.frame(EC_richness, EC_shannon, EC_simpson, EC_evenness)
EC_diversity_metrics$SampleID <- rownames(EC_diversity_metrics)

# add the metadata back
EC_diversity_metrics <- left_join(EC_diversity_metrics, metadata) %>%
  filter(chem.code != "Frozen Starting Community")

svg("../results/diversity-EC-richness.svg", width = 4, height = 4)
ggplot(EC_diversity_metrics, aes(x = complexity, y = EC_richness, fill = as.factor(Oxytetracycline))) +
  geom_beeswarm(shape = 21) +
  scale_fill_manual(values = c("black", "orange")) +
  geom_smooth(method = lm, col = "black") +
  labs(x = "Number of chemicals", y = "Enzyme Richness") +
  theme_bw() +
  theme(legend.position = "none",
        aspect.ratio = 1)
dev.off()

summary(lmer(scale(EC_richness) ~ complexity + Oxytetracycline + (1|community.name), data = EC_diversity_metrics))

svg("../results/diversity-EC-shannon.svg", width = 4, height = 4)
ggplot(EC_diversity_metrics, aes(x = complexity, y = EC_shannon, fill = as.factor(Oxytetracycline))) +
  geom_beeswarm(shape = 21) +
  scale_fill_manual(values = c("black", "orange")) +
  geom_smooth(method = lm, col = "black") +
  labs(x = "Number of chemicals", y = "Enzyme Shannon Diversity") +
  theme_bw() +
  theme(legend.position = "none",
        aspect.ratio = 1)

summary(lmer(scale(EC_shannon) ~ complexity + Oxytetracycline + (1|community.name), data = EC_diversity_metrics))

### now MetaCyc pathways ###
# abundance of MetaCyc pathways in metagenome
MC_abundance <- as.data.frame(fread("../data/Earlham/picrust/pathways_out/path_abun_unstrat_descrip.tsv.gz"))

MC_matrix <- data.matrix(MC_abundance[,3:375])
dimnames(MC_matrix)[[1]] <- MC_abundance$description
MC_matrix_t <- t(MC_matrix)

MC_richness <- specnumber(MC_matrix_t)
MC_shannon <- diversity(MC_matrix_t)
MC_simpson <- diversity(MC_matrix_t, "simpson")

MC_evenness <- MC_shannon/log(MC_richness)

MC_diversity_metrics <- data.frame(MC_richness, MC_shannon, MC_simpson, MC_evenness)
MC_diversity_metrics$SampleID <- rownames(MC_diversity_metrics)
MC_diversity_metrics <- left_join(MC_diversity_metrics, metadata) %>%
  filter(chem.code != "Frozen Starting Community")

ggplot(MC_diversity_metrics, aes(x = complexity, y = MC_richness, fill = as.factor(Oxytetracycline))) +
  geom_beeswarm(shape = 21) +
  scale_fill_manual(values = c("black", "orange")) +
  geom_smooth(method = lm, col = "black") +
  labs(x = "Number of chemicals", y = "Pathway Richness") +
  theme_bw() +
  theme(legend.position = "none",
        aspect.ratio = 1)

summary(lmer(scale(MC_richness) ~ complexity + Oxytetracycline + (1|community.name), data = MC_diversity_metrics))

ggplot(MC_diversity_metrics, aes(x = complexity, y = MC_shannon, fill = as.factor(Oxytetracycline))) +
  geom_beeswarm(shape = 21) +
  scale_fill_manual(values = c("black", "orange")) +
  geom_smooth(method = lm, col = "black") +
  labs(x = "Number of chemicals", y = "Pathway Shannon Diversity") +
  theme_bw() +
  theme(legend.position = "none",
        aspect.ratio = 1)

summary(lmer(scale(MC_shannon) ~ complexity + Oxytetracycline + (1|community.name), data = MC_diversity_metrics))
