#################################################################
# Ploting principal coordinates analysis
# and performing PERMANOVA
#
# T. Smith 2025
#################################################################

rm(list = ls())

setwd("~/Documents/chems-vs-comms/code/")

# load some useful packages
library(qiime2R)
library(tidyverse)
library(dplyr)
library(vegan)
library(ggrepel)

# load in the sequencing data (ASV table)
SVs <- read_qza("../data/Earlham/rarefied_table.qza")$data

# transform matrix for vegan
rarefied_OTUs_t <- t(SVs)

# load in the metadata
metadata <- read.csv("../data/Earlham/metadata.tsv", sep = "\t") %>%
  rename(SampleID = sample.id)

metadata[metadata$chem.code == "",]$chem.code <- "Control"

# load in the taxonomy
taxonomy <- read_qza("../data/Earlham/taxonomy.qza")$data %>% parse_taxonomy()

# load in the phylogeny
tree <- read_qza("../data/Earlham/rooted_tree.qza")$data

# create a plotting theme so we can actually read the axes
main_theme <- theme_bw() + theme(axis.text = element_text(size = 16),
                                 axis.title = element_text(size = 20),
                                 strip.text = element_text(size = 16),
                                 legend.title = element_text(size = 16),
                                 legend.text = element_text(size = 16))

# try PCoA
pco <- wcmdscale(vegdist(rarefied_OTUs_t), eig = TRUE)

## get PCoA scores
scrs <- vegan::scores(pco, choices = 1:2)

# make it a datafrane
pco.df <- data.frame(scrs)
pco.df$SampleID <- rownames(pco.df)

# re-add the metadata
pco.df <- left_join(pco.df, metadata)

# plot PCoA
svg("../results/PCoA-plot.svg", width = 8, height = 6)
ggplot(pco.df %>% filter(chem.code != "Frozen Starting Community"), aes(x = Dim1, y = Dim2, fill = as.character(community.number))) + 
  geom_point(aes(shape = as.character(Oxytetracycline)), size = 3) +
  scale_shape_manual(values = c(21, 24)) + 
  scale_fill_viridis_d() +
  labs(x = "PCo 1", y = "PCo 2", fill = "Community", shape = "Oxytetracycline") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  main_theme +
  theme(aspect.ratio = 1)
dev.off()

# supplementary material, plot with chemicals labelled and communities separated
svg("../results/PCoA-by-community.svg", width = 24, height = 16)
ggplot(pco.df %>% filter(chem.code != "Frozen Starting Community"), aes(x = Dim1, y = Dim2, fill = as.character(complexity))) + 
  geom_point(shape = 21, size = 3) +
  geom_label_repel(aes(label = chem.code),
                   segment.color = 'grey50',
                   max.overlaps = 20) +
  scale_fill_viridis_d() +
  labs(x = "PCo 1", y = "PCo 2", fill = "# Chemicals") +
  facet_wrap(~community.number, scales = "free") +
  main_theme +
  theme(aspect.ratio = 1)
dev.off()

# As expected, replicate treatments largely produce the same outcome
# Oxytetracycline has the clearest impact on community composition
# But it's not consistent between communities - they are shifted in different
# directions in the trait space (and, not all towards the middle).

### ---- PERMANOVA ---- ###

# need the OTU table with OTUs as columns, samples as rows.
# then need the metadata table in the same order with samples as rows.

# need to chop it down to not include those frozen starting community sequences.
rarefied_metadata <- metadata %>% filter(SampleID %in% rownames(rarefied_OTUs_t),
                                         chem.code != "Frozen Starting Community")

rarefied_OTUs_t <- rarefied_OTUs_t[rarefied_metadata$SampleID,]

# rownames(rarefied_OTUs_t) == rarefied_metadata$SampleID
# sorted

# now try the adonis permanova
# number of chemicals
basic.model <- adonis2(rarefied_OTUs_t ~ community.number + complexity,
        data = rarefied_metadata, method = "bray", by = "terms")
basic.model

interactions.model <- adonis2(rarefied_OTUs_t ~ community.name * (Amoxicillin + Chlorothalonil + Diflufenican +
                                              Glyphosate + Imidacloprid + Metaldehyde + Oxytetracycline + Tebuconazole),
        data = rarefied_metadata, method = "bray", by = "terms")
interactions.model
