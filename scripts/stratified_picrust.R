# Load packages 
library(tidyverse)
library(readr)
library(dplyr)
library(ggplot2)
library(phyloseq)

# Load MetaCyc stratified pathway contribution data
strat_mc_data_1 <- data.table::fread("path_abun_contrib_1.tsv")

norm_calc <- strat_mc_data_1 %>%
  group_by(sample, `function`) %>%
  mutate(norm_taxon_function_contrib = taxon_function_abun / sum(taxon_function_abun)) %>%
ungroup()

test = norm_calc %>%
  filter(sample == "13414.1.gill.R1.fastq.gz", `function` == "1CMET2-PWY")

average_contribution <- norm_calc %>%
  group_by(`function`, taxon) %>%
  mutate(avg_contrib = mean(norm_taxon_function_contrib)) %>%
  ungroup() %>%
  select(`function`, taxon, avg_contrib) %>%
  unique()

sig_contrib <- average_contribution %>%
  filter(avg_contrib >= 0.33)

# Define Energy Pathways of interest
energy_pathways <- c(
  "ANAGLYCOLYSIS-PWY",
  "GLYCOLYSIS",
  "GLYCOLYSIS-E-D",
  "FERMENTATION-PWY",
  "METH-ACETATE-PWY",
  "FAO-PWY",
  "COA-PWY",
  "COA-PWY-1",
  "GLUCONEO-PWY",
  "PENTOSE-P-PWY",
  "NONOXIPENT-PWY",
  "GLYOXYLATE-BYPASS"
)

mc_pathways <- sig_contrib %>%
  filter(`function` %in% energy_pathways) %>%
  unique()


         
# View MetaCyc function names
unique(strat_mc_data$`function`)


# Filter for the energy pathways
filtered_pathways <- strat_mc_data %>%
  filter(`function` %in% energy_pathways)

# Normalize contributions within each (function x sample)
# prop_contrib = proportion of reads to the pathway from each taxon

normalized_pathways <- filtered_pathways %>%
  group_by(`function`, sample) %>%
  mutate(
    prop_contrib = taxon_rel_function_abun /
      sum(taxon_rel_function_abun, na.rm = TRUE)
  ) %>%
  ungroup()

# Group by function and taxon, calculate contribution of reads 
avg_contribution <- normalized_pathways %>%
  group_by(`function`, taxon) %>%
  summarize(
    mean_prop = mean(norm_taxon_function_contrib, na.rm = TRUE),
    .groups = "drop"
  )

# Identify major contributors
# Filter for top contributors (>=33% of reads)
top_contributors <- avg_contribution %>%
  filter(mean_prop >= 0.33, na.rm = TRUE)

# Load Taxonomy Data Names
tax <- read_tsv("data_taxonomy.tsv")

# Merge Taxa names and Feature ID
merged_taxa_names <- top_contributors %>%
  left_join(tax, by = c("taxon" = "Feature ID"))

# See counts of pathways with taxa
taxa_pathway_counts_2 <- merged_taxa_names %>%
  group_by(Taxon) %>%
  summarize(n_pathways = n_distinct(`function`)) %>%
  arrange(desc(n_pathways))

# Filter for n_pathways >0
contributors_only <- taxa_pathway_counts_2 %>%
  filter(n_pathways > 0)

# See how many unique taxa are left 
unique_taxa_count <- length(unique(contributors_only$Taxon))
unique_taxa_count

