# Differential abundance with ANCOM-BC2
library(tidyverse)
library(phyloseq)
library(readxl)

# Load object
ps = readRDS('ps_filtered.rds')

library(ANCOMBC)

# Tax_glom, convert up to the Genus level
ps_glom = tax_glom(ps, 'Genus')

# Filter to hindgut samples only 
ps_genus_final = subset_samples(ps_glom, sample_type == "hindgut")

# Subset your phyloseq object by swim performance group
performance.fast = subset_samples(ps_genus_final, swim_performance %in% c('accelerator', 'crusier sprinter', 'manoeuvrer'))
performance.slow = subset_samples(ps_genus_final, swim_performance %in% c('flow refuging', 'burrowing'))
performance.moderate = subset_samples(ps_genus_final, swim_performance == 'generalist')

# Set random seed
set.seed(421)

# Run ANCOMBC differential abundance analysis
out = ancombc2(data = ps_genus_final, 
               fix_formula = 'swim_performance', 
               p_adj_method = 'BH',
               prv_cut = 0.1)

# Extract stats table from ANCOMBC 
statistical_table = out$res

# Write ANCOMBC results into an Excel file
writexl::write_xlsx(statistical_table,'ANCOMBC Results.xlsx')

# Read ANCOMBC results
ANCOMBC_results <- read_excel("ANCOMBC Results.xlsx")

# Important taxonomy table containing taxa names
tax <- read_tsv("data_taxonomy.tsv")

# Merge ANCOMBC feature IDs with taxa names
merged_ancombc <- ANCOMBC_results %>%
  left_join(tax, by = c("taxon" = "Feature ID"))
