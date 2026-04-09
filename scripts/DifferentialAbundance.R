# This differential abundance analysis is performed using Maaslin2 and was used to see if there are any differential abundant taxa in Slow and Fast groups. 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Differential abundance with Maaslin2
library(tidyverse)
library(phyloseq)
library(readxl)
 if(!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
 BiocManager::install("Maaslin2")
library(Maaslin2)

# Load object
ps = readRDS('../datasets/ps_filtered.rds')

# Tax_glom if desired. DON'T rarefy or use a compositional transformation!
ps_glom = tax_glom(ps, 'Genus')

# Filter to hindgut samples only 
ps_genus_final = subset_samples(ps_glom, sample_type == "hindgut")

# Split into fast / moderate / slow
sample_data(ps_genus_final)$swim_group <- 
  fct_collapse(sample_data(ps_genus_final)$swim_performance, 
               fast =c("accelerator", "crusier sprinter", "manoeuvrer"), 
               moderate=c("generalist"), 
               slow =c("flow refuging", "burrowing"), 
               )

set.seed(421)
out = Maaslin2(
  # Your count table and metadata, extracted into data frames
  input_data = data.frame(ps_genus_final@otu_table), 
  input_metadata = data.frame(ps_genus_final@sam_data), 
  # Maaslin2 forces you to specify an output directory, where it automatically saves some outputs. I usually just delete these results because I don't like the formatting.
  output = 'to_delete', 
  # List out your explanatory variable(s) - ex. c('var1','var2')
  fixed_effects = c('swim_group'),
  # Optional: specify a reference group for your categorical variable(s). 
  # This will be the group that all other groups are compared to. 
  # Ex. if you have a variable called 'Treatment' with two groups 'Control' and 'Intervention', you would set reference = 'Treatment,Control' to assess Intervention relative to Control.
  reference = 'swim_group,moderate', # We used the moderate group as a reference group so we can make a fair comparison between slow and fast.
  # TSS is Total Sum Scaled, AKA relative abundance.
  normalization = 'TSS',
  # Arcsine transformation is recommended for relative abundance data.
  transform = 'AST',
  # Set abundance, prevalence, Padj value thresholds
  min_abundance = 0.001, # 0.001 means 0.1% relative abundance
  min_prevalence = 0.1, # 0.1 means 10% of samples
  max_significance = 0.05, # Padj threshold
  # These plots increase run time and are usually poorly formatted, so you might as well stop them from being generated
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)

statistical_table = out$results

# Export results table
writexl::write_xlsx(statistical_table,'Module 14 - Maaslin2 Results.xlsx')

# Loading in results as an object, then converting feature names to taxa
maaslin_results <- read_excel('Module 14 - Maaslin2 Results.xlsx')
maaslin_results$feature <- sub("^X", "", maaslin_results$feature)
tax <- read_tsv("data_taxonomy.tsv")

merged_maaslin <- maaslin_results %>%
  left_join(tax, by = c("feature" = "Feature ID"))
