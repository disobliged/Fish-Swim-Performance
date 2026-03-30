### This R script will perform a core microbiome analysis of fish hindgut samples by assigned swim performance group (Fast, Moderate, Slow)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

# Load object/data
ps = readRDS('ps_filtered.rds')


# Core Microbiome Analysis
  

# Convert to relative abundance, use tax_glom 
ps = transform(ps, 'compositional')
ps_genus = tax_glom(ps, 'Genus')

# Filter to hindgut samples only 
ps_genus_final = subset_samples(ps_genus, sample_type == "hindgut")

# Subset your phyloseq object by swim performance group (Fast, Moderate, Slow)
performance.fast = subset_samples(ps_genus_final, swim_performance %in% c('accelerator', 'crusier sprinter', 'manoeuvrer'))
performance.slow = subset_samples(ps_genus_final, swim_performance %in% c('flow refuging', 'burrowing'))
performance.moderate = subset_samples(ps_genus_final, swim_performance == 'generalist')

# Find core members, Detection Level = 0.001, Prevalence Level = 0.10
ASVs_fast = core_members(performance.fast, detection=0.001, prevalence = 0.1)
ASVs_slow = core_members(performance.slow, detection=0.001, prevalence = 0.1)
ASVs_moderate = core_members(performance.moderate, detection=0.001, prevalence = 0.1)

# Venn Diagrams Fast vs. Slow Swimming Performance
fast_v_slow <- ggVennDiagram(list(ASVs_fast, ASVs_slow),
              set_size = 6,
              category.names = c('Fast','Slow') )
ggsave("fast_v_slow.png", fast_v_slow)

# Fast vs. Moderate
fast_v_moderate <- ggVennDiagram(list(ASVs_fast, ASVs_moderate), 
                          set_side = 6, 
                          category.names = c('Fast', 'Moderate'))
ggsave("fast_v_moderate.png", fast_v_moderate)

# Moderate vs. Slow

moderate_v_slow <- ggVennDiagram(list(ASVs_moderate, ASVs_slow), 
              set_side = 6, 
              category.names = c('Moderate', 'Slow'))
ggsave("moderate_v_slow.png", moderate_v_slow)

# Fast vs. Moderate vs. Slow

fast_v_moderate_v_slow <- ggVennDiagram(list(ASVs_fast, ASVs_moderate, ASVs_slow), 
              set_side = 6, 
              category.names = c('Fast', 'Moderate', 'Slow'))
ggsave("fast_v_moderate_v_slow.png", fast_v_moderate_v_slow)


# View the ASVs
tax_table(prune_taxa(ASVs_fast, ps_genus))
tax_table(prune_taxa(ASVs_slow, ps_genus))
tax_table(prune_taxa(ASVs_moderate, ps_genus))

# Make a list
swim_performance_list <- list("Fast" = ASVs_fast, "Slow" = ASVs_slow, "Moderate" = ASVs_moderate)
