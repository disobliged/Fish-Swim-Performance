# This script calculates how much each taxon contributes to each pathway.
# The idea is to figure out if there are any taxa that contribute a large 
#   proportion of the reads to a given pathway, which would suggest that they 
#   are the main drivers of that pathway in the community.
# If there are no clear high-contributing taxa, that would suggest that the 
#   pathway is more evenly distributed across the community, and there are no 
#   clear 'drivers' of that pathway.

# Credit: Izaak Yip generated the stratified PICRUSt2 data used here.
#         Sam Donato drafted the first version of this script. 
#         Dr. Avril Metcalfe-Roach completed the second version of this script. 
#         Sam Donato visualized the taxonomic bar plot. (March 2026)

# Stratified PICRUSt2 reference: https://github-wiki-see.page/m/picrust/picrust2/wiki/PICRUSt2-Tutorial-(v2.5.0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load packages 
library(tidyverse)
library(readr)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(readxl)

# Define Pathways of interest
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

# PART 1 - Getting the data ready
# * Aggregate to the desired taxonomic level 
#   CURRENTLY FAMILY LEVEL
# * Convert to relative abundance to control for sequencing depth

# Load Fish Metadata 
metadata <- read_excel("fish_metadata_2.xlsx")

# Load stratified pathway contribution data (MetaCyc used)
strat_mc_data_1 <- data.table::fread("path_abun_contrib.tsv")

# Rename column name from #SampleID to sample 
metadata <- metadata %>%
  rename(sample = `#SampleID`)

# Left join sample names from metadata to stratified MetaCyc data
strat_mc_data_1 <- data.table::fread("path_abun_contrib.tsv") %>%
  left_join(metadata, by = "sample")

# Filter for hindgut samples
strat_hindgut <- strat_mc_data_1 %>%
  filter(sample_type == "hindgut")

# Filter for fast, moderate, slow groups
strat_fast <- strat_hindgut %>%
  filter(swim_performance %in% c("accelerator", "crusier sprinter", "manoeuvrer"))

strat_slow <- strat_hindgut %>%
  filter(swim_performance %in% c("flow refuging", "burrowing"))

strat_moderate <- strat_hindgut %>%
  filter(swim_performance == "generalist")

# Convert reads per sample to relative abundance to control for sequencing depth
strat_rel_fast = strat_fast %>% 
  group_by(sample) %>%
  mutate(rel_abun = taxon_function_abun/sum(taxon_function_abun)) %>%
  ungroup()

strat_rel_moderate = strat_moderate %>% 
  group_by(sample) %>%
  mutate(rel_abun = taxon_function_abun/sum(taxon_function_abun)) %>%
  ungroup()

strat_rel_slow = strat_slow %>% 
  group_by(sample) %>%
  mutate(rel_abun = taxon_function_abun/sum(taxon_function_abun)) %>%
  ungroup()


# Only keep the pathways of interest
strat_filt_fast = strat_rel_fast %>% 
  # Only include pathways of interest
  filter(`function` %in% energy_pathways)

strat_filt_moderate = strat_rel_moderate %>% 
  # Only include pathways of interest
  filter(`function` %in% energy_pathways)

strat_filt_slow = strat_rel_slow %>% 
  # Only include pathways of interest
  filter(`function` %in% energy_pathways)


# Aggregate reads to the Family level
# Load tax data, keep just the Family info
tax <- read_tsv("data_taxonomy.tsv") %>%  # Load Taxonomy Data Names
  separate(`Taxon`, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";") %>%
  select(taxon = `Feature ID`, Family) %>% 
  mutate(Family = ifelse(is.na(Family),taxon, Family)) # If no Family name, keep the feature ID for now

# Aggregate the ASVs to the Family level
strat_Family_fast = strat_filt_fast %>%
  left_join(tax) %>%
  group_by(Family, `function`, sample, swim_performance) %>% 
  # This is the 'amount' of the function contributed by each Family
  summarize(rel_abun = sum(rel_abun)) %>% 
  ungroup()

# Aggregate the ASVs to the Family level
strat_Family_moderate = strat_filt_moderate %>%
  left_join(tax) %>%
  group_by(Family, `function`, sample, swim_performance) %>% 
  # This is the 'amount' of the function contributed by each Family
  summarize(rel_abun = sum(rel_abun)) %>% 
  ungroup()

# Aggregate the ASVs to the Family level
strat_Family_slow = strat_filt_slow %>%
  left_join(tax) %>%
  group_by(Family, `function`, sample, swim_performance) %>% 
  # This is the 'amount' of the function contributed by each Family
  summarize(rel_abun = sum(rel_abun)) %>% 
  ungroup()


# PART 2 - Calculating the contributions of each taxon to each function

# Add together the counts from each sample for each taxon and function
# combine all the relative abundances for a given taxon and pathway across all samples to get an average contribution of that taxon to that pathway
strat_avg_fast = strat_Family_fast %>% 
  group_by(Family, `function`) %>% 
  summarize(total_contrib_per_Family = sum(rel_abun)) %>% 
  ungroup()

strat_avg_moderate = strat_Family_moderate %>% 
  group_by(Family, `function`) %>% 
  summarize(total_contrib_per_Family = sum(rel_abun)) %>% 
  ungroup()

strat_avg_slow = strat_Family_slow %>% 
  group_by(Family, `function`) %>% 
  summarize(total_contrib_per_Family = sum(rel_abun)) %>% 
  ungroup()


# Now we just have three columns: the Family, the function, and the function's abundance that belongs to that taxon.
# Next, we figure out the proportion of each function that belongs to each taxon.
strat_prop_fast = strat_avg_fast %>% 
  group_by(`function`) %>%
  mutate(prop_tax_reads_per_function = total_contrib_per_Family/sum(total_contrib_per_Family)) %>% 
  ungroup() %>% 
  select(-total_contrib_per_Family)

strat_prop_moderate = strat_avg_moderate %>% 
  group_by(`function`) %>%
  mutate(prop_tax_reads_per_function = total_contrib_per_Family/sum(total_contrib_per_Family)) %>% 
  ungroup() %>% 
  select(-total_contrib_per_Family) 

strat_prop_slow = strat_avg_slow %>% 
  group_by(`function`) %>%
  mutate(prop_tax_reads_per_function = total_contrib_per_Family/sum(total_contrib_per_Family)) %>% 
  ungroup() %>% 
  select(-total_contrib_per_Family) # We don't need the total contribution per Family anymore, just the proportion

# Sanity check. At this point, the abundances for each function should add up to 1, 
# and the individual values represent proportions.
strat_prop_fast %>% group_by(`function`) %>% 
  summarize(sum = sum(prop_tax_reads_per_function)) %>% 
  ungroup() # Every pathway adds up to 1! Good!

strat_prop_moderate %>% group_by(`function`) %>% 
  summarize(sum = sum(prop_tax_reads_per_function)) %>% 
  ungroup()

strat_prop_slow %>% group_by(`function`) %>% 
  summarize(sum = sum(prop_tax_reads_per_function)) %>% 
  ungroup()

# PART 3 - get a feel for the data

# Filter for taxa that contribute at least 33% of the reads to a given pathway
strat_cutoff_fast = strat_prop_fast %>% 
  filter(prop_tax_reads_per_function >= 0.33) # 1

strat_cutoff_moderate = strat_prop_moderate %>% 
  filter(prop_tax_reads_per_function >= 0.33) # None!

strat_cutoff_slow = strat_prop_slow %>% 
  filter(prop_tax_reads_per_function >= 0.33) # None!

# Top hits:
head(strat_prop_fast %>% arrange(-prop_tax_reads_per_function))
head(strat_prop_moderate %>% arrange(-prop_tax_reads_per_function))
head(strat_prop_slow %>% arrange(-prop_tax_reads_per_function))

# Histogram of the proportions, just to see the distributions. 
# No clear high-contributing taxa - a taxonomic bar plot-style figure might work best.
hist(log10(strat_prop_fast$prop_tax_reads_per_function),breaks=50)
hist(log10(strat_prop_moderate$prop_tax_reads_per_function),breaks=50)
hist(log10(strat_prop_slow$prop_tax_reads_per_function),breaks=50)

# Part 4 - Visualizing taxonomic bar plot 

# Create the 5% cutoff
threshold <- 0.05  

# Group the <5% reads by function and family, sum all the rare taxa 
strat_prop_grouped_fast <- strat_prop_fast %>%
  group_by(`function`) %>%
  mutate(Family = ifelse(prop_tax_reads_per_function < threshold,
                         "<5% of reads",
                         Family)) %>%
  ungroup() %>%
  group_by(`function`, Family) %>%
  summarize(prop_tax_reads_per_function = sum(prop_tax_reads_per_function),
            .groups = "drop")

strat_prop_grouped_moderate <- strat_prop_moderate %>%
  group_by(`function`) %>%
  mutate(Family = ifelse(prop_tax_reads_per_function < threshold,
                         "<5% of reads",
                         Family)) %>%
  ungroup() %>%
  group_by(`function`, Family) %>%
  summarize(prop_tax_reads_per_function = sum(prop_tax_reads_per_function),
            .groups = "drop")

strat_prop_grouped_slow <- strat_prop_slow %>%
  group_by(`function`) %>%
  mutate(Family = ifelse(prop_tax_reads_per_function < threshold,
                         "<5% of reads",
                         Family)) %>%
  ungroup() %>%
  group_by(`function`, Family) %>%
  summarize(prop_tax_reads_per_function = sum(prop_tax_reads_per_function),
            .groups = "drop")

# Plot without a legend for visibility
plot_avg_fast <- strat_prop_grouped_fast %>%
  mutate(`function` = str_to_title(`function`)) %>%
  ggplot(aes(`function`, prop_tax_reads_per_function, fill = Family)) +
  labs(y = "Proportion of Family taxa reads per function") +
  geom_col(position = "stack", color = NA) +
  theme_classic(base_size = 22) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
       legend.position = "none")+
  xlab(NULL)

plot_avg_fast

plot_avg_moderate <- strat_prop_grouped_moderate %>%
  mutate(`function` = str_to_title(`function`)) %>%
  ggplot(aes(`function`, prop_tax_reads_per_function, fill = Family)) +
  labs(y = "Proportion of Family taxa reads per function") +
  geom_col(position = "stack", color = NA) +
  theme_classic(base_size = 22) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")+
  xlab(NULL)

plot_avg_moderate 

plot_avg_slow <- strat_prop_grouped_slow %>%
  mutate(`function` = str_to_title(`function`)) %>%
  ggplot(aes(`function`, prop_tax_reads_per_function, fill = Family)) +
  labs(y = "Proportion of Family taxa reads per function") +
  geom_col(position = "stack", color = NA) +
  theme_classic(base_size = 22) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")+
  xlab(NULL)

plot_avg_slow

# Plot with a legend, as well as convert an unclassified Family taxa 

strat_prop_grouped_fast <- strat_prop_grouped_fast %>%
  mutate(Family = ifelse(Family == "3a4580bd0eb64ee4701c414a8ba63ee2",
                         "Unclassified",
                         Family))

plot_avg_fast <- strat_prop_grouped_fast %>%
  mutate(`function` = str_to_title(`function`)) %>%
  ggplot(aes(`function`, prop_tax_reads_per_function, fill = Family)) +
  labs(y = "Fast Fish Proportion of Family taxa reads per function") +
  geom_col(position = "stack", color = NA) +
  theme_classic(base_size = 22) +
  theme(strip.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 16)) +
  xlab(NULL)

plot_avg_fast
ggsave("plot_avg_fast.jpg", plot_avg)

plot_avg_moderate <- strat_prop_grouped_moderate %>%
  mutate(`function` = str_to_title(`function`)) %>%
  ggplot(aes(`function`, prop_tax_reads_per_function, fill = Family)) +
  labs(y = "Moderate Fish Proportion of Family taxa reads per function") +
  geom_col(position = "stack", color = NA) +
  theme_classic(base_size = 22) +
  theme(strip.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 16)) +
  xlab(NULL)

plot_avg_moderate
  
plot_avg_slow <- strat_prop_grouped_slow %>%
  mutate(`function` = str_to_title(`function`)) %>%
  ggplot(aes(`function`, prop_tax_reads_per_function, fill = Family)) +
  labs(y = "Slow Fish Proportion of Family taxa reads per function") +
  geom_col(position = "stack", color = NA) +
  theme_classic(base_size = 22) +
  theme(strip.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 16)) +
  xlab(NULL)

plot_avg_slow
