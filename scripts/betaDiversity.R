#loading libraries
library(tidyverse)
library(phyloseq)
library(vegan)

#loading ps
ps = readRDS('../datasets/ps_filtered.rds') # TODO check file path

# Subset hindgut samples 
ps_hindgut = subset_samples(ps, sample_type == "hindgut") # 97 hindgut samples 

### Remove the two obvious outliers ###
# 1/ 13414.100.hg.R1.fastq.gz has MDS1 of -8542.9979 when included in bray-curtis 
# 1/ 13414.100.hg.R1.fastq.gz has MDS1 of 0.97450528 when included in jaccard
# 2/ 13414.34.hg.R1.fastq.gz has MDS2 of -3.4208732669 when included in bray-curtis
# 2/ 13414.34.hg.R1.fastq.gz has MDS2 of 1.724540e-04 when included in jaccard
to_remove <- c("13414.100.hg.R1.fastq.gz", "13414.34.hg.R1.fastq.gz")
# ps <- subset_samples(ps, !(SampleID %in% to_remove)) # can't use this bcs SampleID not a column name
ps_rm_outlier <- prune_samples(!(sample_names(ps_hindgut) %in% to_remove), ps_hindgut) # 95 samples left after removing 2 


# Rarefying data 
# (sample depth selected from alpha diversity)
set.seed(1)
# For ps with outliers removed
psrare_rm_outlier = ps_rm_outlier %>% rarefy_even_depth(sample.size = 2525, rngseed = 1)
# For ps with all samples
# psrare = ps_hindgut %>% rarefy_even_depth(sample.size = 2525, rngseed = 1)


## Calculating metrics ##
# 1/ Bray-Curtis
ps_bray = phyloseq::distance(psrare_rm_outlier, method = "bray")
# 2/ (Jaccard) 
ps_jaccard <- phyloseq::distance(psrare_rm_outlier, method = "jaccard", binary = TRUE)
# View(as.matrix(ps_bray)) # for human interpretation
# View(as.matrix(ps_jaccard)) # for human interpretation


## MDS scaling ##
set.seed(1)
mds_bray = metaMDS(ps_bray) 
mds_jaccard = metaMDS(ps_jaccard) 

# Extracting data
mds_data_bray = mds_bray$points %>% as.data.frame %>%  
  merge(sample_data(psrare_rm_outlier), by='row.names', sort=F) %>%
  # Sorting data by swim performance group (Fast, Moderate, Slow)
  mutate(speed_category = ifelse(swim_performance %in% c('flow refuging', 'burrowing'), 'slow',
                                 ifelse(swim_performance %in% c('generalist'), 'moderate', 'fast')), .before = swim_performance)

mds_data_jaccard = mds_jaccard$points %>% as.data.frame %>%  
  merge(sample_data(psrare_rm_outlier), by='row.names', sort=F) %>%
  mutate(speed_category = ifelse(swim_performance %in% c('flow refuging', 'burrowing'), 'slow',
                                 ifelse(swim_performance %in% c('generalist'), 'moderate', 'fast')), .before = swim_performance)


## Plotting ##

# 1/ Bray-Curtis 
p_bray = mds_data_bray %>%
  ggplot(aes(MDS1,MDS2,color = speed_category)) + 
  geom_point() + #scatterplot 
  stat_ellipse() + # 95% confidence interval around median
  labs(color = "Swim Performance") +
  theme_classic(base_size=18) +
  scale_color_manual(values = c("fast" = "#8E7CA6", 
                                "moderate" = "#E07A9A", 
                                "slow" = "#fcba65")) #colour palette
p_bray

# 2/ Jaccard 
p_jaccard = mds_data_jaccard %>%
  ggplot(aes(MDS1,MDS2,color = speed_category)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(color = "Swim Performance") +
  theme_classic(base_size=18) +
  scale_color_manual(values = c("fast" = "#8E7CA6", 
                                "moderate" = "#E07A9A", 
                                "slow" = "#fcba65")) #colour palette
p_jaccard


### PERMANOVA Tests ###
metadata = sample_data(mds_data_bray) %>% data.frame() #%>% na.omit() #or use colSums(is.na(metadata)) 

# Running adonis on data to get Single variable comparison
stats_univar_bray = adonis2(ps_bray ~ speed_category, data = metadata, na.action = na.omit)
# p-value = 0.005 significant
# R2 = 0.07235082 our model is not the best at explaining variance in our diversity
stats_univar_jaccard = adonis2(ps_jaccard ~ speed_category, data = metadata, na.action = na.omit)
# p-value = 0.05 borderline significant
# R2 = 0.05831928 our model is not the best at explaining variance in our diversity


## Comparing Two Speeds ##
# 1/ Compare fast and slow
# 2/ Compare fast and moderate
# 3/ Compare moderate and slow
to_compare_fast = c('accelerator', 'crusier sprinter', 'manoeuvrer')
to_compare_mod = c('generalist')
to_compare_slow = c('flow refuging', 'burrowing')

# 1/ fast vs slow
samples_to_keep_fast_slow = metadata %>% 
  filter(swim_performance %in% c(to_compare_fast, to_compare_slow)) %>% 
  select(Row.names) %>% 
  pull() #to turn output to char vector
# Subsetting distance matrix to only include fast & slow
ps_bray_sub_fast_slow = as.matrix(ps_bray)[samples_to_keep_fast_slow,samples_to_keep_fast_slow] %>% as.dist()
ps_jacc_sub_fast_slow = as.matrix(ps_jaccard)[samples_to_keep_fast_slow,samples_to_keep_fast_slow] %>% as.dist()
# Running PERMANOVA on fast and slow
stats_bray_fast_slow = adonis2(ps_bray_sub_fast_slow ~ swim_performance, #metadata has speed_category, ps doesn't
                               data = metadata %>% filter(swim_performance %in% c(to_compare_fast, to_compare_slow)))
stats_jacc_fast_slow = adonis2(ps_jacc_sub_fast_slow ~ swim_performance, 
                               data = metadata %>% filter(swim_performance %in% c(to_compare_fast, to_compare_slow)))
# View(stats_bray_fast_slow)
# View(stats_jacc_fast_slow)
# p=0.728 fast vs slow are not significant in both Bray and Jaccard


# 2/ fast vs moderate
samples_to_keep_fast_mod = metadata %>% 
  filter(swim_performance %in% c(to_compare_fast, to_compare_mod)) %>% 
  select(Row.names) %>% 
  pull() 
# Subsetting distance matrix to only include fast & moderate
ps_bray_sub_fast_mod = as.matrix(ps_bray)[samples_to_keep_fast_mod,samples_to_keep_fast_mod] %>% as.dist()
ps_jacc_sub_fast_mod = as.matrix(ps_jaccard)[samples_to_keep_fast_mod,samples_to_keep_fast_mod] %>% as.dist()
# Running PERMANOVA on fast and moderate
stats_bray_fast_mod = adonis2(ps_bray_sub_fast_mod ~ swim_performance,
                              data = metadata %>% filter(swim_performance %in% c(to_compare_fast, to_compare_mod)))
stats_jacc_fast_mod = adonis2(ps_jacc_sub_fast_mod ~ swim_performance, 
                              data = metadata %>% filter(swim_performance %in% c(to_compare_fast, to_compare_mod)))
View(stats_bray_fast_mod)
View(stats_jacc_fast_mod)
# fast vs moderate are not significant in both Bray and Jaccard


  
# 3/ moderate & slow
samples_to_keep_mod_slow = metadata %>% 
  filter(swim_performance %in% c(to_compare_mod, to_compare_slow)) %>% 
  select(Row.names) %>% 
  pull() 
# Subsetting distance matrix to only include moderate & slow
ps_bray_sub_mod_slow = as.matrix(ps_bray)[samples_to_keep_mod_slow,samples_to_keep_mod_slow] %>% as.dist()
ps_jacc_sub_mod_slow = as.matrix(ps_jaccard)[samples_to_keep_mod_slow,samples_to_keep_mod_slow] %>% as.dist()
# Running PERMANOVA on moderate and slow
stats_bray_mod_slow = adonis2(ps_bray_sub_mod_slow ~ swim_performance, 
                              data = metadata %>% filter(swim_performance %in% c(to_compare_mod, to_compare_slow)))
stats_jacc_mod_slow = adonis2(ps_jacc_sub_mod_slow ~ swim_performance, 
                              data = metadata %>% filter(swim_performance %in% c(to_compare_mod, to_compare_slow)))
View(stats_bray_mod_slow)
View(stats_jacc_mod_slow)
# moderate and slow are not significant in both Bray and Jaccard



# Combining stats into a table
stats = bind_rows('Univariate' = stats_univar_bray %>% as.data.frame() %>% 
                    rownames_to_column('Variable') %>%
                    cbind(Metric = c('Bray')),
                  'Fast v. Slow' = stats_bray_fast_slow %>% as.data.frame() %>% 
                    rownames_to_column('Variable') %>%
                    cbind(Metric = c('Bray')),
                  'Fast v. Moderate' = stats_bray_fast_mod %>% as.data.frame() %>% 
                    rownames_to_column('Variable') %>%
                    cbind(Metric = c('Bray')),
                  'Moderate v. Slow' = stats_bray_mod_slow %>% as.data.frame() %>% 
                    rownames_to_column('Variable') %>%
                    cbind(Metric = c('Bray')),
                  
                  'Univariate' = stats_univar_jaccard %>% as.data.frame() %>% 
                    rownames_to_column('Variable') %>%
                    cbind(Metric = c('Jaccard')),
                  'Fast v. Slow' = stats_jacc_fast_slow %>% as.data.frame() %>% 
                    rownames_to_column('Variable') %>%
                    cbind(Metric = c('Jaccard')),
                  'Fast v. Moderate' = stats_jacc_fast_mod %>% as.data.frame() %>% 
                    rownames_to_column('Variable') %>%
                    cbind(Metric = c('Jaccard')),
                  'Moderate v. Slow' = stats_jacc_mod_slow %>% as.data.frame() %>% 
                    rownames_to_column('Variable') %>%
                    cbind(Metric = c('Jaccard')),
                  .id = 'Comparison') # Naming the first column 'Comparison' 

stats_clean = stats %>% rename(Pval = `Pr(>F)`) %>% #renaming column to Pval
  filter(!is.na(Pval)) %>% # Removing irrelevant lines 
  select(-Variable) %>%
  relocate(Metric, .after = Comparison)
#View(stats_clean) 

### Saving Output
# Save plot
ggsave('../plots/beta_diversity_plot_bray.jpeg',
       plot = p_bray,
       height= 4, width = 6)

ggsave('../plots/beta_diversity_plot_jaccard.jpeg',
       plot = p_jaccard,
       height= 4, width = 6)

# Save statistics as an Excel file
writexl::write_xlsx(list('Beta Diversity Stats' = stats_clean),
                    '../datasets/beta_diversity_stats.xlsx')

