# Consolidated Bash Scripts
Author: IY <br>
Created: Apr 5, 2026 <br>
Last Modified: Apr 5, 2026

This document contains the Qiime2 and PICRUSt2 commands we have used to build our project. Common Bash functions will not be included here. Line breaks are included in our code for clarity, but are not neccessary for the function to run.

## Qiime2

Qiime was mainly used to generate our initial analyses before we did more complex analyses in R. All ```.qzv``` files were viewed with

### Importing the Fish Dataset with a Manifest

```
# A directory for our Qiime analyses was made prior to importing.
# Using a manifest helped us simultaneously generate a demultiplexed file.
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path ../fish/fish_manifest.tsv \
  --output-path demuxFishSeqs.qza
```

### Demultiplexed Visualization

```
qiime demux summarize
  --i-data demuxFishSeqs.qza \
  --o-visualization demuxFish.qzv
```

### Denoising

```
# Denoise
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demuxFishSeqs.qza \
  --p-trim-left 0 \
  --p-trunc-len 228 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats fishDenoiseStats.qza

# Visualize the denoised data
qiime metadata tabulate \
  --m-input-file fishDenoiseStats.qza \
  --o-visualization fishDenoiseStats.qzv
```

### Clustering

```
# Visualize ASVs stats
qiime feature-table summarize \
  --i-table fishTable.qza \
  --o-visualization fishTable.qzv \
  --m-sample-metadata-file ../fish/fish_metadata.txt
  
qiime feature-table tabulate-seqs \
  --i-data repFishSeqs.qza \
  --o-visualization repFishSeqs.qzv
```

### Taxonomic Analysis

```
# We used the silva-138-99-nb-classifier.qza classifier for our dataset
qiime feature-classifier extract-reads \
  --i-sequences /datasets/silva_ref_files/silva-138-99-seqs.qza \
  --p-f-primer GTGYCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACNVGGGTWTCTAAT \
  --p-trunc-len 228 \
  --o-reads fishRefSeqsTrimmed.qza
```

### Taxonomy-Based Filtering

```
qiime taxa filter-table \
  --i-table fishTable.qza \
  --i-taxonomy fishTaxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table fishTableFiltered.qza
```

### Visualize the Filtered Table

```
qiime feature-table summarize \
  --i-table fishTableFiltered.qza \
  --o-visualization fishTableFiltered.qzv \
  --m-sample-metadata-file ../fish/fish_metadata.txt
```

### Tree Generation

```
# Generate a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences repFishSeqs.qza \
  --o-alignment alignedRepFishSeqs.qza \
  --o-masked-alignment maskedAlignedRepFishSeqs.qza \
  --o-tree unrootedFishTree.qza \
  --o-rooted-tree rootedFishTree.qza
```

### Alpha Rarefaction
```
# Alpha-rarefaction
qiime diversity alpha-rarefaction \
  --i-table fishTable.qza \
  --i-phylogeny rootedFishTree.qza \
  --p-max-depth 2525 \
  --m-metadata-file ../fish/fish_metadata.txt \
  --o-visualization fishAlphaRarefaction.qzv
```

### Calculate Alpha Diversity Metrics
```
# Calculate the core metrics
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rootedFishTree.qza \
  --i-table fishTableFiltered.qza \
  --p-sampling-depth 2525 \
  --m-metadata-file ../fish/fish_metadata.txt \
  --output-dir fishCoreMetricsResults

# Calculate alpha-group-significance

qiime diversity alpha-group-significance \
  --i-alpha-diversity fishCoreMetricsResults/shannon_vector.qza \
  --m-metadata-file ../fish/fish_metadata.txt \
  --o-visualization fishCoreMetricsResults/shannonFish.qzv
```

## PICRUSt2

We used PICRUSt2 to provide us with pathway analysis. We used both unstratified and stratified PICRUSt2 outputs.

### Data preparation
```
# Filter the fishTable.qza file to remove features with 5 or lower counts for efficiency.
qiime feature-table filter-features \
  --i-table fishTableFiltered.qza \
  --p-min-frequency 5 \
  --o-filtered-table fishFreqFiltered.qza

# Make the PICRUSt2 directory and export the data
mkdir picrust

qiime tools export \
   --input-path fishFreqFiltered.qza \
   --output-path picrust

qiime tools export \
   --input-path repFishSeqs.qza \
   --output-path picrust

# Switch from the Qiime2 to the PICRUSt2 environment

conda deactivate
conda activate picrust2
```

### Unstratified PICRUSt2 Analysis
```
picrust2_pipeline.py \
   -s picrust/dna-sequences.fasta \
   -i picrust/feature-table.biom \
   -o picrust_output
```

### Stratified PICRUSt2 Analysis
```
picrust2_pipeline.py \
-s rep_seqs.fasta \
-i table.biom \
-o picrust_strat_output \
--stratified # The stratified flag is appended to the command for the stratified analysis to run
```

## References

We are thankful for the bioinformatics community to provide these tools to help build our project and have cited their tools below.

### Qiime2
>Bolyen E, Rideout JR, Dillon MR, Bokulich NA, Abnet CC, Al-Ghalith GA, Alexander H, Alm EJ, Arumugam M, Asnicar F, Bai Y, Bisanz JE, Bittinger K, Brejnrod A, Brislawn CJ, Brown CT, Callahan BJ, Caraballo-Rodríguez AM, Chase J, Cope EK, Da Silva R, Diener C, Dorrestein PC, Douglas GM, Durall DM, Duvallet C, Edwardson CF, Ernst M, Estaki M, Fouquier J, Gauglitz JM, Gibbons SM, Gibson DL, Gonzalez A, Gorlick K, Guo J, Hillmann B, Holmes S, Holste H, Huttenhower C, Huttley GA, Janssen S, Jarmusch AK, Jiang L, Kaehler BD, Kang KB, Keefe CR, Keim P, Kelley ST, Knights D, Koester I, Kosciolek T, Kreps J, Langille MGI, Lee J, Ley R, Liu YX, Loftfield E, Lozupone C, Maher M, Marotz C, Martin BD, McDonald D, McIver LJ, Melnik AV, Metcalf JL, Morgan SC, Morton JT, Naimey AT, Navas-Molina JA, Nothias LF, Orchanian SB, Pearson T, Peoples SL, Petras D, Preuss ML, Pruesse E, Rasmussen LB, Rivers A, Robeson MS, Rosenthal P, Segata N, Shaffer M, Shiffer A, Sinha R, Song SJ, Spear JR, Swafford AD, Thompson LR, Torres PJ, Trinh P, Tripathi A, Turnbaugh PJ, Ul-Hasan S, van der Hooft JJJ, Vargas F, Vázquez-Baeza Y, Vogtmann E, von Hippel M, Walters W, Wan Y, Wang M, Warren J, Weber KC, Williamson CHD, Willis AD, Xu ZZ, Zaneveld JR, Zhang Y, Zhu Q, Knight R, and Caporaso JG. 2019. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology 37: 852–857. https://doi.org/10.1038/s41587-019-0209-9

### Qiime2View
>Bolyen E, Rideout JR, Dillon MR, Bokulich NA, Abnet CC, Al-Ghalith GA, Alexander H, Alm EJ, Arumugam M, Asnicar F, Bai Y, Bisanz JE, Bittinger K, Brejnrod A, Brislawn CJ, Brown CT, Callahan BJ, Caraballo-Rodríguez AM, Chase J, Cope EK, Da Silva R, Diener C, Dorrestein PC, Douglas GM, Durall DM, Duvallet C, Edwardson CF, Ernst M, Estaki M, Fouquier J, Gauglitz JM, Gibbons SM, Gibson DL, Gonzalez A, Gorlick K, Guo J, Hillmann B, Holmes S, Holste H, Huttenhower C, Huttley GA, Janssen S, Jarmusch AK, Jiang L, Kaehler BD, Kang KB, Keefe CR, Keim P, Kelley ST, Knights D, Koester I, Kosciolek T, Kreps J, Langille MGI, Lee J, Ley R, Liu YX, Loftfield E, Lozupone C, Maher M, Marotz C, Martin BD, McDonald D, McIver LJ, Melnik AV, Metcalf JL, Morgan SC, Morton JT, Naimey AT, Navas-Molina JA, Nothias LF, Orchanian SB, Pearson T, Peoples SL, Petras D, Preuss ML, Pruesse E, Rasmussen LB, Rivers A, Robeson MS, Rosenthal P, Segata N, Shaffer M, Shiffer A, Sinha R, Song SJ, Spear JR, Swafford AD, Thompson LR, Torres PJ, Trinh P, Tripathi A, Turnbaugh PJ, Ul-Hasan S, van der Hooft JJJ, Vargas F, Vázquez-Baeza Y, Vogtmann E, von Hippel M, Walters W, Wan Y, Wang M, Warren J, Weber KC, Williamson CHD, Willis AD, Xu ZZ, Zaneveld JR, Zhang Y, Zhu Q, Knight R, and Caporaso JG. 2019. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology 37: 852–857. https://doi.org/10.1038/s41587-019-0209-9
Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP. 2016. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13:581–583.
McDonald D, Clemente JC, Kuczynski J, Rideout JR, Stombaugh J, Wendel D, Wilke A, Huse S, Hufnagle J, Meyer F, Knight R, Caporaso JG. 2012. The Biological Observation Matrix (BIOM) format or: how I learned to stop worrying and love the ome-ome. GigaScience 1:7.
Price MN, Dehal PS, Arkin AP. 2010. FastTree 2–approximately maximum-likelihood trees for large alignments. PloS one 5:e9490.
Lane D. 1991. 16S/23S rRNA sequencing, p. 115–175. In Stackebrandt, E, Goodfellow, M (eds.), Nucleic Acid Techniques in Bacterial Systematics. John Wiley, New York.
Katoh K, Standley DM. 2013. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular biology and evolution 30:772–780.

### PICRUSt2
>Douglas, G.M., Maffei, V.J., Zaneveld, J.R. et al. PICRUSt2 for prediction of metagenome functions. Nat Biotechnol 38, 685–688 (2020). https://doi.org/10.1038/s41587-020-0548-6
Pierre Barbera, Alexey M Kozlov, Lucas Czech, Benoit Morel, Diego Darriba, Tomáš Flouri, Alexandros Stamatakis; EPA-ng: Massively Parallel Evolutionary Placement of Genetic Sequences, Systematic Biology, syy054, https://doi.org/10.1093/sysbio/syy054
Louca S. 2017. Castor: Efficient phylogenetics on large trees. CRAN: Contributed Packages.
Ye Y, Doak TG. 2009. A parsimony approach to biological pathway reconstruction/inference for genomes and Metagenomes. PLoS Computational Biology 5. 
