## Construct feature table ##

# Import
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path paired-end_sequences/lane_1 \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path paired-end_sequences/demux-paired-end.qza
#Sequences from the RTSF are in CASAVA-1.8 format.

qiime demux summarize \
  --i-data paired-end_sequences/demux-paired-end.qza \
  --o-visualization paired-end_sequences/demux-paired-end.qzv


# Quality control and feature table construction with DADA2
mkdir feature_table

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end_sequences/demux-paired-end.qza \
  --p-trunc-len-f 235\
  --p-trunc-len-r 235\
  --p-trim-left-f 19 \
  --p-trim-left-r 18 \
  --p-n-threads 40 \
  --verbose \
  --o-denoising-stats feature_table/dada2-stats.qza \
  --o-representative-sequences feature_table/rep-seqs.qza \
  --o-table feature_table/table.qza
#Reads are truncated to 235 to remove lower quality positions towards the ends. Additionally, the forward and reverse primer are trimmed to prevent the ambiguous nucleotides from causing a large number of false-positive chimeras to be identified.
#This step also removes chimeras.

qiime metadata tabulate \
  --m-input-file feature_table/dada2-stats.qza \
  --o-visualization feature_table/dada2-stats.qzv


# Visualize feature table data summaries
qiime feature-table summarize \
  --i-table feature_table/table.qza \
  --m-sample-metadata-file sample-metadata.tsv \
  --o-visualization feature_table/table.qzv

qiime feature-table tabulate-seqs \
  --i-data feature_table/rep-seqs.qza \
  --o-visualization feature_table/rep-seqs.qzv
  
  
  ## Taxonomic analysis ##
#For more info see:  https://docs.qiime2.org/2018.11/tutorials/moving-pictures/
mkdir taxonomy

qiime feature-classifier classify-sklearn \
  --i-classifier feature_classifier/gg13_8_99_799F-1193R.qza \
  --i-reads feature_table/rep-seqs.qza \
  --p-n-jobs 40 \
  --o-classification taxonomy/taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy/taxonomy.qza \
  --o-visualization taxonomy/taxonomy.qzv


# Generate taxa bar plots
qiime taxa barplot \
  --i-table feature_table/table.qza \
  --i-taxonomy taxonomy/taxonomy.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization taxonomy/taxa-bar-plots.qzv


## Filtering ##

# Filter  mitochondria and chloroplast from feature table
qiime taxa filter-table \
  --i-table feature_table/table.qza \
  --i-taxonomy taxonomy/taxonomy.qza \
  --p-exclude mitochondria,Chloroplast \
  --o-filtered-table feature_table/table_noM_noC.qza


# Remove singletons
qiime feature-table filter-features \
  --i-table feature_table/table_noM_noC.qza \
  --p-min-samples 2 \
  --o-filtered-table feature_table/filtered_table.qza


# Remove rep sequences corresponding to mitochondria and chloroplast
qiime taxa filter-seqs \
  --i-sequences feature_table/rep-seqs.qza \
  --i-taxonomy taxonomy/taxonomy.qza \
  --p-exclude mitochondria,Chloroplast \
  --o-filtered-sequences feature_table/filtered_rep-seqs.qza


# Visualize data summaries
qiime feature-table tabulate-seqs \
  --i-data feature_table/filtered_rep-seqs.qza \
  --o-visualization feature_table/filtered_rep-seqs.qzv

qiime feature-table summarize \
  --i-table feature_table/filtered_table.qza \
  --m-sample-metadata-file sample-metadata.tsv \
  --o-visualization feature_table/filtered_table.qzv


# Generate taxa bar plots of filtered data
qiime taxa barplot \
  --i-table feature_table/filtered_table.qza \
  --i-taxonomy taxonomy/taxonomy.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization taxonomy/taxa-bar-plots_filtered.qzv


## Table parsing ##

# Create table for holoxenic/SynCom/input analyses
qiime feature-table filter-samples \
  --i-table feature_table/filtered_table.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-where "treatment IN ('SynCom-fresh', 'MSU19', 'input-SynCom', 'input-MSU19')" \
  --o-filtered-table feature_table/filtered_table_HO-SC.qza

qiime feature-table summarize \
  --i-table feature_table/filtered_table_HO-SC.qza \
  --m-sample-metadata-file sample-metadata.tsv \
  --o-visualization feature_table/filtered_table_HO-SC.qzv

# Create table for SynCom cryoprotection analyses
qiime feature-table filter-samples \
  --i-table feature_table/filtered_table.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-where "treatment IN ('SynCom-fresh', 'SynCom-glycerol-ice', 'SynCom-glycerol-37C', 'SynCom-DMSO-ice', 'SynCom-DMSO-37C')" \
  --o-filtered-table feature_table/filtered_table_CP.qza

qiime feature-table summarize \
  --i-table feature_table/filtered_table_CP.qza \
  --m-sample-metadata-file sample-metadata.tsv \
  --o-visualization feature_table/filtered_table_CP.qzv


## Generate a tree for phylogenetic diversity analyses ##
qiime alignment mafft \
  --i-sequences feature_table/filtered_rep-seqs.qza \
  --o-alignment feature_table/filtered_rep-seqs_aligned.qza

qiime alignment mask \
  --i-alignment feature_table/filtered_rep-seqs_aligned.qza \
  --o-masked-alignment feature_table/filtered_rep-seqs_aligned_masked.qza

qiime phylogeny fasttree \
  --i-alignment feature_table/filtered_rep-seqs_aligned_masked.qza \
  --o-tree taxonomy/unrooted-tree.qza

qiime phylogeny midpoint-root \
  --i-tree taxonomy/unrooted-tree.qza \
  --o-rooted-tree taxonomy/rooted-tree.qza


## Holoxenic/SynCom diversity analysis ##

#Alpha and beta diversity calculations
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny taxonomy/rooted-tree.qza \
  --i-table feature_table/filtered_table_HO-SC.qza \
  --p-sampling-depth 1705 \
  --m-metadata-file sample-metadata.tsv \
  --output-dir diversity_HO-SC
#p-sampling-depth using "Interactive Sample Detail" tab in filtered_table.qzv.


# Simpson diversity index calculations
qiime diversity alpha \
  --i-table feature_table/filtered_table_HO-SC.qza \
  --p-metric simpson \
  --o-alpha-diversity diversity_HO-SC/simpson_vector.qza


# Visualize alpha diversity significance testing
qiime diversity alpha-group-significance \
  --i-alpha-diversity diversity_HO-SC/simpson_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization diversity_HO-SC/simpson-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity diversity_HO-SC/observed_otus_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization diversity_HO-SC/observed_otus-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity diversity_HO-SC/faith_pd_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization diversity_HO-SC/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity diversity_HO-SC/evenness_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization diversity_HO-SC/evenness-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity diversity_HO-SC/shannon_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization diversity_HO-SC/shannon-group-significance.qzv


# Alpha rarefaction plotting
qiime diversity alpha-rarefaction \
  --i-table feature_table/filtered_table_HO-SC.qza \
  --i-phylogeny taxonomy/rooted-tree.qza \
  --p-max-depth 15000 \
  --p-steps 50 \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization diversity_HO-SC/alpha-rarefaction-plot_15000.qzv
#The value for --p-max-depth was determined by reviewing the “Frequency per sample” information presented in the filtered_table.qza file that was created above. 

qiime diversity alpha-rarefaction \
  --i-table feature_table/filtered_table_HO-SC.qza \
  --i-phylogeny taxonomy/rooted-tree.qza \
  --p-max-depth 5000 \
  --p-steps 50 \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization diversity_HO-SC/alpha-rarefaction-plot_5000.qzv
#The value for --p-max-depth was determined by reviewing the “Frequency per sample” information presented in the filtered_table.qza file that was created above. 


# Visualize beta diversity significance testing
qiime diversity beta-group-significance \
  --i-distance-matrix diversity_HO-SC/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column treatment \
  --o-visualization diversity_HO-SC/unweighted-unifrac-treatment-group-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix diversity_HO-SC/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column treatment \
  --o-visualization diversity_HO-SC/weighted-unifrac-treatment-group-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix diversity_HO-SC/jaccard_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column treatment \
  --o-visualization diversity_HO-SC/jaccard-treatment-group-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix diversity_HO-SC/bray_curtis_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column treatment \
  --o-visualization diversity_HO-SC/bray_curtis-treatment-group-significance.qzv \
  --p-pairwise


## SynCom cryoprotection diversity analysis ##

#Alpha and beta diversity calculations
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny taxonomy/rooted-tree.qza \
  --i-table feature_table/filtered_table_CP.qza \
  --p-sampling-depth 3170 \
  --m-metadata-file sample-metadata.tsv \
  --output-dir diversity_CP
#p-sampling-depth using "Interactive Sample Detail" tab in filtered_table.qzv.


# Simpson diversity index calculations
qiime diversity alpha \
  --i-table feature_table/filtered_table_CP.qza \
  --p-metric simpson \
  --o-alpha-diversity diversity_CP/simpson_vector.qza


# Alpha diversity significance testing
qiime diversity alpha-group-significance \
  --i-alpha-diversity diversity_CP/simpson_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization diversity_CP/simpson-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity diversity_CP/observed_otus_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization diversity_CP/observed_otus-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity diversity_CP/faith_pd_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization diversity_CP/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity diversity_CP/evenness_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization diversity_CP/evenness-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity diversity_CP/shannon_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization diversity_CP/shannon-group-significance.qzv


# Visualize alpha rarefaction plotting
qiime diversity alpha-rarefaction \
  --i-table feature_table/filtered_table_CP.qza \
  --i-phylogeny taxonomy/rooted-tree.qza \
  --p-max-depth 15000 \
  --p-steps 50 \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization diversity_CP/alpha-rarefaction-plot_15000.qzv
#The value for --p-max-depth was determined by reviewing the “Frequency per sample” information presented in the filtered_table.qza file that was created above. 

qiime diversity alpha-rarefaction \
  --i-table feature_table/filtered_table_CP.qza \
  --i-phylogeny taxonomy/rooted-tree.qza \
  --p-max-depth 5000 \
  --p-steps 50 \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization diversity_CP/alpha-rarefaction-plot_5000.qzv
#The value for --p-max-depth was determined by reviewing the “Frequency per sample” information presented in the filtered_table.qza file that was created above. 


# Visualize beta diversity significance testing
qiime diversity beta-group-significance \
  --i-distance-matrix diversity_CP/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column treatment \
  --o-visualization diversity_CP/unweighted-unifrac-treatment-group-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix diversity_CP/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column treatment \
  --o-visualization diversity_CP/weighted-unifrac-treatment-group-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix diversity_CP/jaccard_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column treatment \
  --o-visualization diversity_CP/jaccard-treatment-group-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix diversity_CP/bray_curtis_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column treatment \
  --o-visualization diversity_CP/bray_curtis-treatment-group-significance.qzv \
  --p-pairwise
