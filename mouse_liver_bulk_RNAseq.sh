#!/bin/bash
# Set the working directory
setwd("path/to/your/project_folder")
dir <- getwd()
dir
 
# Load required library
library("sleuth")
library("dplyr")
library(ggplot2)
library(tidyr)
 
 
sample_id <- dir(file.path("Kallisto_Output"))
sample_id
# A list of paths to the kallisto results indexed by the sample IDs is collated with
kal_dirs <- file.path(dir, "Kallisto_Output", sample_id)
kal_dirs
 
#Here we are creating a column to label the samples. h4 is just the name of the label, it does not mean anything
#h4 <- c("H4", "H4", "H4", "C", "C")
h4 <- c("C", "C", "T", "T", "C", "C", "T", "T")
 
# The next step is to load an auxillary table that describes the experimental 
# design and the relationship between the kallisto directories and the samples:
 
s2c <- as.data.frame(cbind(sample_id, h4))
s2c
 
# Now the directories must be appended in a new column to the table describing the experiment.
# This column must be labeled path, otherwise sleuth will report an error. 
# This is to ensure that samples can be associated with kallisto quantifications.
s2c <- dplyr::mutate(s2c, path = kal_dirs)
# It is important to check that the pairings are correct
s2c
colnames(s2c) <- c("sample", "h4", "path")
s2c
 
group <- paste(s2c$h4,sep=".")
group
 
s2c <- dplyr::mutate(s2c, group=group)
s2c
 
#12.26.24 I keep getting an error when looking at the model set-up, so to fix it, we must re-level so that WT is set as baseline
# We must do this before setting up sleuth object
 
s2c$group <- factor(s2c$group, levels = c("C", "T"))
s2c$group <- relevel(s2c$group, ref = "C")
 
 
# Next, the “sleuth object” can be constructed. This object will store not only the information about the experiment,
# but also details of the model to be used for differential testing, and the results.
# It is prepared and used with four commands that
 
#annotate transcript to gene
 
# collect gene names
mouse_transcripts_df <- read.table("path/to/transcripts_to_genes_mouse.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
 
colnames(mouse_transcripts_df) <- c("target_id", "ens_gene", "ext_gene")
head(mouse_transcripts_df)
 
# (1) load the kallisto processed data into the object 
so <- sleuth_prep(s2c, target_mapping = mouse_transcripts_df, gene_mode = TRUE, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE, num_cores = 1, verbose = TRUE)
 
# 1.7.25 Creating PCA plot
plot_pca(so, color_by = 'group')
plot_pca(so, color_by = 'sample')
plot_pc_variance(so) # this will tell you the pca variance percentages
 
plot_group_density(so, use_filtered = TRUE, units = "est_counts",
                   trans = "log", grouping = setdiff(colnames(so$sample_to_covariates),
                                                     "sample"), offset = 1)
 
plot_sample_heatmap(so, use_filtered = TRUE, color_high = "dodgerblue",
                    color_low = "white", x_axis_angle = 50,
                    annotation_cols = setdiff(colnames(so$sample_to_covariates), "sample"),
                    cluster_bool = TRUE)
 
# (2) estimate parameters for the sleuth response error measurement (full) model 
so <- sleuth_fit(so, ~group, 'full')
 
# (3) estimate parameters for the sleuth reduced model 
so <- sleuth_fit(so, ~1, 'reduced')
 
# (4) perform differential analysis (testing) using the likelihood ratio test.
so <- sleuth_lrt(so, 'reduced', 'full')
 
# The models that have been fit can always be examined with the models() function.
models(so)
tests(so)
 
##### WALD TEST
# Get the results for the comparison of interest (run these one a
#so <- sleuth_wt(so, which_beta=as.character('groupDMSO'))
so <- sleuth_wt(so, which_beta=as.character('groupT'))
#so <- sleuth_wt(so, which_beta =as.character('groupGLP_5uM'))
 
#12.26.24 to see the results of the volcano plot, we will use this command
# since we want both the significant and non-significant to show, we will use show_all = T
 
results_group_mouse_wt_liver <- sleuth_results(so, "groupT", test_type = "C", show_all = T)

write.csv(results_group_mouse_wt_liver , "2026_4_27_CCL4_WT.csv")
 
# 1.7.25 Creating PCA plot
png("PCA_plot.png", width = 2000, height = 1500, res = 300)
plot_pca(so, color_by = "group")
dev.off()
 
plot_group_density(so, use_filtered = TRUE, units = "est_counts",
                   trans = "log", grouping = setdiff(colnames(so$sample_to_covariates),
                                                     "sample"), offset = 1)
 
png("sample_heatmap_plot.png", width = 2000, height = 1500, res = 300)
plot_sample_heatmap(so, use_filtered = TRUE, color_high = "dodgerblue",
                    color_low = "white", x_axis_angle = 50,
                    annotation_cols = setdiff(colnames(so$sample_to_covariates), "sample"),
                    cluster_bool = TRUE)
