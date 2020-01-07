setwd("~/sandbox/git/signature-tools/tests")

library(signature.tools.lib)


cluster_means <- signature.tools.lib::readTable("../results/pancan_extraction_clustering/2018_09_26_subs_partitions_manually_updated/partitions_mean.tsv")
cos_sim_means <- computeCorrelation(cluster_means)
cos_sim_COSMIC <- computeCorrelation(COSMIC30_subs_signatures)

all_pancan_sigs <- signature.tools.lib::readTable("../results/2018_05_21_pancan_signatures_subs.tsv")
cos_sim_pancan <- computeCorrelation(all_pancan_sigs)

#sig3 and 8
cos_sim_means[8,10]
cos_sim_COSMIC[3,8]
cos_sim_pancan["Breast741_11","Breast741_5"]
cos_sim_pancan["Ovary_7","Ovary_5"]
cos_sim_pancan["Pancreas_5","Pancreas_11"]

#not sure this is a good direction


