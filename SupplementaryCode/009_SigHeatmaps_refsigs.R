setwd("/Volumes/GoogleDrive/My Drive/NikZainal/projects/SignatureTools/tests")

library(signature.tools.lib)


NewCOSMIC <- read.table("../data/sigProfiler_SBS_signatures.csv", sep = ",",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
subs_sigs <- readTable("../results/2018_05_21_pancan_signatures_subs.tsv")
rownames(NewCOSMIC) <- rownames(subs_sigs)
NewCOSMIC <- NewCOSMIC[,3:ncol(NewCOSMIC)]
colnames(NewCOSMIC) <- paste0("C",substr(colnames(NewCOSMIC),4,6))
cluster_means <- signature.tools.lib::readTable("../results/pancan_extraction_clustering/2018_09_26_subs_partitions_manually_updated/partitions_mean.tsv")
subs_sigs_similarities_new <- computeCorrelationOfTwoSetsOfSigs(cluster_means,NewCOSMIC)
writeTable(subs_sigs_similarities_new,file = "../results/2018_05_21_pancan_signatures_subs_mean_vs_COSMICnew.tsv")

subs_sigs_similarities <- readTable("../results/2018_05_21_pancan_signatures_subs_mean_vs_COSMICnew.tsv")
jpg_file <- "../results/2018_05_21_pancan_signatures_subs_mean_vs_COSMICnew_inverted.jpg"

subs_sigs_similarities <- t(subs_sigs_similarities)

plot_matrix <- subs_sigs_similarities
plot_matrix[subs_sigs_similarities>0.9] <- 1
plot_matrix[subs_sigs_similarities>0.85 & subs_sigs_similarities<=0.9] <- 0.5
plot_matrix[subs_sigs_similarities<=0.85] <- 0



subs_sigs_names <- c("1","18","17","MMR1","2",
                      "MMR2","5","8","13","3",
                      "30","MIXED1","4","N1","PLATINUM",
                      "MIXED2","33","22","10","36",
                      "7","MIXED3","16","19","N2",
                      "9","N3","52","11","38",
                      "51","N4","N5","N6","N7",
                      "N8","N9","N10","N11","24",
                      "N12")

colnames(plot_matrix) <- paste0("Group ",colnames(subs_sigs_similarities),": RefSig",subs_sigs_names)

jpeg(filename = jpg_file,width = 1400,height = 2000,res = 200)
col<- colorRampPalette(c("white","blue", "red"))(51)
par(mar=c(10.5,4,4.5,2))
image(as.matrix(t(plot_matrix)),col = col,zlim = c(0,1), axes=F,main="Cosine Similarity between Substitution\nReference Signatures and COSMIC")
legend(x="topleft",legend = c("cos.sim.>0.9"),fill = "red",cex = 0.9,bty = "n",inset = c(0,-0.04),xpd = TRUE)
legend(x="topright",legend = c("0.9>=cos.sim.>0.85"),fill = "blue",cex = 0.9,bty = "n",inset = c(0,-0.04),xpd = TRUE)
extrabit <- 1/(ncol(subs_sigs_similarities)-1)/2
abline(v=seq(0-extrabit,1+extrabit,length.out = ncol(subs_sigs_similarities)+1),col="grey",lty=1)
extrabit_h <- 1/(nrow(subs_sigs_similarities)-1)/2
h_seq <-(0:(nrow(subs_sigs_similarities)))*2*extrabit_h - extrabit_h
# for (i in 1:nrow(subs_table)){
#   h_seq <- c(h_seq, h_seq[length(h_seq)] + extrabit_h*2*subs_table$nsig[i])
# }
# par(xpd=TRUE)
# for (i in 1:(length(h_seq)-1)){
#   text(-0.04,h_seq[i]+(h_seq[i+1]-h_seq[i])/2,labels=subs_table$Organ[i],adj = c(1, NA),cex=0.8)
# }
# par(xpd=FALSE)
# abline(h=h_seq,col="grey",lty=1)
# for (i in 1:nrow(subs_table)){
#   h_seq <- c(h_seq, h_seq[length(h_seq)] + extrabit_h*2*subs_table$nsig[i])
# }
abline(h=h_seq,col="grey",lty=1)
par(xpd=TRUE)
for (i in 1:(nrow(subs_sigs_similarities)-1)){
  text(-0.04,h_seq[1:(length(h_seq)-1)]+extrabit_h,labels=row.names(plot_matrix),adj = c(1, NA),cex=0.8)
}
par(xpd=FALSE)
#axis(2,at = seq(0,1,length.out = ncol(subs_sigs_similarities)),labels = colnames(subs_sigs_similarities),las=1,cex.lab=0.8)
axis(1,at = seq(0,1,length.out = ncol(subs_sigs_similarities)),labels = colnames(plot_matrix),las=2,cex.axis=0.8)
dev.off()
