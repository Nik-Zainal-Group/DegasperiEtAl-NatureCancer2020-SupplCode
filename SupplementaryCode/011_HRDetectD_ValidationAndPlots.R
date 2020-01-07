runLocal <- TRUE

if(runLocal) {
  setwd("/Volumes/GoogleDrive/My Drive/NikZainal/projects/PanCan_HRDetect/bin")
  nparallel <- 5
}else{
  nparallel <- 10
}

source("../lib/applyHRDetectDavies2017.R")
source("../lib/logisticRegression.R")
source("../lib/logisticRegression.glm.R")
source("../lib/plotHRDetect_withCI.R")

doPlotHRDetect <- TRUE

results_dir <- "../results/HRDetectTraining_741_2018_10_16/"
training_data_file <- paste0(results_dir,"2018_10_16_trainingData741.rData")
load(training_data_file)

out_dir <- paste0(results_dir,"HRDetectDaviesOnRef/")
dir.create(out_dir,showWarnings = FALSE,recursive = TRUE)

#get a ROC curve on Breast Training data also for Davies2017
BRCAprob <- applyHRDetectDavies2017(xMatrix,
                                    c("del.mh.prop",
                                      "Ref_subs_10",
                                      "Ref_rearr_9",
                                      "Ref_rearr_6",
                                      "hrd",
                                      "Ref_subs_8"),
                                    attachContributions=TRUE)
ROC.Davies2017 <- roc(yLabel,BRCAprob$Probability)
LogLoss.Davies2017 <- computeLogLoss(yLabel,BRCAprob$Probability)


#------apply trained model on additional datasets
validation.datasets <- paste0(results_dir,c("validation/scanb_validation_table.tsv",
                                            "validation/ovary_validation_table.tsv",
                                            "validation/pancreatic_validation_table.tsv"))
validation.datasets.names <- c("SCANB","Ovary","Pancreatic")
validation.ROCs.HRDetectDavies2017 <- list()
validation.LogLoss.HRDetectDavies2017 <- list()

#Davies2017
if (length(validation.datasets)>0){
  for (i in 1:length(validation.datasets)){
    dataset <- validation.datasets[i]
    d <- validation.datasets.names[i]
    val.data <- read.table(dataset,
                           sep = "\t",
                           as.is = TRUE,
                           header = TRUE,
                           check.names = FALSE)
    BRCAdef.prob <- applyHRDetectDavies2017(val.data,features_names = c("del.mh.prop","Ref_subs_10","Ref_rearr_9","Ref_rearr_6","hrd","Ref_subs_8"))
    validation.ROCs.HRDetectDavies2017[[d]] <- roc(val.data$BRCA_deficient,as.vector(BRCAdef.prob))
    validation.LogLoss.HRDetectDavies2017[[d]] <- computeLogLoss(val.data$BRCA_deficient,as.vector(BRCAdef.prob))
  }  
}



#other version of AUC plots-------------------

#get some info
Breast741_T <- sum(yLabel)
Breast741_F <- sum(!yLabel)
validation_table_scanb <- read.table(validation.datasets[1],
                                     sep = "\t",
                                     as.is = TRUE,
                                     header = TRUE,
                                     check.names = FALSE)
SCANB_T <- sum(validation_table_scanb$BRCA_deficient)
SCANB_F <- sum(!validation_table_scanb$BRCA_deficient)
validation_table_ovary <- read.table(validation.datasets[2],
                                     sep = "\t",
                                     as.is = TRUE,
                                     header = TRUE,
                                     check.names = FALSE)
ovary_T <- sum(validation_table_ovary$BRCA_deficient)
ovary_F <- sum(!validation_table_ovary$BRCA_deficient)
validation_table_pancreatic <- read.table(validation.datasets[3],
                                          sep = "\t",
                                          as.is = TRUE,
                                          header = TRUE,
                                          check.names = FALSE)
pancreatic_T <- sum(validation_table_pancreatic$BRCA_deficient)
pancreatic_F <- sum(!validation_table_pancreatic$BRCA_deficient)

#plotting AUROC of all the models
colours_plot <- c("red","blue","green")
set_mar <- c(3.5, 4, 2, 2.5)
set_mgp <- c(2,0.8,0)
jpeg(filename = paste0(out_dir,"AUCs_HRDetectDaviesOnRefSig.jpg"),
     width = 550,
     height = 600,
     res = 150)
par(mfrow=c(1,1),oma=c(3.5,0,0,0),mar=c(4,4,3.5,2.5),mgp = c(2,0.8,0))
#HRDetect Davies2017
if (length(validation.datasets)>0){
  d <- validation.datasets.names[1]
  plot(ROC.Davies2017,mar=set_mar,mgp = set_mgp,col = "black")
  AUCs <- c(ROC.Davies2017$auc)
  if (length(validation.datasets)>1){
    for(i in 1:length(validation.datasets.names)){
      d <- validation.datasets.names[i]
      plot(validation.ROCs.HRDetectDavies2017[[d]],add=TRUE,col = colours_plot[i])
      AUCs <- c(AUCs,validation.ROCs.HRDetectDavies2017[[d]]$auc)
    }
  }
}
legend("bottomright",legend = sprintf("%0.3f",AUCs),title = "AUC",bty = "n",cex = 0.83,
       lty=c(rep(1,length(validation.datasets.names)+1)), col = c("black",colours_plot),lwd = 2)
title(main="HRDetect (Davies2017)",line = 2,cex.main=0.9)

#add overall title and legend
#title(paste0("Performance of HRDetect compared to single features (Breast741 dataset)"), outer=TRUE)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
set_legend <- c(paste0("Breast (T=",Breast741_T,"/F=",Breast741_F,")"),
                paste0("SCANB (T=",SCANB_T,"/F=",SCANB_F,")"),
                paste0("Ovary (T=",ovary_T,"/F=",ovary_F,")"),
                paste0("Pancreatic (T=",pancreatic_T,"/F=",pancreatic_F,")"))
legend("bottom", legend = set_legend[3:4], xpd = TRUE, horiz = TRUE, 
       inset = c(0,0.02), bty = "n", lty=c(1), col = c(colours_plot[2:3]), cex = 0.8,lwd = 2)
legend("bottom", legend = set_legend[1:2], xpd = TRUE, horiz = TRUE, 
       inset = c(0,0.08), bty = "n", lty=c(1), col = c("black",colours_plot[1]), cex = 0.8,lwd = 2)
dev.off()


#-------------------------------------------
#Tissue wise HRDetect scores (no data uncertainty considered)
#-------------------------------------------

#features names for HRDetect Davies et al. 2017
features_names <- c("del.mh.prop","Ref_subs_10","Ref_rearr_9","Ref_rearr_6","hrd","Ref_subs_8")

#get the organs to test
HRDetect_files_dir <- "../data/HRDetect_data_files/"
samples_files <- list.files(path = HRDetect_files_dir)
organs <- unique(sapply(samples_files,function(x) paste0(strsplit(x,"[.]")[[1]][2:2],collapse = "")))
organs <- setdiff(organs,"Breast")

scores_tables <- list()
full_table <- NULL

groups <- organs
for (group in groups){ 
  message("Organ: ",group)
  
  hrdet_table <- read.table(paste0("../data/HRDetect_data_files/HRDetect_data_table.",group,".vRefSig.tsv"),
                            sep = "\t",header = TRUE,check.names = FALSE,as.is = TRUE)
  hrdet_table[is.na(hrdet_table)] <- 0
  scores_tables[[group]] <- applyHRDetectDavies2017(hrdet_table,features_names = features_names,attachContributions = TRUE)

  #build up the full table
  tmptable <- cbind(hrdet_table[,features_names],scores_tables[[group]][,"Probability",drop=FALSE])
  colnames(tmptable) <- c("del.mh.prop","Ref.Sig.3","Ref.Sig.R3","Ref.Sig.R5","HRD-LOH","Ref.Sig.8","HRDetect single score")
  if (group=="Breast741") group <- "Breast"
  tmptable$organ <- group
  full_table <- rbind(full_table,tmptable)
}



#-------------------------------------------
#Tissue wise HRDetect scores (with data uncertainty considered)
#-------------------------------------------

nboots_hr <- 1000
outdir <- "../results/HRDetectTraining_741_2018_10_16/HRDetect_with_Confidence/"
paste0(outdir,"bootstrap_tables_nboots",nboots_hr,".rData")
#bootstrap_tables_file
bootstrap_tables_file <- paste0(outdir,"bootstrap_tables_nboots",nboots_hr,".rData")
load(bootstrap_tables_file)


#features names for HRDetect Davies et al. 2017
features_names <- c("del.mh.prop","Ref_subs_10","Ref_rearr_9","Ref_rearr_6","hrd","Ref_subs_8")

#bootstrap_HRDetect results file
bootstrap_scores_file <- paste0(out_dir,"bootstrap_HRDetectDaviesOnRef_scores_nboots",nboots_hr,".rData")

if(!file.exists(bootstrap_scores_file)){
  current_scores_table <- list()
  for (organ in organs){
    message("Organ: ",organ)
    #for all bootstrapped tables
    current_scores_table[[organ]] <- list()
    current_scores_table[[organ]][["Davies2017"]] <- NULL
    for (i in 1:nboots_hr){
      cat(".")
      #run the models
      current_table <- bootstrap_tables[[organ]][[i]]
      current_table[is.na(current_table)] <- 0
      #HRDetect Davies et al. 2017
      current_scores <- t(applyHRDetectDavies2017(current_table,features_names,attachContributions = FALSE))
      current_scores_table[[organ]][["Davies2017"]] <- rbind(current_scores_table[[organ]][["Davies2017"]],current_scores)
    }
    cat("\n")
  }
  #save the bootstrap scores
  save(file = bootstrap_scores_file,current_scores_table,organs,nboots_hr)
}else{
  load(bootstrap_scores_file)
  message("Bootstraps results loaded from file")
}

#update the table
newcols <- NULL
for (organ in organs){
  q_5_50_95 <- t(apply(current_scores_table[[organ]][["Davies2017"]],2,function(x) quantile(x,c(0.05,0.5,0.95))))
  colnames(q_5_50_95) <- paste("HRDetect bootstrap score percentile",colnames(q_5_50_95))
  newcols <- rbind(newcols,q_5_50_95)
}
full_table <- cbind(full_table,newcols[row.names(full_table),])
write.table(full_table,file = paste0(out_dir,"HRDetect_summary_table",nboots_hr,".tsv"),
            quote = FALSE,col.names = TRUE,row.names = TRUE,sep = "\t")


#-------------------------------------------
#Plotting all
#-------------------------------------------

#load BRCA deficient samples
brca_status_table_v3 <- read.table("../data/BRCA_status/pcawg_samples_BRCA_status_v3_20171115.csv",header = TRUE,sep = "\t",as.is = TRUE)
biallelic_samples <- unique(brca_status_table_v3$sample_id[brca_status_table_v3$biallelic=="yes"])
biallelic_samples <- union(biallelic_samples,rownames(summaryDF)[summaryDF$BRCA_deficient==TRUE])

if(doPlotHRDetect){
  for (organ in organs){
    message("Organ: ",organ)
    #HRDetect Davies et al. 2017
    #boxplot(current_scores_table[,order(apply(current_scores_table,2,median))])
    q_5_50_95 <- t(apply(current_scores_table[[organ]][["Davies2017"]],2,function(x) quantile(x,c(0.05,0.5,0.95))))
    filename <- paste0(out_dir,"HRDetect_Score_bootstrap_",organ,"_Davies2017OnRefSigs.jpg")
    par_title <- paste0("HRDetect Score Davies et al. 2017, ",organ,", median and 5%-95% quantiles interval (n=",nboots_hr,")")
    plot_HRDetect_withCI(filename,par_title,q_5_50_95,samplesBRCAdef = biallelic_samples,single_score = scores_tables[[organ]])
  }
}

#plot also the training data alone:
organ <- "Breast741"
q_5_50_95 <- t(apply(current_scores_table[[organ]][["Davies2017"]],2,function(x) quantile(x,c(0.05,0.5,0.95))))
q_5_50_95 <- q_5_50_95[row.names(xMatrix),]
single_score_tab <- scores_tables[[organ]]
single_score_tab <- single_score_tab[row.names(xMatrix),]
sens <- sum(yLabel & q_5_50_95[,"5%"]>0.5)/sum(yLabel)
spec <- sum(!yLabel & q_5_50_95[,"5%"]<=0.5)/sum(!yLabel)
sens_single <- sum(yLabel & single_score_tab[,"Probability"]>0.7)/sum(yLabel)
spec_single <- sum(!yLabel & single_score_tab[,"Probability"]<=0.7)/sum(!yLabel)
filename <- paste0(out_dir,"HRDetect_Score_bootstrap_BreastTrainingData_Davies2017OnRefSigs.jpg")
par_title <- paste0("HRDetect Davies et al. 2017, Breast Training data, median and 5%-95% quantiles interval (n=",nboots_hr,")")
plot_HRDetect_withCI(filename,par_title,q_5_50_95,samplesBRCAdef = biallelic_samples,single_score = single_score_tab)

#---------------------------------
#--- plot some individual smaples of interest to explain CI
#----------------------------------

# library(vioplot)

#Let's find some BC samples that have barely passed the 0.7 threshold but have large CI
selection <- single_score_tab[,"Probability"] > 0.7 & q_5_50_95[,"5%"] <= 0.5
BRCA_flag <- yLabel[selection]
sample_names <- rownames(single_score_tab)[selection]

organ <- "Breast741"
samples_to_plot <- c("PD11346a","PD24186a","fca6150f-d555-a29e-e040-11ac0d4873b2","PD22355a")
#samples_to_plot_lab <- paste0(substr(samples_to_plot,1,11),c("","","...",""))
samples_to_plot_lab <-c("PD11346a","PD24186a","PD13622a","PD22355a")
# vioplot(current_scores_table[[organ]][["Davies2017"]][,samples_to_plot[1]],
#         current_scores_table[[organ]][["Davies2017"]][,samples_to_plot[2]],
#         current_scores_table[[organ]][["Davies2017"]][,samples_to_plot[3]],
#         current_scores_table[[organ]][["Davies2017"]][,samples_to_plot[4]],col = "#9999FF",
#         names = samples_to_plot,rectCol = "grey",pchMed = 4,colMed = "black")

filename <- paste0(out_dir,"HRDetect_Score_Examples_Davies2017OnRefSigs.jpg")
jpeg(filename = filename,width = 1500,height = 900,res = 250)
par(mar=c(7,4,2,2))
#Maybe the usual boxplot is better
hrdet_score <- current_scores_table[[organ]][["Davies2017"]][,samples_to_plot]
boxplot(hrdet_score,main="Distribution of HRDetect score in selected samples",
        ylab="HRDetect score",las=3, names=samples_to_plot_lab,ylim=c(0,1))
thisvalues <- hrdet_score
# take the x-axis indices and add a jitter, proportional to the N in each level
for(i in 1:ncol(hrdet_score)){
  myjitter<-jitter(rep(i, length(thisvalues[,i])), amount=0.2)
  points(myjitter, thisvalues[,i], pch=20, col=rgb(1,0,0,.1))
}
points(scores_tables[[organ]][samples_to_plot,"Probability"],pch=23,col="black",bg="yellow")
dev.off()


#-------------------------------------------
# Compare Median score with Sig3 (Ref10)
#-------------------------------------------

filename <- paste0(out_dir,"HRDetect_vs_Sig3_Boxplot_Davies2017OnRefSigs.jpg")
jpeg(filename = filename,width = 2100,height = 900,res = 150)
par(mfrow=c(3,7))
for (organ in organs){
  message("Organ: ",organ)
  
  #table with no uncertainty, to get the consensus median Ref10 exposure
  hrdet_table <- read.table(paste0("../data/HRDetect_data_files/HRDetect_data_table.",organ,".vRefSig.tsv"),
                            sep = "\t",header = TRUE,check.names = FALSE,as.is = TRUE)
  hrdet_table[is.na(hrdet_table)] <- 0
  #compute median HRDetect score
  q_5_50_95 <- t(apply(current_scores_table[[organ]][["Davies2017"]],2,function(x) quantile(x,c(0.05,0.5,0.95))))
  
  #groups
  g1 <- hrdet_table[q_5_50_95[,"50%"]<0.5,"Ref_subs_10"]
  g2 <- hrdet_table[q_5_50_95[,"50%"]>=0.5,"Ref_subs_10"]
  howmanyg1 <- length(g1)
  howmanyg2 <- length(g2)
  
  boxplot(g1,
       g2,
       xlab = "HRDetect Median Score",ylab = "Ref_subs_10 (mutations)",
       main = paste0("",organ),cex.main = 0.9,pch=1,names = c("<0.5",">=0.5"),ylim=c(0,max(g1,g2)))
  myjitter <- jitter(rep(1, howmanyg1), amount=0.2)
  points(myjitter, g1, pch=20, col=rgb(1,0,0,.1))
  myjitter <- jitter(rep(2, howmanyg2), amount=0.2)
  points(myjitter, g2, pch=20, col=rgb(1,0,0,.5))
  
  
}
dev.off()


#-------------------------------------------
# How many high confidence do we have in each organ?
#-------------------------------------------

#calculate the high confidence in percent
#compare with number of samples with high HRDetect score (>0.7 as we historically do)
#compare with number of samples that contain Signature 3

data_names <- c("HRDetect high confidence","HRDetect high single score","RefSig 3 present","RefSig 3 >= 100 muts","RefSig 3 >= 1000 muts")
barplot_table <- matrix(NA,nrow = length(organs),ncol = length(data_names),dimnames = list(organs,data_names))

nsamples <- c()
for (organ in organs){
  message("Organ: ",organ)
  
  #table with no uncertainty, to get the consensus median Ref10 exposure
  hrdet_table <- read.table(paste0("../data/HRDetect_data_files/HRDetect_data_table.",organ,".vRefSig.tsv"),
                            sep = "\t",header = TRUE,check.names = FALSE,as.is = TRUE)
  hrdet_table[is.na(hrdet_table)] <- 0
  #compute median HRDetect score
  q_5_50_95 <- t(apply(current_scores_table[[organ]][["Davies2017"]],2,function(x) quantile(x,c(0.05,0.5,0.95))))
  
  
  barplot_table[organ,"HRDetect high confidence"] <- sum(q_5_50_95[,"5%"]>0.5)/nrow(q_5_50_95)*100
  barplot_table[organ,"HRDetect high single score"] <- sum(scores_tables[[organ]][,"Probability"]>0.7)/nrow(scores_tables[[organ]])*100
  barplot_table[organ,"RefSig 3 present"] <- sum(hrdet_table[,"Ref_subs_10"]>0)/nrow(hrdet_table)*100
  barplot_table[organ,"RefSig 3 >= 100 muts"] <- sum(hrdet_table[,"Ref_subs_10"]>=100)/nrow(hrdet_table)*100
  barplot_table[organ,"RefSig 3 >= 1000 muts"] <- sum(hrdet_table[,"Ref_subs_10"]>=1000)/nrow(hrdet_table)*100
  nsamples <- c(nsamples,nrow(hrdet_table))
}
names(nsamples) <- organs

data_to_plot <- barplot_table[,c(1,2,3,5)]

jpeg(filename = paste0(out_dir,"HRDetect_andSig3_pancan.jpg"),
     width = 3000,
     height = 900,
     res = 200)
par(xpd=FALSE,mar=c(7.5,4.1,5.1,14.5))
colours_lab <- c("white","lightgrey","orange","blue","black","green","red","yellow","purple","pink","lightblue")[1:ncol(data_to_plot)]
barplot(t(data_to_plot),
        ylim = c(0,100),beside = TRUE,las = 2,names.arg = rep("",nrow(data_to_plot)))
abline(h=c(seq(10,90,10)),lty=3, col="lightgrey")
barplot(t(data_to_plot),
        main = "Percentage of HRDetect high score samples and Signature 3 samples across organs",
        beside = TRUE,las = 2,
        col = colours_lab,
        ylab = "Samples (%)",
        legend.text = colnames(data_to_plot),args.legend = list(x = "right", bty='n',inset=c(-0.25,0)),add = TRUE)
par(xpd=TRUE)
startpoint <- 3
for (i in nsamples){
  text(startpoint,110,i,col = "red")
  startpoint <- startpoint + 5
}
text(startpoint + 15,110,"total samples",col = "red")
rect(1,105,135,115,border = "red")
par(xpd=FALSE)
dev.off()

#same but in boxplot
jpeg(filename = paste0(out_dir,"HRDetect_andSig3_pancan_Boxplot.jpg"),
     width = 1000,
     height = 1100,
     res = 190)
par(xpd=FALSE,mar=c(11,4.1,5.5,2.1))
boxplot(data_to_plot,ylim=c(0,120),
        las=2,
        border=rgb(0,0,0.8,1),
        ylab = "Samples (%)",boxlwd = 2,whisklwd = 2,yaxt="n",
        names = c("HRDetect\nhigh confidence","HRDetect\nhigh single score",
                  "RefSig 3\npresent","RefSig 3\n>= 1000 muts"),
        main = "Percentage of HRDetect high score samples\n and Signature 3 samples across organs")
axis(2,at=seq(0,100,20),labels = seq(0,100,20))
for(i in 1:ncol(data_to_plot)){
  # thislevel<-mylevels[i]
  thisvalues <- data_to_plot[,i]
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter<-jitter(rep(i, length(thisvalues)), amount=0.2)
  points(myjitter, thisvalues, pch=20, col=rgb(1,0,0,.4)) 
}
#p-values
#
pg1 <- 1
pg2 <- 3
pheight <- 105
lines(x=c(pg1,pg1,pg2,pg2),y=c(pheight-3,pheight,pheight,pheight-3))
g1 <- data_to_plot[,pg1]
g2 <- data_to_plot[,pg2]
t.res <- t.test(x = g1,y = g2,alternative = "two.sided",paired = TRUE)
gcor <- cor(g1,g2)
text(x = (pg1+pg2)/2,y = pheight+3,labels=sprintf("p=%.1e",t.res$p.value),cex = 0.6,font = 2)
#
pg1 <- 1
pg2 <- 4
pheight <- 117
lines(x=c(pg1,pg1,pg2,pg2),y=c(pheight-3,pheight,pheight,pheight-3))
g1 <- data_to_plot[,pg1]
g2 <- data_to_plot[,pg2]
t.res <- t.test(x = g1,y = g2,alternative = "two.sided",paired = TRUE)
gcor <- cor(g1,g2)
text(x = (pg1+pg2)/2,y = pheight+3,labels=sprintf("p=%.1e",t.res$p.value),cex = 0.6,font = 2)
#
pg1 <- 1
pg2 <- 2
pheight <- 85
lines(x=c(pg1,pg1,pg2,pg2),y=c(pheight-3,pheight,pheight,pheight-3))
g1 <- data_to_plot[,pg1]
g2 <- data_to_plot[,pg2]
t.res <- t.test(x = g1,y = g2,alternative = "two.sided",paired = TRUE)
gcor <- cor(g1,g2)
text(x = (pg1+pg2)/2,y = pheight+3,labels=sprintf("p=%.1e",t.res$p.value),cex = 0.6,font = 2)
#legend
par(xpd=TRUE)
text(x = 2,y = 130,"p is the p-value of paired two-sided t-test",cex = 0.9)
#text(x = 2.5,y = -130,"model")
par(xpd=FALSE)
dev.off()


