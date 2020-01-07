runLocally <- TRUE

if (runLocally){
  setwd("~/sandbox/git/signature-tools/tests")
  resultdir <- "../results/testsSimulatedData/"
  ngenomes <- 30
  nparallel <- 5
}else{
  args = commandArgs(trailingOnly=TRUE)
  ngenomes <- as.numeric(args[1])
  resultdir <- "~/scratch116/testsSimulatedData/"
  nparallel <- 10
}

source("../lib/SignatureExtractionLib.R")
source("../lib/SignatureFitLib.R")
library("lpSolve")

recompute_Efit_table <- FALSE
fit_method <- "bootstrap" #bootstrap or bleedingfilter

n <- 5
nsig_range <- 8:12
nmf_algorithms <- c("brunet","lee","nsNMF")
nmf_algorithms_labels <- c("Lee KLD","Lee Frobenius","nsNMF")
#nmf_algorithms <- nmf_algorithms[1:2]
#nmf_algorithms_labels <- nmf_algorithms_labels[1:2]

filter_labels <- c("filter","no filter")
clustering_plotlabels <- c("HC","CM","PAM")
clustering_labels <- c("HC","MC","PAM") #hierarchical clustering or matched clustering
nsettings <- length(clustering_labels)*length(filter_labels)

#Cosmic Signatures
cosmic30 <- read.table("../data/COSMIC_signatures.txt",sep = "\t",header = TRUE,as.is = TRUE,check.names = FALSE)
cosmic30 <- sortCatalogue(cosmic30)
signatures <- c(1,2,3,5,6,8,12,13,17,18)
P <- cosmic30[,signatures]

#parameters
nsignatures <- ncol(P)
nsignatures_per_sample <- 5
min_mutations <- 1000
max_mutations <- 50000

#set number of bootstraps
nboots_withFilter <- 20
nboots_withoutFilter <- 1000

alphas <- c(-1,0.005,0.01,0.02,0.04)
# alpha_sig <- alphas[2]
alphas_KLD <- c(-1,1,2,4)

#Signature Fit with Bootstraps
thresholds_perc <- c(0,1,2,5,10)
nboots_fit <- 100
p.value <- 0.05

#change directory
resultdir <- paste0(resultdir,"sigPerSample",nsignatures_per_sample,"_samples",ngenomes,"/")

#read the results from files (stats and signatures)

all_labels <- paste(rep(nmf_algorithms_labels,each = nsettings*n),
                    rep(rep(filter_labels,each = n),length(nmf_algorithms)*length(clustering_labels)),
                    rep(rep(clustering_labels,each = n*length(filter_labels)),length(nmf_algorithms)),
                    paste0("rep",rep(1:n,length(nmf_algorithms)*nsettings)),
                    sep=" ")

signatures_similarity_table <- matrix(NA,nrow = nsignatures, ncol = length(nmf_algorithms)*nsettings*n)
row.names(signatures_similarity_table) <- paste0("S",1:nsignatures)
colnames(signatures_similarity_table) <- all_labels

#ASW
ASW_table <- matrix(NA,nrow = length(nsig_range), ncol = length(nmf_algorithms)*nsettings*n)
row.names(ASW_table) <- as.character(nsig_range)
colnames(ASW_table) <- all_labels

if (fit_method=="bleedingfilter"){
  #Exposure fit (BF CosSim)
  Efit_RMSE_table <- matrix(NA,nrow = length(alphas), ncol = length(nmf_algorithms)*nsettings*n)
  row.names(Efit_RMSE_table) <- as.character(alphas)
  colnames(Efit_RMSE_table) <- all_labels
}else if(fit_method=="bootstrap"){
  Efit_RMSE_table <- matrix(NA,nrow = length(thresholds_perc), ncol = length(nmf_algorithms)*nsettings*n)
  row.names(Efit_RMSE_table) <- as.character(thresholds_perc)
  colnames(Efit_RMSE_table) <- all_labels
}

if (fit_method=="bleedingfilter"){
  Efit_file <- paste0(resultdir,"ExposureFit_CosSim.tsv")
}else if(fit_method=="bootstrap"){
  Efit_file <- paste0(resultdir,"ExposureFit_bootstraps.tsv")
}

# #Exposure fit (BF KLD)
# Efit_RMSE_table_KLD <- matrix(NA,nrow = length(alphas_KLD), ncol = length(nmf_algorithms)*nsettings*n)
# row.names(Efit_RMSE_table_KLD) <- as.character(alphas_KLD)
# colnames(Efit_RMSE_table_KLD) <- all_labels
# 
# Efit_file_KLD <- paste0(resultdir,"ExposureFit_KLD.tsv")

#Error
Error_table <- matrix(NA,nrow = length(nsig_range), ncol = length(nmf_algorithms)*length(filter_labels)*n)
row.names(Error_table) <- as.character(nsig_range)
colnames(Error_table) <- paste(rep(nmf_algorithms_labels,each = length(filter_labels)*n),
                               rep(rep(filter_labels,each = n),length(nmf_algorithms)),
                               paste0("rep",rep(1:n,length(nmf_algorithms)*length(filter_labels))),
                               sep=" ")

#Classification of mutations and sigantures
listOfClassificationMetrics <- c("FP_mut","FP_sig","FN_mut","FN_sig")
listOfClassificationMetrics_labels <- c("FP rate of mutation assignment","FP rate of signature assignment","FN rate of mutation assignment","FN rate of signature assignment")
listOfClassificationMetricsSensSpec <- c("Specificity_mut","Specificity_sig","Sensitivity_mut","Sensitivity_sig")
listOfClassificationMetricsSensSpec_labels <- c("Specificity of mutation assignment","Specificity of signature assignment to samples","Sensitivity of mutation assignment","Sensitivity of signature assignment to samples")
listOfClassificationMetrics_tables <- list()
listOfClassificationMetrics_files <- list()

if (fit_method=="bleedingfilter"){
  for (m in listOfClassificationMetrics){
    listOfClassificationMetrics_tables[[m]] <- matrix(NA,nrow = length(alphas), ncol = length(nmf_algorithms)*nsettings*n)
    row.names(listOfClassificationMetrics_tables[[m]]) <- as.character(alphas)
    colnames(listOfClassificationMetrics_tables[[m]]) <- all_labels
    listOfClassificationMetrics_files[[m]] <- paste0(resultdir,"ExposureClassification_",m,".tsv")
  }
}else if(fit_method=="bootstrap"){
  for (m in listOfClassificationMetrics){
    listOfClassificationMetrics_tables[[m]] <- matrix(NA,nrow = length(thresholds_perc), ncol = length(nmf_algorithms)*nsettings*n)
    row.names(listOfClassificationMetrics_tables[[m]]) <- as.character(thresholds_perc)
    colnames(listOfClassificationMetrics_tables[[m]]) <- all_labels
    listOfClassificationMetrics_files[[m]] <- paste0(resultdir,"ExposureClassification_",m,"_bootstraps.tsv")
  }
}

#Classification of mutations and signatures, individual signatures
listOfClassificationMetrics_IndSig_tables <- list()
listOfClassificationMetrics_IndSig_files <- list()

if (fit_method=="bleedingfilter"){
  for(alpha in alphas){
    listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]] <- list()
    listOfClassificationMetrics_IndSig_files[[as.character(alpha)]] <- list()
    for (m in listOfClassificationMetrics){
      listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][[m]] <- matrix(NA,nrow = length(1:nsignatures), ncol = length(nmf_algorithms)*nsettings*n)
      row.names(listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][[m]]) <- as.character(1:nsignatures)
      colnames(listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][[m]]) <- all_labels
      listOfClassificationMetrics_IndSig_files[[as.character(alpha)]][[m]] <- paste0(resultdir,"ExposureClassification_IndSig_",m,"_a",alpha,".tsv")
    }
  }
}else if(fit_method=="bootstrap"){
  for(alpha in thresholds_perc){
    listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]] <- list()
    listOfClassificationMetrics_IndSig_files[[as.character(alpha)]] <- list()
    for (m in listOfClassificationMetrics){
      listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][[m]] <- matrix(NA,nrow = length(1:nsignatures), ncol = length(nmf_algorithms)*nsettings*n)
      row.names(listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][[m]]) <- as.character(1:nsignatures)
      colnames(listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][[m]]) <- all_labels
      listOfClassificationMetrics_IndSig_files[[as.character(alpha)]][[m]] <- paste0(resultdir,"ExposureClassification_IndSig_",m,"_a",alpha,"_bootstraps.tsv")
    }
  }
}


#just change alphas if necessary
if (fit_method=="bleedingfilter"){
  #do nothing, alphas are alphas
}else if(fit_method=="bootstrap"){
  alphas <- thresholds_perc
}

pos <- 1
pos_error <- 1
for (algo in nmf_algorithms){ 
  message("algo ",algo)
  for(cl in clustering_labels){
    message("clustering ",cl)
    
    # if (cl=="MC"){
    #   clusteringLabel <- "MC"
    # }else{ #HC
    #   clusteringLabel <- ""
    # }
    
    #---With Filter---
    for (i in 1:n){
    
      #read catalogues
      cat_file <- paste0(resultdir,"simulated_catalogue_",i,".tsv")
      E_file <- paste0(resultdir,"simulated_E_",i,".tsv")
      P_file <- paste0(resultdir,"simulated_P_",i,".tsv")
      if (file.exists(cat_file)){
        catalogue <- read.table(file = cat_file,sep = "\t",header = TRUE,check.names = FALSE,as.is = TRUE)
        E <- read.table(file = E_file,sep = "\t",header = TRUE,check.names = FALSE,as.is = TRUE)
        P <- read.table(file = P_file,sep = "\t",header = TRUE,check.names = FALSE,as.is = TRUE)
        message("catalogue ",i," loaded from file")
      } else {
        message("catalogue",i,"not found")
      }
      
      outFilePath_withFilterBest=paste0(resultdir,"with_filterBest_",algo,"_",i,"/")
      
      sig_file <- paste0(outFilePath_withFilterBest,"round_1/sig_",nsignatures,"/Sigs_plot_SimulatedData_ns",nsignatures,"_nboots",nboots_withFilter,"_",cl,".tsv")
      sigs <- read.table(sig_file,sep = "\t",check.names = FALSE,header = TRUE,as.is = TRUE)
      sigs.cor <- as.matrix(computeCorrelationOfTwoSetsOfSigs(sigs,P))
      res.match <- lp.assign(1 - sigs.cor)
      cosmic_match <- c()
      for(m in 1:nsignatures){
        cosmic_match <- c(cosmic_match,which.max(res.match$solution[,m]))
        signatures_similarity_table[m,pos] <- sigs.cor[m,which.max(res.match$solution[m,])]     
      }
      names(cosmic_match) <- colnames(sigs.cor)
      
      #sigFit and compute accuracy
      if(recompute_Efit_table | !file.exists(Efit_file)){
        for (alpha in alphas){
          message("threshold ",alpha)
          if (fit_method=="bleedingfilter"){
            E_hat <- SignatureFit(catalogue,signature_data_matrix = sigs[,cosmic_match],method = "KLD",bf_method = "CosSim",alpha = alpha)
          }else if(fit_method=="bootstrap"){
            file_store <- paste0(outFilePath_withFilterBest,"round_1/sig_",nsignatures,"/SigFit_withBootstrap_Summary_m","KLD","_bfm","CosSim","_alpha",-1,"_tr",alpha,"_p",p.value,"_",cl,".rData")
            if(file.exists(file_store)){
              load(file_store)
              message("Bootstrap Signature Fits loaded from file")
            }else{
              res <- SignatureFit_withBootstrap(catalogue,signature_data_matrix = sigs[,cosmic_match],method = "KLD",nboot = nboots_fit,threshold_percent = alpha,threshold_p.value = p.value,verbose = FALSE,nparallel = nparallel)
              save(file = file_store,res)
            }
            E_hat <- res$E_median_filtered
          }
          Efit_RMSE_table[as.character(alpha),pos] <- sqrt(sum((E-E_hat)^2)/(ncol(E)*nrow(E)))
          
          #classification metrics
          Ediff <- round(E)-E_hat
          #compute proportions relative to each sample number of mutations
          Ediff <- Ediff/matrix(rep(apply(E,2,sum),nsignatures),nrow = nsignatures,byrow = TRUE)
          #proportion of mutations per sample (i.e. divide by ngenomes)
          listOfClassificationMetrics_tables[["FP_mut"]][as.character(alpha),pos] <- sum(abs(Ediff[Ediff<0]))/ngenomes
          listOfClassificationMetrics_tables[["FN_mut"]][as.character(alpha),pos] <- sum(abs(Ediff[Ediff>0]))/ngenomes
          #proportion of mutations per sample (i.e. divide by ngenomes), proportion w.r.t actual positives and negatives in E
          listOfClassificationMetrics_tables[["FP_sig"]][as.character(alpha),pos] <- (sum(E==0 & E_hat>0)/ngenomes)/(nsignatures-nsignatures_per_sample)
          listOfClassificationMetrics_tables[["FN_sig"]][as.character(alpha),pos] <- (sum(E>0 & E_hat==0)/ngenomes)/(nsignatures_per_sample)
          
          #classification of individual signatures 
          FP_mut <- rep(0,nsignatures)
          FN_mut <- rep(0,nsignatures)
          for (s in 1:nsignatures) {
            if (any(Ediff[s,]<0)) FP_mut[s] <- sum(abs(Ediff[s,Ediff[s,]<0]))
            if (any(Ediff[s,]>0)) FN_mut[s] <- sum(abs(Ediff[s,Ediff[s,]>0]))
          }
          listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][["FP_mut"]][,pos] <- FP_mut/ngenomes
          listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][["FN_mut"]][,pos] <- FN_mut/ngenomes
          #need to know total number of times that a signature is not in E across the genomes
          sig_negatives <- apply(E==0,1,sum)
          listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][["FP_sig"]][,pos] <- apply(E==0 & E_hat>0,1,sum)/sig_negatives
          #need to know total number of times that a signature is in E across the genomes
          sig_positives <- apply(E>0,1,sum)
          listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][["FN_sig"]][,pos] <- apply(E>0 & E_hat==0,1,sum)/sig_positives
          
        }
        # for (alpha in alphas_KLD){
        #   E_hat <- SignatureFit(catalogue,signature_data_matrix = sigs[,cosmic_match],method = "KLD",bf_method = "KLD",alpha = alpha)
        #   Efit_RMSE_table_KLD[as.character(alpha),pos] <- sqrt(sum((E-E_hat)^2)/(ncol(E)*nrow(E)))
        # }
      }
      perf_file <- paste0(outFilePath_withFilterBest,"Sigs_OverallMetrics_SimulatedData_nboots",nboots_withFilter,".tsv")
      perf_table <- read.table(perf_file,sep = "\t",check.names = FALSE,header = TRUE,as.is = TRUE)
      if(cl=="MC"){
        ASW_table[,pos] <- perf_table$ave.SilWid.MC
      }else if(cl=="HC"){
        ASW_table[,pos] <- perf_table$ave.SilWid.hclust
      }else if(cl=="PAM"){
        ASW_table[,pos] <- perf_table$ave.SilWid.PAM
      }
      #get error just once
      if(cl==clustering_labels[1]){
        if(algo=="brunet" | algo=="nsNMF"){
          Error_table[,pos_error] <- perf_table$ave.KLD.orig
        }else if(algo=="lee"){
          Error_table[,pos_error] <- perf_table$ave.RMSE.orig
        }
        pos_error <- pos_error + 1
      }
      pos <- pos+1
    }
    
    #---Without Filter---
    for (i in 1:n){
      
      #read catalogues
      cat_file <- paste0(resultdir,"simulated_catalogue_",i,".tsv")
      E_file <- paste0(resultdir,"simulated_E_",i,".tsv")
      P_file <- paste0(resultdir,"simulated_P_",i,".tsv")
      if (file.exists(cat_file)){
        catalogue <- read.table(file = cat_file,sep = "\t",header = TRUE,check.names = FALSE,as.is = TRUE)
        E <- read.table(file = E_file,sep = "\t",header = TRUE,check.names = FALSE,as.is = TRUE)
        P <- read.table(file = P_file,sep = "\t",header = TRUE,check.names = FALSE,as.is = TRUE)
        message("catalogue ",i," loaded from file")
      } else {
        message("catalogue",i,"not found")
      }
      
      outFilePath_withoutFilterBest=paste0(resultdir,"without_filterBest_",algo,"_",i,"/")
      
      sig_file <- paste0(outFilePath_withoutFilterBest,"round_1/sig_",nsignatures,"/Sigs_plot_SimulatedData_ns",nsignatures,"_nboots",nboots_withoutFilter,"_",cl,".tsv")
      sigs <- read.table(sig_file,sep = "\t",check.names = FALSE,header = TRUE,as.is = TRUE)
      sigs.cor <- as.matrix(computeCorrelationOfTwoSetsOfSigs(sigs,P))
      res.match <- lp.assign(1 - sigs.cor)
      cosmic_match <- c()
      for(m in 1:nsignatures){
        cosmic_match <- c(cosmic_match,which.max(res.match$solution[,m]))
        signatures_similarity_table[m,pos] <- sigs.cor[m,which.max(res.match$solution[m,])]     
      }
      names(cosmic_match) <- colnames(sigs.cor)
      
      #sigFit and compute accuracy
      if(recompute_Efit_table | !file.exists(Efit_file)){
        for (alpha in alphas){
          message("threshold ",alpha)
          if (fit_method=="bleedingfilter"){
            E_hat <- SignatureFit(catalogue,signature_data_matrix = sigs[,cosmic_match],method = "KLD",bf_method = "CosSim",alpha = alpha)
          }else if(fit_method=="bootstrap"){
            file_store <- paste0(outFilePath_withoutFilterBest,"round_1/sig_",nsignatures,"/SigFit_withBootstrap_Summary_m","KLD","_bfm","CosSim","_alpha",-1,"_tr",alpha,"_p",p.value,"_",cl,".rData")
            if(file.exists(file_store)){
              load(file_store)
              message("Bootstrap Signature Fits loaded from file")
            }else{
              res <- SignatureFit_withBootstrap(catalogue,signature_data_matrix = sigs[,cosmic_match],method = "KLD",nboot = nboots_fit,threshold_percent = alpha,threshold_p.value = p.value,verbose = FALSE,nparallel = nparallel)
              save(file = file_store,res)
            }
            E_hat <- res$E_median_filtered
          }
          Efit_RMSE_table[as.character(alpha),pos] <- sqrt(sum((E-E_hat)^2)/(ncol(E)*nrow(E)))
          
          #classification metrics
          Ediff <- round(E)-E_hat
          #compute proportions relative to each sample number of mutations
          Ediff <- Ediff/matrix(rep(apply(E,2,sum),nsignatures),nrow = nsignatures,byrow = TRUE)
          #proportion of mutations per sample (i.e. divide by ngenomes)
          listOfClassificationMetrics_tables[["FP_mut"]][as.character(alpha),pos] <- sum(abs(Ediff[Ediff<0]))/ngenomes
          listOfClassificationMetrics_tables[["FN_mut"]][as.character(alpha),pos] <- sum(abs(Ediff[Ediff>0]))/ngenomes
          #proportion of mutations per sample (i.e. divide by ngenomes), proportion w.r.t actual positives and negatives in E
          listOfClassificationMetrics_tables[["FP_sig"]][as.character(alpha),pos] <- (sum(E==0 & E_hat>0)/ngenomes)/(nsignatures-nsignatures_per_sample)
          listOfClassificationMetrics_tables[["FN_sig"]][as.character(alpha),pos] <- (sum(E>0 & E_hat==0)/ngenomes)/(nsignatures_per_sample)
          
          #classification of individual signatures 
          FP_mut <- rep(0,nsignatures)
          FN_mut <- rep(0,nsignatures)
          for (s in 1:nsignatures) {
            if (any(Ediff[s,]<0)) FP_mut[s] <- sum(abs(Ediff[s,Ediff[s,]<0]))
            if (any(Ediff[s,]>0)) FN_mut[s] <- sum(abs(Ediff[s,Ediff[s,]>0]))
          }
          listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][["FP_mut"]][,pos] <- FP_mut/ngenomes
          listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][["FN_mut"]][,pos] <- FN_mut/ngenomes
          #need to know total number of times that a signature is not in E across the genomes
          sig_negatives <- apply(E==0,1,sum)
          listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][["FP_sig"]][,pos] <- apply(E==0 & E_hat>0,1,sum)/sig_negatives
          #need to know total number of times that a signature is in E across the genomes
          sig_positives <- apply(E>0,1,sum)
          listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][["FN_sig"]][,pos] <- apply(E>0 & E_hat==0,1,sum)/sig_positives
          
        }
        # for (alpha in alphas_KLD){
        #   E_hat <- SignatureFit(catalogue,signature_data_matrix = sigs[,cosmic_match],method = "KLD",bf_method = "KLD",alpha = alpha)
        #   Efit_RMSE_table_KLD[as.character(alpha),pos] <- sqrt(sum((E-E_hat)^2)/(ncol(E)*nrow(E)))
        # }
      }
      perf_file <- paste0(outFilePath_withoutFilterBest,"Sigs_OverallMetrics_SimulatedData_nboots",nboots_withoutFilter,".tsv")
      perf_table <- read.table(perf_file,sep = "\t",check.names = FALSE,header = TRUE,as.is = TRUE)
      if(cl=="MC"){
        ASW_table[,pos] <- perf_table$ave.SilWid.MC
      }else if(cl=="HC"){
        ASW_table[,pos] <- perf_table$ave.SilWid.hclust
      }else if(cl=="PAM"){
        ASW_table[,pos] <- perf_table$ave.SilWid.PAM
      }
      #get error just once
      if(cl==clustering_labels[1]){
        if(algo=="brunet" | algo=="nsNMF"){
          Error_table[,pos_error] <- perf_table$ave.KLD.orig
        }else if(algo=="lee"){
          Error_table[,pos_error] <- perf_table$ave.RMSE.orig
        }
        pos_error <- pos_error + 1
      }
      pos <- pos+1
    }
  }
}

if(recompute_Efit_table | !file.exists(Efit_file)){
  write.table(Efit_RMSE_table,file = Efit_file,quote = FALSE,sep = "\t",col.names = TRUE,row.names = TRUE)
  for (m in listOfClassificationMetrics){
    write.table(listOfClassificationMetrics_tables[[m]],
                file = listOfClassificationMetrics_files[[m]],
                quote = FALSE,sep = "\t",col.names = TRUE,row.names = TRUE)
  }
  for(alpha in alphas){
    for (m in listOfClassificationMetrics){
      write.table(listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][[m]],
                  file = listOfClassificationMetrics_IndSig_files[[as.character(alpha)]][[m]],
                  quote = FALSE,sep = "\t",col.names = TRUE,row.names = TRUE)
    }
  }
}else{
  Efit_RMSE_table <- read.table(Efit_file,sep = "\t",stringsAsFactors = FALSE,header = TRUE,as.is = TRUE,check.names = FALSE)
  listOfClassificationMetrics_tables <- list()
  for (m in listOfClassificationMetrics){
    listOfClassificationMetrics_tables[[m]] <- read.table(listOfClassificationMetrics_files[[m]],
                                                          sep = "\t",stringsAsFactors = FALSE,
                                                          header = TRUE,as.is = TRUE,check.names = FALSE)
  }
  for(alpha in alphas){
    for (m in listOfClassificationMetrics){
      listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][[m]] <- read.table(listOfClassificationMetrics_IndSig_files[[as.character(alpha)]][[m]],
                                                            sep = "\t",stringsAsFactors = FALSE,
                                                            header = TRUE,as.is = TRUE,check.names = FALSE)
    }
  }
}

#stats mean and standard error of mean
all_labels <- paste(rep(nmf_algorithms_labels,each = nsettings),
                    paste(rep(filter_labels,length(nmf_algorithms)*length(clustering_labels)),
                          rep(rep(clustering_labels,each = length(filter_labels)),length(nmf_algorithms)),sep = " "),
                    sep=" ")
err_labels <-  paste(rep(nmf_algorithms_labels,each = length(filter_labels)),
                     paste(rep(filter_labels,length(nmf_algorithms)),sep = " "),
                     sep=" ")

#ASW mean
ASW_table_mean <- matrix(0,nrow = length(nsig_range), ncol = length(nmf_algorithms)*nsettings)
colnames(ASW_table_mean) <- all_labels
row.names(ASW_table_mean) <- nsig_range
for (i in 1:n){
  ASW_table_mean <- ASW_table_mean + ASW_table[,seq(i,ncol(ASW_table),n)]
}
ASW_table_mean <- ASW_table_mean/n
#ASW standard error
ASW_table_se <- matrix(0,nrow = length(nsig_range), ncol = length(nmf_algorithms)*nsettings)
colnames(ASW_table_se) <- all_labels
row.names(ASW_table_se) <- nsig_range
for (i in 1:n){
  ASW_table_se <- ASW_table_se + (ASW_table[,seq(i,ncol(ASW_table),n)] - ASW_table_mean)^2
}
ASW_table_se <- sqrt(ASW_table_se/n)/sqrt(n)


#Efit mean
Efit_table_mean <- matrix(0,nrow = length(alphas), ncol = length(nmf_algorithms)*nsettings)
colnames(Efit_table_mean) <- all_labels
row.names(Efit_table_mean) <- alphas
for (i in 1:n){
  Efit_table_mean <- Efit_table_mean + Efit_RMSE_table[,seq(i,ncol(Efit_RMSE_table),n)]
}
colnames(Efit_table_mean) <- all_labels
Efit_table_mean <- Efit_table_mean/n
#Efit standard error
Efit_table_se <- matrix(0,nrow = length(alphas), ncol = length(nmf_algorithms)*nsettings)
colnames(Efit_table_se) <- all_labels
row.names(Efit_table_se) <- alphas
for (i in 1:n){
  Efit_table_se <- Efit_table_se + (Efit_RMSE_table[,seq(i,ncol(Efit_RMSE_table),n)] - Efit_table_mean)^2
}
Efit_table_se <- sqrt(Efit_table_se/n)/sqrt(n)

#Error mean
Error_table_mean <- matrix(0,nrow = length(nsig_range), ncol = length(nmf_algorithms)*length(filter_labels))
colnames(Error_table_mean) <- err_labels
row.names(Error_table_mean) <- nsig_range
for (i in 1:n){
  Error_table_mean <- Error_table_mean + Error_table[,seq(i,ncol(Error_table),n)]
}
Error_table_mean <- Error_table_mean/n
#Error standard error
Error_table_se <- matrix(0,nrow = length(nsig_range), ncol = length(nmf_algorithms)*length(filter_labels))
colnames(Error_table_se) <- err_labels
row.names(Error_table_se) <- nsig_range
for (i in 1:n){
  Error_table_se <- Error_table_se + (Error_table[,seq(i,ncol(Error_table),n)] - Error_table_mean)^2
}
Error_table_se <- sqrt(Error_table_se/n)/sqrt(n)
#normalise error tables
for (i in 1:length(nmf_algorithms)){
  positions <- 1:length(filter_labels) + (i-1)*length(filter_labels)
  norm_matrix <- matrix(rep(max(Error_table_mean[,positions]),nrow(Error_table_mean[,positions])*length(filter_labels)),nrow = nrow(Error_table_mean),byrow = TRUE)
  Error_table_mean[,positions] <- Error_table_mean[,positions]/norm_matrix
  Error_table_se[,positions] <- Error_table_se[,positions]/norm_matrix
}


#Classification
listOfClassificationMetrics_tables_mean <- list()
listOfClassificationMetrics_tables_se <- list()
for (m in listOfClassificationMetrics){
  #Classification mean
  listOfClassificationMetrics_tables_mean[[m]] <- matrix(0,nrow = length(alphas), ncol = length(nmf_algorithms)*nsettings)
  colnames(listOfClassificationMetrics_tables_mean[[m]]) <- all_labels
  row.names(listOfClassificationMetrics_tables_mean[[m]]) <- alphas
  for (i in 1:n){
    listOfClassificationMetrics_tables_mean[[m]] <- listOfClassificationMetrics_tables_mean[[m]] + listOfClassificationMetrics_tables[[m]][,seq(i,ncol(listOfClassificationMetrics_tables[[m]]),n)]
  }
  colnames(listOfClassificationMetrics_tables_mean[[m]]) <- all_labels
  listOfClassificationMetrics_tables_mean[[m]] <- listOfClassificationMetrics_tables_mean[[m]]/n
  #classification standard error
  listOfClassificationMetrics_tables_se[[m]] <- matrix(0,nrow = length(alphas), ncol = length(nmf_algorithms)*nsettings)
  colnames(listOfClassificationMetrics_tables_se[[m]]) <- all_labels
  row.names(listOfClassificationMetrics_tables_se[[m]]) <- alphas
  for (i in 1:n){
    listOfClassificationMetrics_tables_se[[m]] <- listOfClassificationMetrics_tables_se[[m]] + (listOfClassificationMetrics_tables[[m]][,seq(i,ncol(listOfClassificationMetrics_tables[[m]]),n)] - listOfClassificationMetrics_tables_mean[[m]])^2
  }
  listOfClassificationMetrics_tables_se[[m]] <- sqrt(listOfClassificationMetrics_tables_se[[m]]/n)/sqrt(n)
}

#Classification of mutations and sigantures, individual signatures
listOfClassificationMetrics_IndSig_tables_mean <- list()
listOfClassificationMetrics_IndSig_tables_se <- list()
for(alpha in alphas){
  listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]] <- list()
  listOfClassificationMetrics_IndSig_tables_se[[as.character(alpha)]] <- list()
  for (m in listOfClassificationMetrics){
    listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]] <- matrix(0,nrow = length(1:nsignatures), ncol = length(nmf_algorithms)*nsettings)
    colnames(listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]]) <- all_labels
    row.names(listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]]) <- colnames(P)
    for (i in 1:n){
      listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]] <- listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]] + listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][[m]][,seq(i,ncol(listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][[m]]),n)]
    }
    colnames(listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]]) <- all_labels
    row.names(listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]]) <- colnames(P)
    listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]] <- listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]]/n
    #classification standard error
    listOfClassificationMetrics_IndSig_tables_se[[as.character(alpha)]][[m]] <- matrix(0,nrow = length(1:nsignatures), ncol = length(nmf_algorithms)*nsettings)
    colnames(listOfClassificationMetrics_IndSig_tables_se[[as.character(alpha)]][[m]]) <- all_labels
    row.names(listOfClassificationMetrics_IndSig_tables_se[[as.character(alpha)]][[m]]) <- colnames(P)
    for (i in 1:n){
      listOfClassificationMetrics_IndSig_tables_se[[as.character(alpha)]][[m]] <- listOfClassificationMetrics_IndSig_tables_se[[as.character(alpha)]][[m]] + (listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][[m]][,seq(i,ncol(listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][[m]]),n)] - listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]])^2
    }
    listOfClassificationMetrics_IndSig_tables_se[[as.character(alpha)]][[m]] <- sqrt(listOfClassificationMetrics_IndSig_tables_se[[as.character(alpha)]][[m]]/n)/sqrt(n)
  }
}

#----plots----

colours_list <- c("red",
                  "green",
                  "skyblue",
                  "purple",
                  "orange",
                  "grey")


file_out <- paste0(resultdir,"2018_03_12_CosineSimilarity_performance.jpg")
jpeg(filename = file_out,width = 2000,height = 800,res = 200)
xlabel <- paste(rep(nmf_algorithms_labels,each = nsettings),
                paste(rep(filter_labels,length(nmf_algorithms)*length(clustering_plotlabels)),
                rep(rep(clustering_plotlabels,each = length(filter_labels)),length(nmf_algorithms)),sep = " "),
                sep="\n")
par(mar=c(6,4,3,7))
#start from replicate 1
i <- 1
selectedCols <-  seq(i,n*nsettings*length(nmf_algorithms),n)
plot(rep(1:(nsettings*length(nmf_algorithms)),each=nsignatures),
     signatures_similarity_table[,selectedCols],
     type = "p",
     pch = 21,
     bg = "red",
     ylim = c(min(signatures_similarity_table),1),
     xlim = c(0.5,nsettings*length(nmf_algorithms)+0.5),
     ylab = "Cosine Similarity",
     xlab = "",
     main = "Similarity of retrieved signature to COSMIC signatures",
     xaxt = "n",
     panel.first = abline(v=1:(nsettings*length(nmf_algorithms)), h=seq(0.1,1,0.1),lty=3, col="lightgrey"))
axis(1, at=1:(nsettings*length(nmf_algorithms)), labels=xlabel,las=2,cex.axis=0.8)
if(n>1){
  for (i in 2:n) {
    selectedCols <-  seq(i,n*nsettings*length(nmf_algorithms),n)
    points(rep(1:(nsettings*length(nmf_algorithms)),each=nsignatures) + ifelse((i-1) %% 2,(i%/%2)*0.1,-(i%/%2)*0.1),
           signatures_similarity_table[,selectedCols],
           type = "p",pch = 21,bg = colours_list[i])
  }
}
legend("right",legend = paste0("replicate ",1:n),
       pch = 21,pt.bg = colours_list[1:n],
       xpd = TRUE,bty = "n",inset = c(-0.15,-0.4))
dev.off()


#plot ASW
current_labels <- paste(paste(rep(filter_labels,length(clustering_plotlabels)),
                              rep(rep(clustering_plotlabels,each = length(filter_labels)),1),sep = " "),
                        sep=" ")
file_out <- paste0(resultdir,"2018_03_12_ASW_performance.jpg")
jpeg(filename = file_out,width = 550*length(nmf_algorithms),height = 700,res = 200)
par(mfrow=c(1,length(nmf_algorithms)),oma=c(2,1,1.5,1))
par(mar=c(5,4,4,1))
for (j in 1:length(nmf_algorithms)){
  plot(nsig_range, ASW_table_mean[,1+(j-1)*nsettings],
       ylim=c(min(ASW_table_mean - ASW_table_se), max(1,max(ASW_table_mean + ASW_table_se))),
       t="l", xlab="n signatures", ylab="Average Silhouette Width",
       main=nmf_algorithms_labels[j],
       lwd=2,col=colours_list[1],panel.first=grid(),cex.lab=1.3,cex.axis=1.3)
  # hack: we draw arrows but with very special "arrowheads"
  if (n>1) arrows(nsig_range, ASW_table_mean[,1+(j-1)*nsettings]-ASW_table_se[,1+(j-1)*nsettings], nsig_range, ASW_table_mean[,1+(j-1)*nsettings]+ASW_table_se[,1+(j-1)*nsettings], length=0.05, angle=90, code=3,lwd=2,col=colours_list[1])
  
  for(i in 2:nsettings){
    lines(nsig_range, ASW_table_mean[,i+(j-1)*nsettings],lwd=2,col=colours_list[i])
    if (n>1) arrows(nsig_range, ASW_table_mean[,i+(j-1)*nsettings]-ASW_table_se[,i+(j-1)*nsettings], nsig_range, ASW_table_mean[,i+(j-1)*nsettings]+ASW_table_se[,i+(j-1)*nsettings], length=0.05, angle=90, code=3,lwd=2,col=colours_list[i])
    
  }

}
title(paste0("Average Silhoutte Width (Mean and SE, n=",n,")"),outer = TRUE)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom",legend = current_labels,horiz = TRUE,
       lwd = 2,col = colours_list[1:nsettings],
       xpd = TRUE,bty = "n",inset = c(0,0),cex = 1.2)
dev.off()


#plot Efit
current_labels <- paste(paste(rep(filter_labels,length(clustering_plotlabels)),
                          rep(rep(clustering_plotlabels,each = length(filter_labels)),1),sep = " "),
                    sep=" ")
if (fit_method=="bleedingfilter"){
  file_out <- paste0(resultdir,"2018_03_27_Efit_performance.jpg")
  xlabel <- "alpha"
}else if(fit_method=="bootstrap"){
  file_out <- paste0(resultdir,"2018_03_27_Efit_performance_bootstraps.jpg")
  xlabel <- "threshold (%)"
}
jpeg(filename = file_out,width = 550*length(nmf_algorithms),height = 700,res = 200)
par(mfrow=c(1,length(nmf_algorithms)),oma=c(2,1,1.5,1))
par(mar=c(5,4,4,1))
for (j in 1:length(nmf_algorithms)){
  plot(1:length(alphas), Efit_table_mean[,1+(j-1)*nsettings],
       ylim=c(min(Efit_table_mean - Efit_table_se), max(1,max(Efit_table_mean + Efit_table_se))),
       t="l", xlab=xlabel, ylab="RMSE",
       main=nmf_algorithms_labels[j],
       lwd=2,col=colours_list[1],panel.first=grid(),xaxt='n',cex.lab=1.3,cex.axis=1.3)
  axis(1, at=1:length(alphas), labels=alphas,cex.axis=1.3)
  # hack: we draw arrows but with very special "arrowheads"
  if (n>1) arrows(1:length(alphas), Efit_table_mean[,1+(j-1)*nsettings]-Efit_table_se[,1+(j-1)*nsettings], 1:length(alphas), Efit_table_mean[,1+(j-1)*nsettings]+Efit_table_se[,1+(j-1)*nsettings], length=0.05, angle=90, code=3,lwd=2,col=colours_list[1])
  
  for(i in 2:nsettings){
    lines(1:length(alphas), Efit_table_mean[,i+(j-1)*nsettings],lwd=2,col=colours_list[i])
    if (n>1) arrows(1:length(alphas), Efit_table_mean[,i+(j-1)*nsettings]-Efit_table_se[,i+(j-1)*nsettings], 1:length(alphas), Efit_table_mean[,i+(j-1)*nsettings]+Efit_table_se[,i+(j-1)*nsettings], length=0.05, angle=90, code=3,lwd=2,col=colours_list[i])
    
  }

}
title(paste0("RMSE of Activity Matrix (Mean and SE, n=",n,")"),outer = TRUE)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom",legend = current_labels,horiz = TRUE,
       lwd = 2,col = colours_list[1:nsettings],
       xpd = TRUE,bty = "n",inset = c(0,0),cex = 1.2)
dev.off()

#plot error
file_out <- paste0(resultdir,"2018_03_12_error_performance.jpg")
jpeg(filename = file_out,width = 700,height = 700,res = 150)
par(mar=c(5,4,4,2))
plot(nsig_range, Error_table_mean[,1],
     ylim=c(min(Error_table_mean - Error_table_se),1.1*max(Error_table_mean + Error_table_se)),
     t="l", xlab="n signatures", ylab="Normalised Mean Error",
     main=paste0("Error (Mean and Standard Error, n=",n,")"),
     lwd=2,col=colours_list[1],panel.first=grid(),cex.lab=1.3,cex.axis=1.3)
# hack: we draw arrows but with very special "arrowheads"
if (n>1) arrows(nsig_range, Error_table_mean[,1]-Error_table_se[,1], nsig_range, Error_table_mean[,1]+Error_table_se[,1], length=0.05, angle=90, code=3,lwd=2,col=colours_list[1])

for(i in 2:length(filter_labels)){
  lines(nsig_range, Error_table_mean[,i],lwd=2,col=colours_list[i])
  if (n>1) arrows(nsig_range, Error_table_mean[,i]-Error_table_se[,i], nsig_range, Error_table_mean[,i]+Error_table_se[,i], length=0.05, angle=90, code=3,lwd=2,col=colours_list[i])
  
}
if(length(nmf_algorithms)>1){
  for (j in 2:length(nmf_algorithms)){
    for(i in 1:length(filter_labels)){
      lines(nsig_range, Error_table_mean[,i+(j-1)*length(filter_labels)],lwd=2,col=colours_list[i+(j-1)*length(filter_labels)])
      if (n>1) arrows(nsig_range, Error_table_mean[,i+(j-1)*length(filter_labels)]-Error_table_se[,i+(j-1)*length(filter_labels)], nsig_range, Error_table_mean[,i+(j-1)*length(filter_labels)]+Error_table_se[,i+(j-1)*length(filter_labels)], length=0.05, angle=90, code=3,lwd=2,col=colours_list[i+(j-1)*length(filter_labels)])
      
    }
  }
}
legend("topright",legend = colnames(Error_table_mean),
       lwd = 2,col = colours_list[1:ncol(Error_table_mean)],
       xpd = TRUE,bty = "n",inset = c(0,0),cex = 1)
dev.off()

#plot Classification Metrics
current_labels <- paste(paste(rep(filter_labels,length(clustering_plotlabels)),
                              rep(rep(clustering_plotlabels,each = length(filter_labels)),1),sep = " "),
                        sep=" ")
ml <-1
for(m in listOfClassificationMetrics){
  if (fit_method=="bleedingfilter"){
    file_out <- paste0(resultdir,"2018_03_27_E_metrics_performance",m,".jpg")
    xlabel <- "alpha"
  }else if(fit_method=="bootstrap"){
    file_out <- paste0(resultdir,"2018_03_27_E_metrics_performance",m,"_bootstraps.jpg")
    xlabel <- "threshold (%)"
  }
  jpeg(filename = file_out,width = 550*length(nmf_algorithms),height = 700,res = 200)
  par(mfrow=c(1,length(nmf_algorithms)),oma=c(2,1,1.5,1))
  par(mar=c(5,4,4,1))
  for (j in 1:length(nmf_algorithms)){
    plot(1:length(alphas), listOfClassificationMetrics_tables_mean[[m]][,1+(j-1)*nsettings],
         ylim=c(min(listOfClassificationMetrics_tables_mean[[m]] - listOfClassificationMetrics_tables_se[[m]]), max(0.5,max(listOfClassificationMetrics_tables_mean[[m]] + listOfClassificationMetrics_tables_se[[m]]))),
         t="l", xlab=xlabel, ylab=paste0(listOfClassificationMetrics_labels[[ml]]),
         main=nmf_algorithms_labels[j],
         lwd=2,col=colours_list[1],panel.first=grid(),xaxt='n',cex.lab=1.3,cex.axis=1.3)
    axis(1, at=1:length(alphas), labels=alphas,cex.axis=1.3)
    # hack: we draw arrows but with very special "arrowheads"
    if (n>1) arrows(1:length(alphas), listOfClassificationMetrics_tables_mean[[m]][,1+(j-1)*nsettings]-listOfClassificationMetrics_tables_se[[m]][,1+(j-1)*nsettings], 1:length(alphas), listOfClassificationMetrics_tables_mean[[m]][,1+(j-1)*nsettings]+listOfClassificationMetrics_tables_se[[m]][,1+(j-1)*nsettings], length=0.05, angle=90, code=3,lwd=2,col=colours_list[1])
    
    for(i in 2:nsettings){
      lines(1:length(alphas), listOfClassificationMetrics_tables_mean[[m]][,i+(j-1)*nsettings],lwd=2,col=colours_list[i])
      if (n>1) arrows(1:length(alphas), listOfClassificationMetrics_tables_mean[[m]][,i+(j-1)*nsettings]-listOfClassificationMetrics_tables_se[[m]][,i+(j-1)*nsettings], 1:length(alphas), listOfClassificationMetrics_tables_mean[[m]][,i+(j-1)*nsettings]+listOfClassificationMetrics_tables_se[[m]][,i+(j-1)*nsettings], length=0.05, angle=90, code=3,lwd=2,col=colours_list[i])
      
    }
    
  }
  title(paste0("Mutation Assignment ",listOfClassificationMetrics_labels[[ml]]," (Mean and SE, n=",n,")"),outer = TRUE)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom",legend = current_labels,horiz = TRUE,
         lwd = 2,col = colours_list[1:nsettings],
         xpd = TRUE,bty = "n",inset = c(0,0),cex = 1.2)
  dev.off()
  ml <- ml+1
}


#plot 1 - Classification Metrics (Sensitivity/Specificity)
current_labels <- paste(paste(rep(filter_labels,length(clustering_plotlabels)),
                              rep(rep(clustering_plotlabels,each = length(filter_labels)),1),sep = " "),
                        sep=" ")
ml <-1
for(m in listOfClassificationMetrics){
  mss <- listOfClassificationMetricsSensSpec[ml]
  if (fit_method=="bleedingfilter"){
    file_out <- paste0(resultdir,"2018_03_27_E_metrics_performance",mss,".jpg")
    xlabel <- "alpha"
  }else if(fit_method=="bootstrap"){
    file_out <- paste0(resultdir,"2018_03_27_E_metrics_performance",mss,"_bootstraps.jpg")
    xlabel <- "threshold (%)"
  }
  jpeg(filename = file_out,width = 550*length(nmf_algorithms),height = 700,res = 200)
  par(mfrow=c(1,length(nmf_algorithms)),oma=c(2,1,1.5,1))
  par(mar=c(5,4,4,1))
  for (j in 1:length(nmf_algorithms)){
    plot(1:length(alphas), 1 - listOfClassificationMetrics_tables_mean[[m]][,1+(j-1)*nsettings],
         ylim=c(min(1 - (listOfClassificationMetrics_tables_mean[[m]] - listOfClassificationMetrics_tables_se[[m]])), max(0.01,max(1- (listOfClassificationMetrics_tables_mean[[m]] + listOfClassificationMetrics_tables_se[[m]])))),
         t="l", xlab=xlabel, ylab=paste0(listOfClassificationMetricsSensSpec_labels[ml]),
         main=nmf_algorithms_labels[j],
         lwd=2,col=colours_list[1],panel.first=grid(),xaxt='n',cex.lab=1.3,cex.axis=1.3)
    axis(1, at=1:length(alphas), labels=alphas,cex.axis=1.3)
    # hack: we draw arrows but with very special "arrowheads"
    if (n>1) arrows(1:length(alphas), 1 - (listOfClassificationMetrics_tables_mean[[m]][,1+(j-1)*nsettings]-listOfClassificationMetrics_tables_se[[m]][,1+(j-1)*nsettings]), 1:length(alphas), 1 - (listOfClassificationMetrics_tables_mean[[m]][,1+(j-1)*nsettings]+listOfClassificationMetrics_tables_se[[m]][,1+(j-1)*nsettings]), length=0.05, angle=90, code=3,lwd=2,col=colours_list[1])
    
    for(i in 2:nsettings){
      lines(1:length(alphas), 1 - listOfClassificationMetrics_tables_mean[[m]][,i+(j-1)*nsettings],lwd=2,col=colours_list[i])
      if (n>1) arrows(1:length(alphas), 1 - (listOfClassificationMetrics_tables_mean[[m]][,i+(j-1)*nsettings]-listOfClassificationMetrics_tables_se[[m]][,i+(j-1)*nsettings]), 1:length(alphas), 1 - (listOfClassificationMetrics_tables_mean[[m]][,i+(j-1)*nsettings]+listOfClassificationMetrics_tables_se[[m]][,i+(j-1)*nsettings]), length=0.05, angle=90, code=3,lwd=2,col=colours_list[i])
      
    }
    
  }
  title(paste0("Mutation Assignment ",listOfClassificationMetricsSensSpec_labels[[ml]]," (Mean and SE, n=",n,")"),outer = TRUE)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom",legend = current_labels,horiz = TRUE,
         lwd = 2,col = colours_list[1:nsettings],
         xpd = TRUE,bty = "n",inset = c(0,0),cex = 1.2)
  dev.off()
  ml <- ml+1
}

#plot Classification Metrics for each Signature
colours_list <- c("red",
                  "green",
                  "skyblue",
                  "purple",
                  "orange",
                  "grey")
current_labels <- paste(paste(rep(filter_labels,length(clustering_plotlabels)),
                              rep(rep(clustering_plotlabels,each = length(filter_labels)),1),sep = " "),
                        sep=" ")
sig_lab <- paste0("Sig.",substr(colnames(P),11,12))
ml <-1
for(m in listOfClassificationMetrics){
  if (fit_method=="bleedingfilter"){
    file_out <- paste0(resultdir,"2018_03_27_E_metrics_performance_Signatures_",m,".jpg")
    xlabel <- "alpha"
  }else if(fit_method=="bootstrap"){
    file_out <- paste0(resultdir,"2018_03_27_E_metrics_performance_Signatures_",m,"_bootstraps.jpg")
    xlabel <- "threshold (%)"
  }
  jpeg(filename = file_out,width = 1000*length(nmf_algorithms),height = 600*length(alphas),res = 300)
  par(mar=c(5,4,4,0),mgp=c(2,0.5,0),mfrow=c(length(alphas),length(nmf_algorithms)),oma=c(2,1,1.5,1))
  for(alpha in alphas){
    for (j in 1:length(nmf_algorithms)){
      bb <- barplot(t(as.matrix(listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]][,1:nsettings+(j-1)*nsettings])),
              beside = TRUE,las = 2,col = colours_list,names.arg = sig_lab,cex.names = 0.9,
              ylim=c(min(0,min(listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]] - listOfClassificationMetrics_IndSig_tables_se[[as.character(alpha)]][[m]])), max(1,max(listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]] + listOfClassificationMetrics_IndSig_tables_se[[as.character(alpha)]][[m]]))),
              main = paste0(nmf_algorithms_labels[j],", ",xlabel,"=",alpha),ylab = paste0(listOfClassificationMetrics_labels[[ml]]))
      if (n>1) arrows(unlist(t(bb)), unlist(listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]][,1:nsettings+(j-1)*nsettings]-listOfClassificationMetrics_IndSig_tables_se[[as.character(alpha)]][[m]][,1:nsettings+(j-1)*nsettings]), 
                      unlist(t(bb)), unlist(listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]][,1:nsettings+(j-1)*nsettings]+listOfClassificationMetrics_IndSig_tables_se[[as.character(alpha)]][[m]][,1:nsettings+(j-1)*nsettings]), length=0.015, angle=90, code=3,lwd=1)
    }
  }
  title(paste0("Mutation Assignment ",listOfClassificationMetrics_labels[[ml]]," (Mean and SE, n=",n,")"),outer = TRUE)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom",legend = current_labels,horiz = TRUE,
         fill = colours_list[1:nsettings],
         xpd = TRUE,bty = "n",inset = c(0,0),cex = 1)
  dev.off()
  ml <- ml+1
}

#plot 1 - Classification Metrics for each Signature (Sensitivity/Specificity)
colours_list <- c("red",
                  "green",
                  "skyblue",
                  "purple",
                  "orange",
                  "grey")
current_labels <- paste(paste(rep(filter_labels,length(clustering_plotlabels)),
                              rep(rep(clustering_plotlabels,each = length(filter_labels)),1),sep = " "),
                        sep=" ")
sig_lab <- paste0("Sig.",substr(colnames(P),11,12))
ml <-1
for(m in listOfClassificationMetrics){
  mss <- listOfClassificationMetricsSensSpec[ml]
  if (fit_method=="bleedingfilter"){
    file_out <- paste0(resultdir,"2018_03_27_E_metrics_performance_Signatures_",mss,".jpg")
    xlabel <- "alpha"
  }else if(fit_method=="bootstrap"){
    file_out <- paste0(resultdir,"2018_03_27_E_metrics_performance_Signatures_",mss,"_bootstraps.jpg")
    xlabel <- "threshold (%)"
  }
  jpeg(filename = file_out,width = 1000*length(nmf_algorithms),height = 600*length(alphas),res = 300)
  par(mar=c(5,4,4,0),mgp=c(2,0.5,0),mfrow=c(length(alphas),length(nmf_algorithms)),oma=c(2,1,1.5,1))
  for(alpha in alphas){
    for (j in 1:length(nmf_algorithms)){
      bb <- barplot(1 - t(as.matrix(listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]][,1:nsettings+(j-1)*nsettings])),
                    beside = TRUE,las = 2,col = colours_list,names.arg = sig_lab,cex.names = 0.9,
                    ylim=c(min(0,min(1 - (listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]] - listOfClassificationMetrics_IndSig_tables_se[[as.character(alpha)]][[m]]))), max(1,max(1 - (listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]] + listOfClassificationMetrics_IndSig_tables_se[[as.character(alpha)]][[m]])))),
                    main = paste0(nmf_algorithms_labels[j],", ",xlabel,"=",alpha),ylab = paste0(listOfClassificationMetricsSensSpec_labels[[ml]]))
      if (n>1) arrows(unlist(t(bb)), 1 - unlist(listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]][,1:nsettings+(j-1)*nsettings]-listOfClassificationMetrics_IndSig_tables_se[[as.character(alpha)]][[m]][,1:nsettings+(j-1)*nsettings]), 
                      unlist(t(bb)), 1 - unlist(listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]][,1:nsettings+(j-1)*nsettings]+listOfClassificationMetrics_IndSig_tables_se[[as.character(alpha)]][[m]][,1:nsettings+(j-1)*nsettings]), length=0.015, angle=90, code=3,lwd=1)
    }
  }
  title(paste0("Mutation Assignment ",listOfClassificationMetricsSensSpec_labels[[ml]]," (Mean and SE, n=",n,")"),outer = TRUE)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom",legend = current_labels,horiz = TRUE,
         fill = colours_list[1:nsettings],
         xpd = TRUE,bty = "n",inset = c(0,0),cex = 1)
  dev.off()
  ml <- ml+1
}

