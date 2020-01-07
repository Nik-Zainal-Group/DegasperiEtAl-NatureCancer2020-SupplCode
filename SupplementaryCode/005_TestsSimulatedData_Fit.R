source("../lib/SignatureExtractionLib.R")
source("../lib/SignatureFitLib.R")
#library("lpSolve")

recompute_Efit_table <- FALSE
nparallel <- 5

#Cosmic Signatures
cosmic30 <- read.table("../data/COSMIC_signatures.txt",sep = "\t",header = TRUE,as.is = TRUE,check.names = FALSE)
cosmic30 <- sortCatalogue(cosmic30)
signatures <- c(1,2,3,5,6,8,12,13,17,18)
P <- cosmic30[,signatures]

#parameters
nsignatures <- ncol(P)

n <- 5
ngenomes <- 30
nsignatures_per_sample <- 5
resultdir <- paste0("../results/testsSimulatedData/sigPerSample",nsignatures_per_sample,"_samples",ngenomes,"/")
# alpha <- -1
# threshold_percent <- 1
# threshold_p.value <- 0.05
# nboot <- 20
# method <- "KLD"
# bf_method <- "CosSim"

fit_methods <- c("Bootstraps KLD","Bootstraps NNLS","KLD","NNLS")
#fit_methods <- fit_methods[1]

#Signature Fit with Bootstraps
thresholds_perc <- c(0,1,2,5,10)
nboots_fit <- 100
p.value <- 0.05

all_labels <- paste(rep(fit_methods,each = n),
                    paste0("rep",rep(1:n,length(fit_methods))),
                    sep=" ")

outdir <- paste0(resultdir,"sigfit_tests/")
dir.create(outdir,showWarnings = FALSE,recursive = TRUE)

Efit_RMSE_table <- matrix(NA,nrow = length(thresholds_perc), ncol = length(fit_methods)*n)
row.names(Efit_RMSE_table) <- as.character(thresholds_perc)
colnames(Efit_RMSE_table) <- all_labels
Efit_file <- paste0(outdir,"ExposureFit_bootstraps.tsv")


#Classification of mutations and sigantures
listOfClassificationMetrics <- c("FP_mut","FP_sig","FN_mut","FN_sig")
listOfClassificationMetrics_labels <- c("FP rate of mutation assignment","FP rate of signature assignment","FN rate of mutation assignment","FN rate of signature assignment")
listOfClassificationMetricsSensSpec <- c("Specificity_mut","Specificity_sig","Sensitivity_mut","Sensitivity_sig")
listOfClassificationMetricsSensSpec_labels <- c("Specificity of mutation assignment","Specificity of signature assignment","Sensitivity of mutation assignment","Sensitivity of signature assignment")
listOfClassificationMetrics_tables <- list()
listOfClassificationMetrics_files <- list()


for (m in listOfClassificationMetrics){
  listOfClassificationMetrics_tables[[m]] <- matrix(NA,nrow = length(thresholds_perc), ncol = length(fit_methods)*n)
  row.names(listOfClassificationMetrics_tables[[m]]) <- as.character(thresholds_perc)
  colnames(listOfClassificationMetrics_tables[[m]]) <- all_labels
  listOfClassificationMetrics_files[[m]] <- paste0(outdir,"ExposureClassification_",m,"_bootstraps.tsv")
}


#Classification of mutations and signatures, individual signatures
listOfClassificationMetrics_IndSig_tables <- list()
listOfClassificationMetrics_IndSig_files <- list()
for(alpha in thresholds_perc){
  listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]] <- list()
  listOfClassificationMetrics_IndSig_files[[as.character(alpha)]] <- list()
  for (m in listOfClassificationMetrics){
    listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][[m]] <- matrix(NA,nrow = length(1:nsignatures), ncol = length(fit_methods)*n)
    row.names(listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][[m]]) <- as.character(1:nsignatures)
    colnames(listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][[m]]) <- all_labels
    listOfClassificationMetrics_IndSig_files[[as.character(alpha)]][[m]] <- paste0(outdir,"ExposureClassification_IndSig_",m,"_a",alpha,"_bootstraps.tsv")
  }
}

pos <- 1
for (im in 1:length(fit_methods)){
  message("method ",fit_methods[im])
  for (i in 1:n){
    iDir <- paste0(outdir,"rep_",i,"/")
    dir.create(iDir,showWarnings = FALSE,recursive = TRUE)
    
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
    
    # E.on.cat <- SignatureFit(cat = as.matrix(catalogue),
    #              signature_data_matrix = as.matrix(P),
    #              method = "KLD",
    #              bf_method = "CosSim",
    #              alpha = alpha)
    
    #sigFit and compute accuracy
    if(recompute_Efit_table | !file.exists(Efit_file)){
      for (alpha in thresholds_perc){
        message("threshold ",alpha)
        
        if("Bootstraps KLD"==fit_methods[im]){
          file_store <- paste0(iDir,"SigFit_withBootstrap_Summary_m","KLD","_bfm","CosSim","_alpha",-1,"_tr",alpha,"_p",p.value,".rData")
          if(file.exists(file_store)){
            load(file_store)
            message("Bootstrap Signature Fits loaded from file")
          }else{
            res <- SignatureFit_withBootstrap(catalogue,signature_data_matrix = P,method = "KLD",nboot = nboots_fit,threshold_percent = alpha,threshold_p.value = p.value,verbose = FALSE,nparallel = nparallel)
            save(file = file_store,res)
          }
          E_hat <- res$E_median_filtered
        }else if("Bootstraps NNLS"==fit_methods[im]){
          file_store <- paste0(iDir,"SigFit_withBootstrap_Summary_m","NNLS","_bfm","CosSim","_alpha",-1,"_tr",alpha,"_p",p.value,".rData")
          if(file.exists(file_store)){
            load(file_store)
            message("Bootstrap Signature Fits loaded from file")
          }else{
            res <- SignatureFit_withBootstrap(catalogue,signature_data_matrix = P,method = "NNLS",nboot = nboots_fit,threshold_percent = alpha,threshold_p.value = p.value,verbose = FALSE,nparallel = nparallel)
            save(file = file_store,res)
          }
          E_hat <- res$E_median_filtered
        }else if(fit_methods[im]=="KLD"){
          E_hat <- SignatureFit(catalogue,signature_data_matrix,method = "KLD",verbose=FALSE)
          E_hat_perc <- E_hat/matrix(apply(E_hat,2,sum),nrow = nrow(E_hat),ncol = ncol(E_hat),byrow = TRUE)*100
          E_hat[E_hat_perc<alpha] <- 0
        }else if(fit_methods[im]=="NNLS"){
          E_hat <- SignatureFit(catalogue,signature_data_matrix,method = "NNLS",verbose=FALSE)
          E_hat_perc <- E_hat/matrix(apply(E_hat,2,sum),nrow = nrow(E_hat),ncol = ncol(E_hat),byrow = TRUE)*100
          E_hat[E_hat_perc<alpha] <- 0
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
    }
    pos <- pos + 1
  }
}




#save or load results

if(recompute_Efit_table | !file.exists(Efit_file)){
  write.table(Efit_RMSE_table,file = Efit_file,quote = FALSE,sep = "\t",col.names = TRUE,row.names = TRUE)
  for (m in listOfClassificationMetrics){
    write.table(listOfClassificationMetrics_tables[[m]],
                file = listOfClassificationMetrics_files[[m]],
                quote = FALSE,sep = "\t",col.names = TRUE,row.names = TRUE)
  }
  for(alpha in thresholds_perc){
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
  for(alpha in thresholds_perc){
    for (m in listOfClassificationMetrics){
      listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][[m]] <- read.table(listOfClassificationMetrics_IndSig_files[[as.character(alpha)]][[m]],
                                                                                          sep = "\t",stringsAsFactors = FALSE,
                                                                                          header = TRUE,as.is = TRUE,check.names = FALSE)
    }
  }
}



#stats mean and standard error of mean

#Efit mean
Efit_table_mean <- matrix(0,nrow = length(thresholds_perc), ncol = length(fit_methods))
colnames(Efit_table_mean) <- fit_methods
row.names(Efit_table_mean) <- thresholds_perc
for (i in 1:n){
  Efit_table_mean <- Efit_table_mean + Efit_RMSE_table[,seq(i,ncol(Efit_RMSE_table),n)]
}
colnames(Efit_table_mean) <- fit_methods
Efit_table_mean <- Efit_table_mean/n
#Efit standard error
Efit_table_se <- matrix(0,nrow = length(thresholds_perc), ncol = length(fit_methods))
colnames(Efit_table_se) <- fit_methods
row.names(Efit_table_se) <- thresholds_perc
for (i in 1:n){
  Efit_table_se <- Efit_table_se + (Efit_RMSE_table[,seq(i,ncol(Efit_RMSE_table),n)] - Efit_table_mean)^2
}
Efit_table_se <- sqrt(Efit_table_se/n)/sqrt(n)


#Classification
listOfClassificationMetrics_tables_mean <- list()
listOfClassificationMetrics_tables_se <- list()
for (m in listOfClassificationMetrics){
  #Classification mean
  listOfClassificationMetrics_tables_mean[[m]] <- matrix(0,nrow = length(thresholds_perc), ncol = length(fit_methods))
  colnames(listOfClassificationMetrics_tables_mean[[m]]) <- fit_methods
  row.names(listOfClassificationMetrics_tables_mean[[m]]) <- thresholds_perc
  for (i in 1:n){
    listOfClassificationMetrics_tables_mean[[m]] <- listOfClassificationMetrics_tables_mean[[m]] + listOfClassificationMetrics_tables[[m]][,seq(i,ncol(listOfClassificationMetrics_tables[[m]]),n)]
  }
  colnames(listOfClassificationMetrics_tables_mean[[m]]) <- fit_methods
  listOfClassificationMetrics_tables_mean[[m]] <- listOfClassificationMetrics_tables_mean[[m]]/n
  #classification standard error
  listOfClassificationMetrics_tables_se[[m]] <- matrix(0,nrow = length(thresholds_perc), ncol = length(fit_methods))
  colnames(listOfClassificationMetrics_tables_se[[m]]) <- fit_methods
  row.names(listOfClassificationMetrics_tables_se[[m]]) <- thresholds_perc
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
    listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]] <- matrix(0,nrow = length(1:nsignatures), ncol = length(fit_methods))
    colnames(listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]]) <- fit_methods
    row.names(listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]]) <- colnames(P)
    for (i in 1:n){
      listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]] <- listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]] + listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][[m]][,seq(i,ncol(listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][[m]]),n)]
    }
    colnames(listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]]) <- fit_methods
    row.names(listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]]) <- colnames(P)
    listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]] <- listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]]/n
    #classification standard error
    listOfClassificationMetrics_IndSig_tables_se[[as.character(alpha)]][[m]] <- matrix(0,nrow = length(1:nsignatures), ncol = length(fit_methods))
    colnames(listOfClassificationMetrics_IndSig_tables_se[[as.character(alpha)]][[m]]) <- fit_methods
    row.names(listOfClassificationMetrics_IndSig_tables_se[[as.character(alpha)]][[m]]) <- colnames(P)
    for (i in 1:n){
      listOfClassificationMetrics_IndSig_tables_se[[as.character(alpha)]][[m]] <- listOfClassificationMetrics_IndSig_tables_se[[as.character(alpha)]][[m]] + (listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][[m]][,seq(i,ncol(listOfClassificationMetrics_IndSig_tables[[as.character(alpha)]][[m]]),n)] - listOfClassificationMetrics_IndSig_tables_mean[[as.character(alpha)]][[m]])^2
    }
    listOfClassificationMetrics_IndSig_tables_se[[as.character(alpha)]][[m]] <- sqrt(listOfClassificationMetrics_IndSig_tables_se[[as.character(alpha)]][[m]]/n)/sqrt(n)
  }
}




#plots

alphas <- thresholds_perc
nmf_algorithms <- c("Signature Fit")
nmf_algorithms_labels <- c("Signature Fit")
nsettings <- length(fit_methods)
#plot Efit
# current_labels <- paste(paste(rep(filter_labels,length(clustering_labels)),
#                               rep(rep(clustering_labels,each = length(filter_labels)),1),sep = " "),
#                         sep=" ")
file_out <- paste0(outdir,"2018_04_30_Efit_performance.jpg")
xlabel <- "threshold (%)"
jpeg(filename = file_out,width = 800*length(nmf_algorithms),height = 600,res = 100)
par(mfrow=c(1,length(nmf_algorithms)),oma=c(2,1,1.5,1))
par(mar=c(5,4,4,3))
for (j in 1:length(nmf_algorithms)){
  plot(1:length(alphas), Efit_table_mean[,1+(j-1)*nsettings],
       ylim=c(min(Efit_table_mean - Efit_table_se), max(1,max(Efit_table_mean + Efit_table_se))),
       t="l", xlab=xlabel, ylab="RMSE",
       main=nmf_algorithms_labels[j],
       lwd=2,col=colours_list[1],panel.first=grid(),xaxt='n')
  axis(1, at=1:length(alphas), labels=alphas)
  # hack: we draw arrows but with very special "arrowheads"
  if (n>1) arrows(1:length(alphas), Efit_table_mean[,1+(j-1)*nsettings]-Efit_table_se[,1+(j-1)*nsettings], 1:length(alphas), Efit_table_mean[,1+(j-1)*nsettings]+Efit_table_se[,1+(j-1)*nsettings], length=0.05, angle=90, code=3,lwd=2,col=colours_list[1])
  
  for(i in 2:nsettings){
    lines(1:length(alphas), Efit_table_mean[,i+(j-1)*nsettings],lwd=2,col=colours_list[i])
    if (n>1) arrows(1:length(alphas), Efit_table_mean[,i+(j-1)*nsettings]-Efit_table_se[,i+(j-1)*nsettings], 1:length(alphas), Efit_table_mean[,i+(j-1)*nsettings]+Efit_table_se[,i+(j-1)*nsettings], length=0.05, angle=90, code=3,lwd=2,col=colours_list[i])
    
  }
  
}
title(paste0("RMSE of Mutation Assignment (Mean and SE, n=",n,")"),outer = TRUE)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom",legend = fit_methods,horiz = TRUE,
       lwd = 2,col = colours_list[1:nsettings],
       xpd = TRUE,bty = "n",inset = c(0,0),cex = 1)
dev.off()


#plot Classification Metrics
# current_labels <- paste(paste(rep(filter_labels,length(clustering_labels)),
#                               rep(rep(clustering_labels,each = length(filter_labels)),1),sep = " "),
#                         sep=" ")
ml <-1
for(m in listOfClassificationMetrics){

  file_out <- paste0(outdir,"2018_04_30_E_metrics_performance",m,".jpg")
  xlabel <- "threshold (%)"

  jpeg(filename = file_out,width = 800*length(nmf_algorithms),height = 600,res = 100)
  par(mfrow=c(1,length(nmf_algorithms)),oma=c(2,1,1.5,1))
  par(mar=c(5,4,4,3))
  for (j in 1:length(nmf_algorithms)){
    plot(1:length(alphas), listOfClassificationMetrics_tables_mean[[m]][,1+(j-1)*nsettings],
         ylim=c(min(listOfClassificationMetrics_tables_mean[[m]] - listOfClassificationMetrics_tables_se[[m]]), max(0.5,max(listOfClassificationMetrics_tables_mean[[m]] + listOfClassificationMetrics_tables_se[[m]]))),
         t="l", xlab=xlabel, ylab=paste0(listOfClassificationMetrics_labels[[ml]]),
         main=nmf_algorithms_labels[j],
         lwd=2,col=colours_list[1],panel.first=grid(),xaxt='n')
    axis(1, at=1:length(alphas), labels=alphas)
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
  legend("bottom",legend = fit_methods,horiz = TRUE,
         lwd = 2,col = colours_list[1:nsettings],
         xpd = TRUE,bty = "n",inset = c(0,0),cex = 1)
  dev.off()
  ml <- ml+1
}


#plot 1 - Classification Metrics (Sensitivity/Specificity)
# current_labels <- paste(paste(rep(filter_labels,length(clustering_labels)),
#                               rep(rep(clustering_labels,each = length(filter_labels)),1),sep = " "),
#                         sep=" ")
ml <-1
for(m in listOfClassificationMetrics){
  mss <- listOfClassificationMetricsSensSpec[ml]

  file_out <- paste0(outdir,"2018_04_30_E_metrics_performance",mss,".jpg")
  xlabel <- "threshold (%)"

  jpeg(filename = file_out,width = 800*length(nmf_algorithms),height = 600,res = 100)
  par(mfrow=c(1,length(nmf_algorithms)),oma=c(2,1,1.5,1))
  par(mar=c(5,4,4,3))
  for (j in 1:length(nmf_algorithms)){
    plot(1:length(alphas), 1 - listOfClassificationMetrics_tables_mean[[m]][,1+(j-1)*nsettings],
         ylim=c(min(1 - (listOfClassificationMetrics_tables_mean[[m]] - listOfClassificationMetrics_tables_se[[m]])), max(0.01,max(1- (listOfClassificationMetrics_tables_mean[[m]] + listOfClassificationMetrics_tables_se[[m]])))),
         t="l", xlab=xlabel, ylab=paste0(listOfClassificationMetricsSensSpec_labels[ml]),
         main=nmf_algorithms_labels[j],
         lwd=2,col=colours_list[1],panel.first=grid(),xaxt='n')
    axis(1, at=1:length(alphas), labels=alphas)
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
  legend("bottom",legend = fit_methods,horiz = TRUE,
         lwd = 2,col = colours_list[1:nsettings],
         xpd = TRUE,bty = "n",inset = c(0,0),cex = 1)
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
# current_labels <- paste(paste(rep(filter_labels,length(clustering_labels)),
#                               rep(rep(clustering_labels,each = length(filter_labels)),1),sep = " "),
#                         sep=" ")
sig_lab <- paste0("Sig.",substr(colnames(P),11,12))
ml <-1
for(m in listOfClassificationMetrics){
  
  file_out <- paste0(outdir,"2018_04_30_E_metrics_performance_Signatures_",m,".jpg")
  xlabel <- "threshold (%)"
  
  jpeg(filename = file_out,width = 1000*length(nmf_algorithms),height = 600*length(alphas),res = 200)
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
  legend("bottom",legend = fit_methods,horiz = TRUE,
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
# current_labels <- paste(paste(rep(filter_labels,length(clustering_labels)),
#                               rep(rep(clustering_labels,each = length(filter_labels)),1),sep = " "),
#                         sep=" ")
sig_lab <- paste0("Sig.",substr(colnames(P),11,12))
ml <-1
for(m in listOfClassificationMetrics){
  mss <- listOfClassificationMetricsSensSpec[ml]

  file_out <- paste0(outdir,"2018_04_30_E_metrics_performance_Signatures_",mss,".jpg")
  xlabel <- "threshold (%)"

  jpeg(filename = file_out,width = 1000*length(nmf_algorithms),height = 600*length(alphas),res = 200)
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
  legend("bottom",legend = fit_methods,horiz = TRUE,
         fill = colours_list[1:nsettings],
         xpd = TRUE,bty = "n",inset = c(0,0),cex = 1)
  dev.off()
  ml <- ml+1
}




