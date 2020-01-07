library(signature.tools.lib)
source("../lib/logisticRegression.R")
source("../lib/logisticRegression.glm.R")
source("../lib/plotHRDetect_withCI.R")

#load conversion matrices
Q_toRef_subs <- read.table(file = "../../signature-tools/results/FitClusterAveragesToPanCanSigs/subs/Conversion_table_subs.tsv",
                           sep = "\t",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
Q_toRef_rearr <- read.table(file = "../../signature-tools/results/FitClusterAveragesToPanCanSigs/rearr/Conversion_table_rearr.tsv",
                            sep = "\t",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
colnames(Q_toRef_subs) <- paste0("Ref_subs_",colnames(Q_toRef_subs))
colnames(Q_toRef_rearr) <- paste0("Ref_rearr_",colnames(Q_toRef_rearr))


#---------functions
## Generate a random replicate of the cataloge 
# This method guarantees the total number of signatures is unchanged
generateRandMuts <- function(x){
  #consider the following method as a replacement
  full_r <- matrix(nrow = dim(x)[1],ncol = dim(x)[2])
  colnames(full_r) <- colnames(x)
  row.names(full_r) <- row.names(x)
  for (i in 1:ncol(x)){
    if(sum(x[,i]>0)){
      samples <- sample(1:nrow(x),size = sum(x[,i]),prob = x[,i]/sum(x[,i]),replace = TRUE)
      r <- unlist(lapply(1:nrow(x),function(p) sum(samples==p)))
    }else{ #no rearrangments found
      r <- x[,i]
    }
    names(r) <- rownames(x)
    full_r[,i] <- r
  }
  return(full_r)
}
#---------

buildBootstrapTables <- TRUE
computeBootstrapScores <- TRUE
plotHRDetectCOSMIC3_all <- TRUE

nboots_hr <- 1000
outdir <- "../results/HRDetectTraining_741_2018_10_16/HRDetect_with_Confidence/"
dir.create(outdir,showWarnings = FALSE,recursive = TRUE)

#Load Breast741 summary DF
load("../results/HRDetectTraining_741_2018_10_16/2018_10_16_trainingData741.rData")

#load models
load("../results/HRDetectTraining_741_2018_10_16/2018_07_26_Breast741_perf_and_finalModel_allfeatures.rData")
load("../results/HRDetectTraining_741_2018_10_16/COSMIC3_model.rData")
load("../results/HRDetectTraining_741_2018_10_16/2018_10_16_Breast741_perf_and_finalModel_individualfeatures.rData")

#load BRCA deficient samples
brca_status_table_v3 <- read.table("../data/BRCA_status/pcawg_samples_BRCA_status_v3_20171115.csv",header = TRUE,sep = "\t",as.is = TRUE)
biallelic_samples <- unique(brca_status_table_v3$sample_id[brca_status_table_v3$biallelic=="yes"])
biallelic_samples <- union(biallelic_samples,rownames(summaryDF)[summaryDF$BRCA_deficient==TRUE])


#get the organs to test
HRDetect_files_dir <- "../data/HRDetect_data_files/"
samples_files <- list.files(path = HRDetect_files_dir)
organs <- unique(sapply(samples_files,function(x) paste0(strsplit(x,"[.]")[[1]][2:2],collapse = "")))
organs <- setdiff(organs,"Breast")

#bootstrap_tables_file
bootstrap_tables_file <- paste0(outdir,"bootstrap_tables_nboots",nboots_hr,".rData")

#prepare the bootstraps in the correct format

if(buildBootstrapTables){
  #build the tables, for each organ/tissue
  bootstrap_tables <- list()
  #load signatures info
  choice_table <- list()
  choice_table[["subs"]] <- read.table("../../signature-tools/data/2018_02_18_fullExtraction_pancan_signatures_subs.txt",
                                       header = TRUE,as.is = TRUE,check.names = FALSE,sep = "\t")
  choice_table[["rearr"]] <- read.table("../../signature-tools/data/2018_02_18_fullExtraction_pancan_signatures_rearr.txt",
                                        header = TRUE,as.is = TRUE,check.names = FALSE,sep = "\t")
  
  columnsCOSMIC <- paste0("Signature.",1:30)
  columnsBreast560rearr <- paste0("RS",1:6)
 
  groups <- organs
  for (group in groups){ 
    message("Organ: ",group)
    bootstrap_tables[[group]] <- list()
    
    hrdet_table <- read.table(paste0("../data/HRDetect_data_files/HRDetect_data_table.",group,".vRefSig.tsv"),
                            sep = "\t",header = TRUE,check.names = FALSE,as.is = TRUE)
    hrdet_table <- hrdet_table[,dataColumnsOther]

    
    for (i in 1:nboots_hr){
      cat(".")
      
      #perturb HRD index and deletions classification
      #use rpois(n,lambda)
      tmp_table <- hrdet_table
      deletions_muts <- t(hrdet_table[,c("del.mh.count","del.rep.count","del.none.count")])
      ndeletions <- apply(deletions_muts,2,sum)
      samples_deletions_muts <- generateRandMuts(deletions_muts)
      samples_ndeletions <- apply(samples_deletions_muts,2,sum)
      samples_deletions_prop <- samples_deletions_muts/matrix(rep(samples_ndeletions,3),nrow = nrow(samples_deletions_muts),ncol = ncol(samples_deletions_muts),byrow = TRUE)
      tmp_table[,c("del.mh.count","del.rep.count","del.none.count")] <- t(samples_deletions_muts)
      tmp_table[,c("del.mh.prop","del.rep.prop","del.none.prop")] <- t(samples_deletions_prop)
      tmp_table[,"hrd"] <- rpois(nrow(hrdet_table),hrdet_table[,"hrd"])
      
      #prepare empty matrices
      tmp_subs <- as.data.frame(matrix(0,nrow = nrow(hrdet_table),ncol = length(dataColumnsSubs)+length(columnsCOSMIC),
                                       dimnames = list(rownames(hrdet_table),c(dataColumnsSubs,columnsCOSMIC))))
      tmp_rearr <- as.data.frame(matrix(0,nrow = nrow(hrdet_table),ncol = length(dataColumnsRearr)+length(columnsBreast560rearr),
                                        dimnames = list(rownames(hrdet_table),c(dataColumnsRearr,columnsBreast560rearr))))
      
      #subs
      subgroups <- c("")
      if(group=="Colorectal" | group=="Skin") subgroups <- c("A","B")
      for (subgroup in subgroups){
        project <- paste0(group,subgroup)
        p <- which(choice_table[["subs"]]$Organ==project)
        scaleYN <- choice_table[["subs"]]$scale[p]
        ns <- choice_table[["subs"]]$nsig[p]
        type_of_extraction <- "subs"
        #load exposures subs bootstraps
        load(paste0("../../signature-tools/results/subs/extraction/brunet/extraction_",project,"_subs_brunet/round_1/sig_",ns,"/SigFit_withBootstrap_Summary_mKLD_bfmCosSim_alpha-1_tr5_p0.05.rData"))
        subs_boots <- res$boot_list
        s1 <- sample(length(subs_boots),size = 1)
        current_subs <- t(subs_boots[[s1]])
        current_subs[is.na(current_subs)] <- 0
        #sparsity correction 5%
        sel <- t(apply(current_subs,1,function(x) x<sum(x)*0.05))
        current_subs[sel] <- 0
        
        #load conversion Q tables and convert subs
        signames <- paste0(project,"_",1:ncol(current_subs))
        colnames(current_subs) <- signames
        converted_subs <- as.matrix(current_subs) %*% as.matrix(Q_toRef_subs[signames,])
        samples_subs <- rownames(converted_subs)[rownames(converted_subs) %in% rownames(tmp_subs)]
        tmp_subs[samples_subs,dataColumnsSubs] <- converted_subs[samples_subs,]

        #add COSMIC columns
        if(group=="Breast741"){
          load(paste0("../../signature-tools/results/subs/extraction/brunet/extraction_Breast741_subs_brunet/cosmic12_fit/SigFit_withBootstrap_Summary_mKLD_bfmCosSim_alpha-1_tr5_p0.05.rData"))
        }else{
          load(paste0("../../signature-tools/results/subs/extraction/brunet/extraction_",project,"_subs_brunet/cosmic_fit/SigFit_withBootstrap_Summary_mKLD_bfmCosSim_alpha-1_tr5_p0.05.rData"))
        }
        subs_boots <- res$boot_list
        s1 <- sample(length(subs_boots),size = 1)
        current_subs <- t(subs_boots[[s1]])
        current_subs[is.na(current_subs)] <- 0
        if(group=="Breast741") colnames(current_subs) <- paste0("Signature.",c(1,2,3,5,6,8,13,17,18,20,26,30))
        #sparsity correction 5%
        sel <- t(apply(current_subs,1,function(x) x<sum(x)*0.05))
        current_subs[sel] <- 0
        samples_subs <- rownames(current_subs)[rownames(current_subs) %in% rownames(tmp_subs)]
        tmp_subs[samples_subs,colnames(current_subs)] <- current_subs[samples_subs,]
        
      }
      #load exposures rearr
      #find whether rearrs are available for this organ/tissue
      if (group %in% choice_table[["rearr"]]$Organ){
        p <- which(choice_table[["rearr"]]$Organ==group)
        scaleYN <- choice_table[["rearr"]]$scale[p]
        ns <- choice_table[["rearr"]]$nsig[p]
        type_of_extraction <- "rearr"
        #load exposures subs bootstraps
        load(paste0("../../signature-tools/results/rearr/extraction/brunet/extraction_",group,"_rearr_brunet/round_1/sig_",ns,"/SigFit_withBootstrap_Summary_mKLD_bfmCosSim_alpha-1_tr5_p0.05.rData"))
        rearr_boots <- res$boot_list
        s1 <- sample(length(rearr_boots),size = 1)
        current_rearr <- t(rearr_boots[[s1]])
        current_rearr[is.na(current_rearr)] <- 0
        #sparsity correction 5%
        sel <- t(apply(current_rearr,1,function(x) x<sum(x)*0.05))
        current_rearr[sel] <- 0
        
        #load conversion Q tables and convert subs
        signames <- paste0(project,"_",1:ncol(current_rearr))
        colnames(current_rearr) <- signames
        converted_rearr <- as.matrix(current_rearr) %*% as.matrix(Q_toRef_rearr[signames,])
        samples_rearr <- rownames(converted_rearr)[rownames(converted_rearr) %in% rownames(tmp_rearr)]
        tmp_rearr[samples_rearr,dataColumnsRearr] <- converted_rearr[samples_rearr,]
          
        #add COSMIC columns
        load(paste0("../../signature-tools/results/rearr/extraction/brunet/extraction_",group,"_rearr_brunet/breast560rearr_fit/SigFit_withBootstrap_Summary_mKLD_bfmCosSim_alpha-1_tr5_p0.05.rData"))
        rearr_boots <- res$boot_list
        s1 <- sample(length(rearr_boots),size = 1)
        current_rearr <- t(rearr_boots[[s1]])
        current_rearr[is.na(current_rearr)] <- 0
        #sparsity correction 5%
        sel <- t(apply(current_rearr,1,function(x) x<sum(x)*0.05))
        current_rearr[sel] <- 0
        samples_rearr <- rownames(current_rearr)[rownames(current_rearr) %in% rownames(tmp_rearr)]
        tmp_rearr[samples_rearr,colnames(current_rearr)] <- current_rearr[samples_rearr,]
        
      }#if no rearrangements were found then tmp_rearr is already set with rearrangements to zero
      bootstrap_tables[[group]][[i]] <- cbind(tmp_table,tmp_subs,tmp_rearr)
      # #write to file
      # write.table(tmp_table,file = paste0("../data/HRDetect_data_files/HRDetect_data_table.",group,".vBreast741.tsv"),
      #             sep = "\t",col.names = TRUE,row.names = TRUE,quote = FALSE)
    }
    cat("\n")
  }
  #save the bootstrap tables
  save(file = bootstrap_tables_file,bootstrap_tables,organs,nboots_hr)
}

#load the bootstrap tables
load(bootstrap_tables_file)

#bootstrap_HRDetect results file
bootstrap_scores_file <- paste0(outdir,"bootstrap_scores_nboots",nboots_hr,".rData")

#features names for HRDetect Davies et al. 2017
features_names <- c("del.mh.prop","Signature.3","RS3","RS5","hrd","Signature.8")


if(computeBootstrapScores){
  current_scores_table <- list()
  for (organ in organs){
    message("Organ: ",organ)
    #for all bootstrapped tables
    current_scores_table[[organ]] <- list()
    current_scores_table[[organ]][["Davies2017"]] <- NULL
    current_scores_table[[organ]][["Retrained"]] <- NULL
    current_scores_table[[organ]][["COSMIC3"]] <- NULL
    current_scores_table[[organ]][["Ref_subs_10"]] <- NULL
    for (i in 1:nboots_hr){
      cat(".")
      #run the models
      current_table <- bootstrap_tables[[organ]][[i]]
      current_table[is.na(current_table)] <- 0
      #HRDetect Davies et al. 2017
      current_scores <- t(applyHRDetectDavies2017(current_table,features_names,attachContributions = FALSE))
      current_scores_table[[organ]][["Davies2017"]] <- rbind(current_scores_table[[organ]][["Davies2017"]],current_scores)
      #Retrained
      current_scores <- t(logisticRegressionPredict(fittedModel = res.model.allfeatures,
                                xMatrix = current_table,
                                attachContributions = FALSE))
      current_scores_table[[organ]][["Retrained"]] <- rbind(current_scores_table[[organ]][["Retrained"]],current_scores)
      #COSMIC 3
      cosmic3table <- current_table[,"Signature.3",drop=FALSE]
      colnames(cosmic3table) <- "COSMIC.3"
      current_scores <- t(logisticRegressionPredict.glm(fittedModel = res.model.COSMIC.3,
                                                    xMatrix = cosmic3table,
                                                    attachContributions = FALSE))
      current_scores_table[[organ]][["COSMIC3"]] <- rbind(current_scores_table[[organ]][["COSMIC3"]],current_scores)
      #Ref_subs_10
      current_scores <- t(logisticRegressionPredict.glm(fittedModel = res.model[["Ref_subs_10"]],
                                                    xMatrix = current_table,
                                                    attachContributions = FALSE))
      current_scores_table[[organ]][["Ref_subs_10"]] <- rbind(current_scores_table[[organ]][["Ref_subs_10"]],current_scores)
      
    }
    cat("\n")
  }
  #save the bootstrap scores
  save(file = bootstrap_scores_file,current_scores_table,organs,nboots_hr)
}

#load the bootstrap scores
load(bootstrap_scores_file)

if(plotHRDetectCOSMIC3_all){
  for (organ in organs){
    message("Organ: ",organ)
    #HRDetect Davies et al. 2017
    #boxplot(current_scores_table[,order(apply(current_scores_table,2,median))])
    q_5_50_95 <- t(apply(current_scores_table[[organ]][["Davies2017"]],2,function(x) quantile(x,c(0.05,0.5,0.95))))
    filename <- paste0(outdir,"HRDetect_Score_bootstrap_",organ,"_Davies2017.jpg")
    par_title <- paste0("HRDetect Score Davies et al. 2017, ",organ,", median and 5%-95% quantiles interval (n=",nboots_hr,")")
    plot_HRDetect_withCI(filename,par_title,q_5_50_95,samplesBRCAdef = biallelic_samples)
    #Retrained
    q_5_50_95 <- t(apply(current_scores_table[[organ]][["Retrained"]],2,function(x) quantile(x,c(0.05,0.5,0.95))))
    filename <- paste0(outdir,"HRDetect_Score_bootstrap_",organ,"_Retrained.jpg")
    par_title <- paste0("HRDetect Score (Retrained), ",organ,", median and 5%-95% quantiles interval (n=",nboots_hr,")")
    plot_HRDetect_withCI(filename,par_title,q_5_50_95,samplesBRCAdef = biallelic_samples)
    #COSMIC3
    q_5_50_95 <- t(apply(current_scores_table[[organ]][["COSMIC3"]],2,function(x) quantile(x,c(0.05,0.5,0.95))))
    filename <- paste0(outdir,"HRDetect_Score_bootstrap_",organ,"_COSMIC3.jpg")
    par_title <- paste0("COSMIC3 Score, ",organ,", median and 5%-95% quantiles interval (n=",nboots_hr,")")
    plot_HRDetect_withCI(filename,par_title,q_5_50_95,samplesBRCAdef = biallelic_samples)
    #Ref_subs_10
    q_5_50_95 <- t(apply(current_scores_table[[organ]][["Ref_subs_10"]],2,function(x) quantile(x,c(0.05,0.5,0.95))))
    filename <- paste0(outdir,"HRDetect_Score_bootstrap_",organ,"_Ref_subs_10.jpg")
    par_title <- paste0("Ref_subs_10 Score, ",organ,", median and 5%-95% quantiles interval (n=",nboots_hr,")")
    plot_HRDetect_withCI(filename,par_title,q_5_50_95,samplesBRCAdef = biallelic_samples)
  }
}

