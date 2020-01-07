#Fork of SigHunter from Sandro Morganella 2017
#Andrea Degasperi, andrea.degasperi@sanger.ac.uk

SignatureExtraction <- function(cat, #matrix with samples as columns and channels as rows.
                                outFilePath, #path were the extraction output files should go. Remember to add "/" at the end of the path
                                blacklist=c(), #list of samples (column names) to ignore
                                nrepeats=10, #how many runs for each bootstrap (if filterBestOfEachBootstrap=TRUE with default params, only at most 10 runs within 0.1% of best will be considered, so nrepeats should be at least 10)
                                nboots=20, #how many bootstrapped catalogues to use
                                clusteringMethod="PAM", #choose among {"HC","PAM","MC"}, hierarchical clustering (HC), partitioning around the medoids (PAM) and  matched clustering (MC)
                                completeLinkageFlag=FALSE, #if clusteringMethod="HC", use complete linkage instead of default average linkage
                                useMaxMatching=TRUE, #if clusteringMethod="MC", use the assignment problem algorithm (match with max similarity) instead of the stable matching algorithm (any stable match)
                                filterBestOfEachBootstrap=TRUE, #if true only at most filterBest_nmaxtokeep of the nrepeats runs that are within filterBest_RTOL*best from the best are kept
                                filterBest_RTOL=0.001, #realtive tolerace from best fit to consider a run as good as the best
                                filterBest_nmaxtokeep=10, #max number of runs that should be kept that are within the relative tolerance from the best
                                nparallel=1, # how many processing units to use
                                nsig=c(3:15), # range of number of signatures to try
                                mut_thr=0, # threshold of mutations to remove empty/almost empty rows and columns
                                type_of_extraction="subs", #choose among {"subs","rearr","generic"}
                                project="extraction", # give a name to your project
                                parallel=FALSE, # set to TRUE to use parallel computation (Recommended)
                                nmfmethod="brunet", #choose among {"brunet","lee","nsNMF"}
                                removeDuplicatesInCatalogue = FALSE, #remove 0.99 cos sim similar samples
                                normaliseCatalogue = FALSE, # scale samples to sum to 1
                                plotCatalogue = FALSE, #also plot the catalogue, this may crash the library if the catalogue is too big, should work up to ~300 samples
                                plotResultsFromAllClusteringMethods=TRUE){ #if TRUE, all clustering methods are used and results are reported and plotted for all of them. If FALSE, only the requested clustering is reported
  
  library(NMF)
  library(foreach)
  library(doParallel)
  library(doMC)
  library(methods)
  source("../lib/matchedClustering.R")
  
  #------ code for debugging, normally should be commented out
  # setwd("~/sandbox/SignatureTools/lib")
  # catFile="../data/pancan_catalogues/project_BRCA.subs.catalogue.txt"
  # cat <- read.table(catFile, sep="\t", header=T, as.is=T, check.names = FALSE)
  # outFilePath="../results/pancan_only_tests/"
  # project <- "BRCA"
  # nboots=4
  # nrepeats=10
  # nparallel=3
  # #nsig=c(2:20)
  # nsig=c(4,5)
  # mut_thr=0
  # cos_sim_thr=0.90
  # merge_thr=0.8
  # type_of_extraction="rearr"
  # parallel=TRUE
  # ns=4
  # catFile=paste0("../data/pancan_catalogues/project_",project,".",type_of_extraction,".catalogue.txt")
  # outFilePath=paste0("../results/newExtr7_",project,"_",type_of_extraction,"/")
  #------ code for debugging, end
  
  #tmp for debug
  if (completeLinkageFlag) tmp_outFilePath <- outFilePath
  
  #remove blacklisted samples
  if(length(blacklist)>0){
    cat <- cat[,setdiff(colnames(cat),blacklist)]
  }
  
  group <- project
  
  registerDoMC(nparallel)
  dir.create(outFilePath,showWarnings = FALSE,recursive = TRUE)
  
  message("\n------------- START COMPUTATION ------------\n")
  
  ## Read the catalogue
  message("[STEP 1]: Reading and Preprocessing the Catalogue...", appendLF=F)
  if (type_of_extraction=="subs") cat <- sortCatalogue(cat)
  if(removeDuplicatesInCatalogue){
    cat <- removeSimilarCatalogueSamples(cat)
  }
  #Save catalogue used
  nsamples <- ncol(cat)
  if(plotCatalogue){
    cat_file <- paste0(outFilePath,"CatalogueUsedAfterPreprocessing_plot_",project,".jpg")
    if (type_of_extraction == "subs"){
      plotSubsSignatures(cat,cat_file,plot_sum = TRUE,overall_title = paste0("Catalogue - ",project," (",nsamples," samples)"))
    }else if (type_of_extraction == "rearr"){
      plotRearrSignatures(cat,cat_file,plot_sum = TRUE,overall_title = paste0("Catalogue - ",project," (",nsamples," samples)"))
    }else if (type_of_extraction == "generic"){
      plotGenericSignatures(cat,cat_file,plot_sum = TRUE,overall_title = paste0("Catalogue - ",project," (",nsamples," samples)"))
    }
  }
  cat_file <- paste0(outFilePath,"CatalogueUsedAfterPreprocessing_plot_",group,".txt")
  write.table(cat,file = cat_file,
              sep = "\t",quote = FALSE,row.names = TRUE,col.names = TRUE)
  
  nrow_cat <- nrow(cat)
  
  #remove channels if necessary
  all_rows_cat <- cat
  #cat <-  preprocessCatalgue(cat, mut_thr)
  channelsRemoved <- FALSE
  if(nrow(all_rows_cat)>nrow(cat)) channelsRemoved <- TRUE
  
  #normalised catalogue for computing replicates later
  ncat <- cat
  if(normaliseCatalogue) ncat <- normaliseSamples(cat)
  
  if(sum(nsig>ncol(cat))>0){
    nsig <- nsig[nsig[length(nsig)]<=ncol(cat)]
  }
  
  message("DONE")
  message("\t> ", ncol(cat), " Samples and ",nrow_cat, " Mutations Correclty Processed")
  

  
  if(mut_thr>0){
    message("\t> ", nrow_cat -nrow(cat), " Mutation(s) Channel(s) Removed (<", mut_thr, ")")
  }
  
  nmuts <- apply(cat, 2, sum)
  
  sample_names <- colnames(cat)
  
  ###### Find Signature ######
  message("\n[STEP 2]: Finding Signatures")
  
  ## Iterate Until Convergence
  complete <- FALSE
  round <- 1
  solutions <- matrix(0, 0, 4)
  colnames(solutions) <- c("Round", "Optim_Nsig", "Estimate", "NumUncorrelated") 
  
  prev_num_uncorr <- ncol(cat)
  
  while(complete==FALSE){
    message("\t> Running round ", round)
    
    ## Create Directory
    outDir <- paste(outFilePath, "round_", round, "/", sep="")
    cmd_res <- system(paste("mkdir -p", outDir))
    
    
    err <- c()
    cos_sim <- c()
    average_corr_smpls <- c()
    sil <- c()
    
    #additional metrics
    ave.RMSE <- c()
    sd.RMSE <- c()
    ave.KLD <- c()
    sd.KLD <- c()
    ave.RMSE.orig <- c()
    sd.RMSE.orig <- c()
    ave.KLD.orig <- c()
    sd.KLD.orig <- c()
    # ave.CosSim.hclust <- c()
    # ave.CosSim.PAM <- c()
    ave.SilWid.hclust <- c()
    ave.SilWid.PAM <- c()
    ave.SilWid.MC <- c()
    cophenetic.corr.hclust <- c()
    # proportion.tooSimilar.Signatures <- c()
    min.MinWCCS.hclust <- c()
    min.MinWCCS.PAM <- c()
    min.MinWCCS.MC <- c()
    max.MaxBCCS.hclust <- c()
    max.MaxBCCS.PAM <- c()
    max.MaxBCCS.MC <- c()
    mmcs_pam <- c()
    mmcs_hclust <- c()
    mmcs_MC <- c()
    
    ## Run Hunter on the specified range of signatures
    for(ns in nsig){
      
      #tmp for debug
      if (completeLinkageFlag) outFilePath <- tmp_outFilePath
      
      outNsDir <- paste(outDir, "sig_", ns, "/", sep="")
      cmd_res <- system(paste("mkdir -p ", outNsDir))
      
      message("\n> Running for ", ns,  " Signatures (", nboots, " bootstraps)")
      strt<-Sys.time()
      
      #---------Bootstraps start here
      
      #Define smoothing matrix S in case nsNMF is used
      if (nmfmethod=="nsNMF") S <- diag(ns)*0.5 + matrix(1,nrow=ns,ncol=ns)*0.5/ns
      
      #bootstraps file:
      bootstraps_file <- paste0(outFilePath,"bootstraps_",group,"_ns",ns,"_nboots",nboots,".Rdata")
      if (file.exists(bootstraps_file)){
        load(bootstraps_file)
        message("bootstraps file loaded")
      }else{
        ## Collect the solutions
        e_boot <- matrix(0, ncol(cat), 0)
        p_boot <- matrix(0, nrow(cat), 0)
        err_boot <- c()
        boot_tracker <- c()
        #boots_list <- list()
        boot_cat <- list()
        ## Run NMF on nboots bootsrapped catalogues
        if(parallel==TRUE){ ## Parallel
          #boots_list_tmp <- list()
          nseq <- ceiling(nboots/nparallel)
          #first get and store catalogues
          for (i in 1:(nseq*nparallel)){
            #boot_cat[[i]] <- apply(cat, 2, generateRandMuts)
            boot_cat[[i]] <- generateRandMuts(cat)
            if(normaliseCatalogue) boot_cat[[i]] <- normaliseSamples(boot_cat[[i]])
          }
          for(tt in 1:nseq){
            message(".",  appendLF=F)	
            boots_list <- foreach(i=1:nparallel) %dopar%{						
              rnd_cat <- boot_cat[[(tt-1)*nparallel + i]]			
              rnd_cat <- preprocessCatalgue(rnd_cat, mut_thr) #remove channels with 0 mutations
              nmf(rnd_cat, rank=ns,nrun=nrepeats,method = nmfmethod,.options="k-p")
            }
            #Extract data already to save memory
            for(i in 1:length(boots_list)){
              nmf_res <- boots_list[[i]]
              best_residual <- residuals(nmf_res)
              residuals_list <- c()
              if(length(nmf_res)==1){
                residuals_list <- residuals(nmf_res)
              }else{
                for (j in 1:length(nmf_res@.Data)){
                  #avoid numerical error (minimum cannot be less than 0)
                  residuals_list <- c(residuals_list,max(0,residuals(nmf_res@.Data[[j]])))
                }
              }
              if(filterBestOfEachBootstrap){
                #filter the best runs
                runsToChooseFrom <- which(best_residual*(1+filterBest_RTOL)>=residuals_list)
                #take at most filterBest_nmaxtokeep
                if (length(runsToChooseFrom)>filterBest_nmaxtokeep){
                  runsToChooseFrom <- sample(runsToChooseFrom,filterBest_nmaxtokeep)
                }
              }else{
                #just take all the runs
                runsToChooseFrom <- 1:length(nmf_res)
              }
              countRuns <- 1
              for (j in 1:length(nmf_res)){
                #keep at most 10 repeats that have residuals close to the best (within 0.1% more than residual)
                #if (best_residual*1.001>residuals_list[j] & countRuns <= 10) {
                if (j %in% runsToChooseFrom) {
                  if(length(nmf_res)==1){
                    #here is where I fix channels by adding back the channels I removed
                    bbb <- basis(nmf_res)
                    coln <- paste0("b",i,"r",j,"s",1:ns)
                    p_boot_tmp <- matrix(0,nrow = nrow(cat),ncol = ns,dimnames = list(rownames(cat),coln))
                    p_boot_tmp[rownames(bbb),] <- bbb
                    if (nmfmethod=="nsNMF"){
                      p_boot <- cbind(p_boot , p_boot_tmp %*% S)
                    }else{
                      p_boot <- cbind(p_boot , p_boot_tmp)
                    }
                    e_boot <- cbind(e_boot , t(coef(nmf_res)))
                    err_boot <- c(err_boot,residuals(nmf_res))
                  }else{
                    #here is where I fix channels by adding back the channels I removed
                    bbb <- basis(nmf_res@.Data[[j]])
                    coln <- paste0("b",i,"r",j,"s",1:ns)
                    p_boot_tmp <- matrix(0,nrow = nrow(cat),ncol = ns,dimnames = list(rownames(cat),coln))
                    p_boot_tmp[rownames(bbb),] <- bbb
                    if (nmfmethod=="nsNMF"){
                      p_boot <- cbind(p_boot , p_boot_tmp %*% S)
                    }else{
                      p_boot <- cbind(p_boot , p_boot_tmp)
                    }
                    e_boot <- cbind(e_boot , t(coef(nmf_res@.Data[[j]])))
                    err_boot <- c(err_boot,residuals(nmf_res@.Data[[j]]))
                  }

                  boot_tracker <- c(boot_tracker,rep((tt-1)*nparallel + i,ns))
                  countRuns <- countRuns + 1
                }
              }
              rm(nmf_res)
              gc()
            }
            rm(boots_list)
            gc()
          }
          #Now you need to cut so that you have only nboots and not more
          selection_of_boots <- boot_tracker <= nboots
          selection_of_boots_runs <- boot_tracker[seq(1,length(boot_tracker),ns)] <= nboots
          p_boot <- p_boot[,selection_of_boots]
          e_boot <- e_boot[,selection_of_boots]
          err_boot <- err_boot[selection_of_boots_runs]
          boot_tracker <- boot_tracker[selection_of_boots]
          # boot_pos <- 1
          # for (s in 1:nseq){
          #   for (p in 1:length(boots_list_tmp[[s]])){
          #     boots_list[[boot_pos]] <- boots_list_tmp[[s]][[p]]
          #     boot_pos <- boot_pos + 1
          #   }
          # }
          #boots_list <- unlist(boots_list_tmp)
        }else{ ## No Parallel	
          for(i in 1:nboots){
            message(".",  appendLF=F)			
            boot_cat[[i]] <- generateRandMuts(cat)
            if(normaliseCatalogue) boot_cat[[i]] <- normaliseSamples(boot_cat[[i]])
            rnd_cat <- preprocessCatalgue(boot_cat[[i]], mut_thr) #remove channels with 0 mutations
            nmf_res <- nmf(rnd_cat, rank=ns,nrun=nrepeats,method = nmfmethod,.options="k-p")
            #extract already to save memory
            best_residual <- residuals(nmf_res)
            residuals_list <- c()
            if(length(nmf_res)==1){
              residuals_list <- residuals(nmf_res)
            }else{
              for (j in 1:length(nmf_res@.Data)){
                #avoid numerical error (minimum cannot be less than 0)
                residuals_list <- c(residuals_list,max(0,residuals(nmf_res@.Data[[j]])))
              }
            }
            if(filterBestOfEachBootstrap){
              #filter the best runs
              runsToChooseFrom <- which(best_residual*(1+filterBest_RTOL)>=residuals_list)
              #take at most filterBest_nmaxtokeep
              if (length(runsToChooseFrom)>filterBest_nmaxtokeep){
                runsToChooseFrom <- sample(runsToChooseFrom,filterBest_nmaxtokeep)
              }
            }else{
              #just take all the runs
              runsToChooseFrom <- 1:length(nmf_res)
            }
            countRuns <- 1
            for (j in 1:length(nmf_res)){
              #keep at most 10 repeats that have residuals close to the best (within 0.1% more than residual)
              #if (best_residual*1.001>residuals_list[j] & countRuns <= 10) {
              if (j %in% runsToChooseFrom) {
                #here is where I fix channels by adding back the channels I removed
                bbb <- basis(nmf_res)
                coln <- paste0("b",i,"r",j,"s",1:ns)
                p_boot_tmp <- matrix(0,nrow = nrow(cat),ncol = ns,dimnames = list(rownames(cat),coln))
                p_boot_tmp[rownames(bbb),] <- bbb
                if(length(nmf_res)==1){
                  if (nmfmethod=="nsNMF"){
                    p_boot <- cbind(p_boot , p_boot_tmp %*% S)
                  }else{
                    p_boot <- cbind(p_boot , p_boot_tmp)
                  }
                  e_boot <- cbind(e_boot , t(coef(nmf_res)))
                  err_boot <- c(err_boot,residuals(nmf_res))
                }else{
                  #here is where I fix channels by adding back the channels I removed
                  bbb <- basis(nmf_res@.Data[[j]])
                  coln <- paste0("b",i,"r",j,"s",1:ns)
                  p_boot_tmp <- matrix(0,nrow = nrow(cat),ncol = ns,dimnames = list(rownames(cat),coln))
                  p_boot_tmp[rownames(bbb),] <- bbb
                  if (nmfmethod=="nsNMF"){
                    p_boot <- cbind(p_boot , p_boot_tmp %*% S)
                  }else{
                    p_boot <- cbind(p_boot , p_boot_tmp)
                  }
                  e_boot <- cbind(e_boot , t(coef(nmf_res@.Data[[j]])))
                  err_boot <- c(err_boot,residuals(nmf_res@.Data[[j]]))
                }
                boot_tracker <- c(boot_tracker,rep(i,ns))
                countRuns <- countRuns + 1
              }
            }
            rm(nmf_res)
            gc()
          }	
        }
        #how many solutions were saved?
        saved_nmf_runs <- ncol(p_boot)/ns
        save(file = bootstraps_file,ns,nboots,nrepeats,e_boot,p_boot,err_boot,cat,all_rows_cat,saved_nmf_runs,boot_tracker,boot_cat)
      }
      
      #now add back the missing channel and reset cat to all channels
      
      if(channelsRemoved){
        lmissing <- setdiff(rownames(all_rows_cat),rownames(cat))
        nmissing <- length(lmissing)
        newrows <- matrix(0,nrow = nmissing,ncol = ncol(p_boot))
        colnames(newrows) <- colnames(p_boot)
        rownames(newrows) <- lmissing
        p_boot <- rbind(p_boot,newrows)[rownames(all_rows_cat),]
        #cat <- all_rows_cat
      }
      
      #debug code below
      # ns <- 6
      # load(bootstraps_file)
      # nboots <- 1000
      # displace <- 0
      # e_boot <- e_boot[,1:(nboots*ns) + displace]
      # p_boot <- p_boot[,1:(nboots*ns) + displace]
      
      # ## Compute the average silhouette grouping all the computed solutions in ns clusters
      #sil <- c(sil, summary(silhouette(pam(p_boot, ns)))$avg.width)
      #above dissimilarity is not based on cosine similarity
      colnames_p <- c()
      sigs_p <- 1:ns
      for (i in 1:saved_nmf_runs){
        colnames_p <- c(colnames_p,paste(rep("sig",ns),sigs_p,rep("_run",ns),rep(i,ns),sep = ""))
      }
      colnames(p_boot) <- colnames_p
      names(boot_tracker) <- colnames_p
      #Load distance matrix if it already exists
      distMatrix_file <- paste0(outFilePath,"distMatrix_",group,"_ns",ns,"_nboots",nboots,".Rdata")
      if (file.exists(distMatrix_file)){
        load(distMatrix_file)
        message("distance matrix file loaded")
      }else{
        distMatrix <- 1 - computeCorrelation_parallel(p_boot,nparallel = nparallel,parallel = TRUE)
        save(file = distMatrix_file,distMatrix)
      }
      
      #change prefix (debug)
      if (completeLinkageFlag) outNsDir <- paste0(outNsDir,"completeLinkage")
      
      if(ns>1 & saved_nmf_runs>1){
        #Hierarchical clustering
        if (completeLinkageFlag) {
          fit_clust <- hclust(as.dist(distMatrix), method="complete") 
        }else{
          fit_clust <- hclust(as.dist(distMatrix), method="average") 
        }
        #Hierarchical clustering partitioning
        cut_res <- cutree(fit_clust,k = ns)
        #PAM partitioning
        clustering_result <- pam(as.dist(distMatrix), ns)
        cut_res.PAM <- clustering_result$clustering
        #Matched Clustering
        mc_file <- paste0(outNsDir,"matchedClustering_",group,"_ns",ns,"_nboots",nboots,".Rdata")
        if (file.exists(mc_file)){
          load(mc_file)
          message("matched clustering result loaded from file")
        }else{
          cut_res_MC <- matchedClustering(distMatrix,ns,maxMatch = useMaxMatching,parallel=TRUE,nparallel=nparallel)
          save(file = mc_file,cut_res_MC)
        }
        #get medoids
        medoids_hclust <- findMedoidsHclust(distMatrix,cut_res)
        medoids_MC <- findMedoidsHclust(distMatrix,cut_res_MC)
        medoids_PAM <- clustering_result$medoids
        
      }else if(ns==1 & saved_nmf_runs>1){
        cut_res <- rep(1,ncol(distMatrix))
        names(cut_res) <- colnames(distMatrix)
        cut_res.PAM <- cut_res
        cut_res_MC <- cut_res
        medoids_hclust <- findMedoidsHclust(distMatrix,cut_res)
        medoids_MC <- medoids_hclust
        medoids_PAM <- medoids_hclust
      }else if(ns>1 & saved_nmf_runs==1){
        cut_res <- seq(1,ncol(distMatrix))
        names(cut_res) <- colnames(distMatrix)
        cut_res.PAM <- cut_res
        cut_res_MC <- cut_res
        medoids_hclust <- findMedoidsHclust(distMatrix,cut_res)
        medoids_MC <- medoids_hclust
        medoids_PAM <- medoids_hclust
      }else if(ns==1 & saved_nmf_runs==1){
        cut_res <- 1
        names(cut_res) <- colnames(distMatrix)
        cut_res.PAM <- cut_res
        cut_res_MC <- cut_res
        medoids_hclust <- colnames(distMatrix)
        medoids_MC <- medoids_hclust
        medoids_PAM <- medoids_hclust
      }
      
      norm_p_boot <- p_boot/matrix(data=rep(apply(p_boot,2,sum),nrow(p_boot)),nrow = nrow(p_boot),byrow = TRUE)
      
      #mean and sd
      #set up mean and sd of signatures
      mean_signatures <- list()
      sd_signatures <- list()
      for (cm in c("HC","PAM","MC")){
        mean_signatures[[cm]] <- matrix(NA,nrow = nrow(p_boot),ncol = ns)
        sd_signatures[[cm]] <- matrix(NA,nrow = nrow(p_boot),ncol = ns)
        colnames(mean_signatures[[cm]]) <- paste0("S",1:ns)
        colnames(sd_signatures[[cm]]) <- paste0("S",1:ns)
        row.names(mean_signatures[[cm]]) <- row.names(p_boot)
        row.names(sd_signatures[[cm]]) <- row.names(p_boot)
      }
      
      for (cm in c("HC","PAM","MC")){
        if(cm=="HC"){
          partitions <- cut_res
        }else if(cm=="PAM"){
          partitions <- cut_res.PAM
        }else if(cm=="MC"){
          partitions <- cut_res_MC
        }
        for(nsi in 1:ns){
          pboot_dim <- dim(norm_p_boot[,partitions==nsi,drop=FALSE])
          if(pboot_dim[2]==1){
            mean_signatures[[cm]][,nsi] <- norm_p_boot[,partitions==nsi]
            sd_signatures[[cm]][,nsi] <- 0
          }else{
            mean_signatures[[cm]][,nsi] <- apply(norm_p_boot[,partitions==nsi],1,mean)
            sd_signatures[[cm]][,nsi] <- apply(norm_p_boot[,partitions==nsi],1,sd)
            sd_signatures[[cm]][is.na(sd_signatures[[cm]][,nsi]),nsi] <- 0
          }
        }
      }
      
      if(ns>1 & saved_nmf_runs>1){
        #Compute Silhouettes
        sil_hclust <- summary(silhouette(cut_res,as.dist(distMatrix)))
        sil_pam <- summary(silhouette(cut_res,as.dist(distMatrix)))
        sil_MC <- summary(silhouette(cut_res_MC,as.dist(distMatrix)))
        plotWithinClusterSilWidth(sil_hclust,sil_pam,sil_MC,outNsDir,group,ns,nboots)
        
        #compute within cluster cos similiraty distance
        cosSimPAM <- withinClusterCosSim(clustering_result$clustering,distMatrix,parallel = TRUE)
        cosSimHClust <- withinClusterCosSim(cut_res,distMatrix,parallel = TRUE)
        cosSimMC <- withinClusterCosSim(cut_res_MC,distMatrix,parallel = TRUE)
        plotWithinClusterCosSim(cosSimHClust,cosSimPAM,cosSimMC,outNsDir,group,ns,nboots)
        
        #compute cophenetic correlation
        coph <- cor(cophenetic(fit_clust),as.dist(distMatrix))
        
        #compute the proportion of the solutions that contain signatures that are too similar
        #propTooSimilar <- computePropTooSimilar(distMatrix,saved_nmf_runs,ns)
        
        #Save metrics
        #ave.CosSim.hclust <- c(ave.CosSim.hclust,mean(cosSimHClust))
        #ave.CosSim.PAM <- c(ave.CosSim.PAM,mean(cosSimPAM))
        ave.SilWid.hclust <- c(ave.SilWid.hclust,sil_hclust$avg.width)
        ave.SilWid.PAM <- c(ave.SilWid.PAM,sil_pam$avg.width)
        ave.SilWid.MC <- c(ave.SilWid.MC,sil_MC$avg.width)
        cophenetic.corr.hclust <- c(cophenetic.corr.hclust,coph)
        #proportion.tooSimilar.Signatures <- c(proportion.tooSimilar.Signatures,propTooSimilar)
        
        #save metrics
        additional_perf_file <- paste0(outNsDir,"Sigs_WithinClusterPerf_",group,"_ns",ns,"_nboots",nboots,".rData")
        save(file = additional_perf_file,sil_pam,sil_hclust,sil_MC,cosSimPAM,cosSimHClust,cosSimMC)
        #plotHierarchicalCluster(fit_clust,outNsDir,group,ns,nboots)
        
        #Max Medoids Cosine Similarity
        mmcs_pam <- c(mmcs_pam,max(medoids_cosSimMatrix(p_boot,medoids_PAM) - diag(ns)))
        mmcs_hclust <- c(mmcs_hclust,max(medoids_cosSimMatrix(p_boot,medoids_hclust) - diag(ns)))
        mmcs_MC <- c(mmcs_MC,max(medoids_cosSimMatrix(p_boot,medoids_MC) - diag(ns)))
        
      }else{
        ave.SilWid.hclust <- c(ave.SilWid.hclust,NA)
        ave.SilWid.PAM <- c(ave.SilWid.PAM,NA)
        ave.SilWid.MC <- c(ave.SilWid.MC,NA)
        cophenetic.corr.hclust <- c(cophenetic.corr.hclust,NA)
        mmcs_pam <- c(mmcs_pam,NA)
        mmcs_hclust <- c(mmcs_hclust,NA)
        mmcs_MC <- c(mmcs_MC,NA)
      }
      
      # rmse_list <- c()
      # kld_list <- c()
      # for (i in 1:saved_nmf_runs){
      #   selection <- ((i-1)*ns+1):(i*ns)
      #   current_p <- p_boot[,selection]
      #   current_e <- e_boot[,selection]
      #   rmse_list <- c(rmse_list,sqrt(sum((cat - current_p %*% t(current_e))^2)/(dim(cat)[1]*dim(cat)[2])))
      #   kld_list <- c(kld_list,KLD(cat,current_p %*% t(current_e)))
      # }
      # # #Save errors
      # ave.RMSE <- c(ave.RMSE,mean(rmse_list))
      # sd.RMSE <- c(sd.RMSE,sd(rmse_list))
      # ave.KLD <- c(ave.KLD,mean(kld_list))
      # sd.KLD <- c(sd.KLD,sd(kld_list))
      # 
      
      #error
      rmse_list <- c()
      kld_list <- c()
      rmse_orig_list <- c()
      kld_orig_list <- c()
      for (i in 1:saved_nmf_runs){
        selection <- ((i-1)*ns+1):(i*ns)
        current_cat <- boot_cat[[boot_tracker[selection][1]]]
        if(channelsRemoved){
          current_p <- p_boot[!row.names(p_boot) %in% lmissing,selection]
        }else{
          current_p <- p_boot[,selection]
        }
        current_e <- e_boot[,selection]
        reconstructed_cat <- current_p %*% t(current_e)
        rmse_list <- c(rmse_list,sqrt(sum((current_cat - reconstructed_cat)^2)/(dim(current_cat)[1]*dim(current_cat)[2])))
        kld_list <- c(kld_list,KLD(current_cat,reconstructed_cat))
        
        rmse_orig_list <- c(rmse_orig_list,sqrt(sum((ncat - reconstructed_cat)^2)/(dim(ncat)[1]*dim(ncat)[2])))
        kld_orig_list <- c(kld_orig_list,KLD(ncat,reconstructed_cat))
      }
      # #Save errors
      ave.RMSE <- c(ave.RMSE,mean(rmse_list))
      sd.RMSE <- c(sd.RMSE,sd(rmse_list))
      ave.KLD <- c(ave.KLD,mean(kld_list))
      sd.KLD <- c(sd.KLD,sd(kld_list))
      # 
      ave.RMSE.orig <- c(ave.RMSE.orig,mean(rmse_orig_list))
      sd.RMSE.orig <- c(sd.RMSE.orig,sd(rmse_orig_list))
      ave.KLD.orig <- c(ave.KLD.orig,mean(kld_orig_list))
      sd.KLD.orig <- c(sd.KLD.orig,sd(kld_orig_list))
      
      #use requested medoids
      medoids_final <- medoids_hclust
      if(clusteringMethod=="PAM") medoids_final <- medoids_PAM
      if(clusteringMethod=="MC") medoids_final <- medoids_MC
      
      #plot signatures and compute similarity to known signatures
      clustering_list <- ""
      clustering_list_tag <- clusteringMethod
      if (plotResultsFromAllClusteringMethods) clustering_list <- c("","_HC","_PAM","_MC")
      if (plotResultsFromAllClusteringMethods) clustering_list_tag <- c(clusteringMethod,"HC","PAM","MC")
      
      for (cli in 1:length(clustering_list)){
        cl <- clustering_list[cli]
        cl_tag <- clustering_list_tag[cli]
        signature_names <- paste0("S",1:length(medoids_final))
        
        if (cl==""){
          signature_data_matrix <- p_boot[,medoids_final,drop=FALSE]/matrix(data=rep(apply(p_boot[,medoids_final,drop=FALSE],2,sum),nrow(p_boot)),nrow = nrow(p_boot),byrow = TRUE)
        }else if(cl=="_HC"){
          signature_data_matrix <- p_boot[,medoids_hclust,drop=FALSE]/matrix(data=rep(apply(p_boot[,medoids_hclust,drop=FALSE],2,sum),nrow(p_boot)),nrow = nrow(p_boot),byrow = TRUE)
        }else if(cl=="_PAM"){
          signature_data_matrix <- p_boot[,medoids_PAM,drop=FALSE]/matrix(data=rep(apply(p_boot[,medoids_PAM,drop=FALSE],2,sum),nrow(p_boot)),nrow = nrow(p_boot),byrow = TRUE)
        }else if(cl=="_MC"){
          signature_data_matrix <- p_boot[,medoids_MC,drop=FALSE]/matrix(data=rep(apply(p_boot[,medoids_MC,drop=FALSE],2,sum),nrow(p_boot)),nrow = nrow(p_boot),byrow = TRUE)
        }
        
        colnames(signature_data_matrix) <- signature_names
        row.names(signature_data_matrix) <- row.names(p_boot)
        subs_file <- paste0(outNsDir,"Sigs_plot_",group,"_ns",ns,"_nboots",nboots,cl,".jpg")
        if (type_of_extraction == "subs"){
          #plotSubsSignatures(signature_data_matrix,subs_file,plot_sum = FALSE,overall_title = paste0("Medoids Signatures when extracting ",ns))
          plotSubsSignatures_withMeanSd(signature_data_matrix,mean_signatures[[cl_tag]],sd_signatures[[cl_tag]],subs_file,plot_sum = FALSE,overall_title = paste0("Medoids Signatures when extracting ",ns))
        }else if (type_of_extraction == "rearr"){
          #plotRearrSignatures(signature_data_matrix,subs_file,plot_sum = FALSE,overall_title = paste0("Medoids Signatures when extracting ",ns))
          plotRearrSignatures_withMeanSd(signature_data_matrix,mean_signatures[[cl_tag]],sd_signatures[[cl_tag]],subs_file,plot_sum = FALSE,overall_title = paste0("Medoids Signatures when extracting ",ns))
        }else if (type_of_extraction == "generic"){
          #plotGenericSignatures(signature_data_matrix,subs_file,plot_sum = FALSE,overall_title = paste0("Medoids Signatures when extracting ",ns))
          plotGenericSignatures_withMeanSd(signature_data_matrix,mean_signatures[[cl_tag]],sd_signatures[[cl_tag]],subs_file,plot_sum = FALSE,overall_title = paste0("Medoids Signatures when extracting ",ns))
        }
        subs_file <- paste0(outNsDir,"Sigs_plot_",group,"_ns",ns,"_nboots",nboots,cl,".tsv")
        write.table(signature_data_matrix,file = subs_file,
                    sep = "\t",quote = FALSE,row.names = TRUE,col.names = TRUE)
        
        #similarity to known signatures
        if(type_of_extraction=="subs"){
          #find the most similar cosmic signatures
          res_cosmic <- findClosestCOSMIC30_withSimilarity(signature_data_matrix) 
          res_cosmicComb <- findClosestCOSMIC30andCombinations_withSimilarity(signature_data_matrix) 
          res_cosmic_table <- data.frame(res_cosmic,res_cosmicComb)
          cosmic_file <- paste0(outNsDir,"Sigs_cosmicSimilar_",group,"_ns",ns,"_nboots",nboots,cl,".tsv")
          write.table(res_cosmic_table,file = cosmic_file,
                      sep = "\t",quote = FALSE,row.names = TRUE,col.names = TRUE)
        }else if (type_of_extraction=="rearr"){
          res_rearrBreast560 <- findClosestRearrSigsBreast560_withSimilarity(signature_data_matrix)
          res_rearrBreast560_table <- data.frame(res_rearrBreast560)
          res_rearrBreast560_file <- paste0(outNsDir,"Sigs_rearrBreast560Similar_",group,"_ns",ns,"_nboots",nboots,cl,".tsv")
          write.table(res_rearrBreast560_table,file = res_rearrBreast560_file,
                      sep = "\t",quote = FALSE,row.names = TRUE,col.names = TRUE)
        }
      }
      
      if (ns>1){
        #More metrics
        #spread of the clusters: minimum within cluster cosine similarity (MinWCCS)
        MinWCCS.hclust <- minWithinClusterCosSim(clustering = cut_res,distMatrix = distMatrix,parallel = TRUE)
        names(MinWCCS.hclust) <- signature_names
        MinWCCS.MC <- minWithinClusterCosSim(clustering = cut_res_MC,distMatrix = distMatrix,parallel = TRUE)
        names(MinWCCS.MC) <- signature_names
        MinWCCS.PAM <- minWithinClusterCosSim(clustering = clustering_result$clustering,distMatrix = distMatrix,parallel = TRUE)
        names(MinWCCS.PAM) <- signature_names
        
        MinWCCS_file <- paste0(outNsDir,"MinWithinClusterCosSim_",group,"_ns",ns,"_nboots",nboots,".tsv")
        write.table(data.frame(MinWCCS.hclust=MinWCCS.hclust,MinWCCS.PAM=MinWCCS.PAM,MinWCCS.MC=MinWCCS.MC),file = MinWCCS_file,
                    sep = "\t",quote = FALSE,row.names = TRUE,col.names = TRUE)
        min.MinWCCS.hclust <- c(min.MinWCCS.hclust,min(MinWCCS.hclust))
        min.MinWCCS.PAM <- c(min.MinWCCS.PAM,min(MinWCCS.PAM))
        min.MinWCCS.MC <- c(min.MinWCCS.MC,min(MinWCCS.MC))
        
        #cluster neighbours: maximum between cluster cosine similarity (MaxBCCS)
        MaxBCCS.hclust <- maxBetweenClustersCosSim(clustering = cut_res,distMatrix = distMatrix,parallel = TRUE)
        colnames(MaxBCCS.hclust) <- signature_names
        row.names(MaxBCCS.hclust) <- signature_names
        MaxBCCS.MC <- maxBetweenClustersCosSim(clustering = cut_res_MC,distMatrix = distMatrix,parallel = TRUE)
        colnames(MaxBCCS.MC) <- signature_names
        row.names(MaxBCCS.MC) <- signature_names
        MaxBCCS.PAM <- maxBetweenClustersCosSim(clustering = clustering_result$clustering,distMatrix = distMatrix,parallel = TRUE)
        colnames(MaxBCCS.PAM) <- signature_names
        row.names(MaxBCCS.PAM) <- signature_names
        
        MaxBCCS.hclust_file <- paste0(outNsDir,"MaxBetweenClusterCosSim.hclust_",group,"_ns",ns,"_nboots",nboots,".tsv")
        write.table(MaxBCCS.hclust,file = MaxBCCS.hclust_file,
                    sep = "\t",quote = FALSE,row.names = TRUE,col.names = TRUE)
        MaxBCCS.PAM_file <- paste0(outNsDir,"MaxBetweenClusterCosSim.PAM_",group,"_ns",ns,"_nboots",nboots,".tsv")
        write.table(MaxBCCS.PAM,file = MaxBCCS.PAM_file,
                    sep = "\t",quote = FALSE,row.names = TRUE,col.names = TRUE)
        MaxBCCS.MC_file <- paste0(outNsDir,"MaxBetweenClusterCosSim.MC_",group,"_ns",ns,"_nboots",nboots,".tsv")
        write.table(MaxBCCS.MC,file = MaxBCCS.MC_file,
                    sep = "\t",quote = FALSE,row.names = TRUE,col.names = TRUE)
        
        max.MaxBCCS.hclust <- c(max.MaxBCCS.hclust,max(MaxBCCS.hclust - diag(nrow(MaxBCCS.hclust)))) 
        max.MaxBCCS.PAM <- c(max.MaxBCCS.PAM,max(MaxBCCS.PAM - diag(nrow(MaxBCCS.PAM))))
        max.MaxBCCS.MC <- c(max.MaxBCCS.MC,max(MaxBCCS.MC - diag(nrow(MaxBCCS.MC)))) 
      }else if (ns==1){
        min.MinWCCS.hclust <- c(min.MinWCCS.hclust,NA)
        min.MinWCCS.PAM <- c(min.MinWCCS.PAM,NA)
        min.MinWCCS.MC <- c(min.MinWCCS.MC,NA)
        max.MaxBCCS.hclust <- c(max.MaxBCCS.hclust,NA) 
        max.MaxBCCS.PAM <- c(max.MaxBCCS.PAM,NA)
        max.MaxBCCS.MC <- c(max.MaxBCCS.MC,NA)
      }      
      #-------------------
      #----- clean this ns
      #-------------------
      rm(p_boot)
      rm(e_boot)
      rm(distMatrix)
      gc()
    }
    
    #Collect metrics
    # overall_metrics <- data.frame(nsig,
    #                               ave.RMSE,
    #                               sd.RMSE,
    #                               ave.SilWid.hclust,
    #                               ave.SilWid.PAM,
    #                               cophenetic.corr.hclust,
    #                               min.MinWCCS.hclust,
    #                               min.MinWCCS.PAM,
    #                               max.MaxBCCS.hclust,
    #                               max.MaxBCCS.PAM,
    #                               mmcs_hclust,
    #                               mmcs_pam,
    #                               ave.KLD,
    #                               sd.KLD)
    overall_metrics <- data.frame(nsig,       #1
                                  ave.RMSE,   #2
                                  sd.RMSE,    #3
                                  ave.RMSE.orig, #4
                                  sd.RMSE.orig,  #5
                                  ave.KLD,       #6
                                  sd.KLD,        #7
                                  ave.KLD.orig,  #8
                                  sd.KLD.orig,   #9
                                  ave.SilWid.hclust, #10
                                  ave.SilWid.PAM,    #11
                                  ave.SilWid.MC,     #12
                                  cophenetic.corr.hclust, #13
                                  min.MinWCCS.hclust, #14
                                  min.MinWCCS.PAM,    #15
                                  min.MinWCCS.MC,     #16
                                  max.MaxBCCS.hclust, #17
                                  max.MaxBCCS.PAM,    #18
                                  max.MaxBCCS.MC,     #19
                                  mmcs_hclust,        #20
                                  mmcs_pam,           #21
                                  mmcs_MC)            #22
    
    #write table with all the metrics
    overall_metrics_file <- paste0(outFilePath,"Sigs_OverallMetrics_",group,"_nboots",nboots,".tsv")
    write.table(overall_metrics,file = overall_metrics_file,
                sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE)
    
    whattoplot_hclust <- c(10,14,17,20,13)
    whattoplot_PAM <- c(11,15,18,21)    
    whattoplot_MC <- c(12,16,19,22)
    
    #plot requested medoids
    whattoplot_final <- whattoplot_hclust
    if(clusteringMethod=="PAM") whattoplot_final <- whattoplot_PAM
    if(clusteringMethod=="MC") whattoplot_final <- whattoplot_MC
    
    for (cl in clustering_list){

      overall_metrics_file <- paste0(outFilePath,"Sigs_OverallMetrics_",group,"_nboots",nboots,cl,".jpg")
      #read file for debugging
      # overall_metrics <- read.table(file = overall_metrics_file,sep = "\t",header = TRUE,as.is = TRUE)
      #specify which columns to plot. nsig and ave.RMSE will be used already, just specify others
      # whattoplot <- c(4,6,7,9,11)
      if(cl==""){
        whattoplot <- whattoplot_final
      }else if(cl=="_HC"){
        whattoplot <- whattoplot_hclust
      }else if(cl=="_PAM"){
        whattoplot <- whattoplot_PAM
      }else if(cl=="_MC"){
        whattoplot <- whattoplot_MC
      }
      plotOverallMetrics(overall_metrics,whattoplot,overall_metrics_file,group,nboots,nmfmethod)
    }
    complete <- TRUE
  }
  

  message("\n\nResults can be found in ", outFilePath)
  
  message("\n------------- COMPUTATION  SUCCESSFULLY COMPLETED ------------\n")
  
}

#------------------------------------------------
#------- Auxiliary functions ---------------------
#------------------------------------------------


## Generate a random replicate of the cataloge 
# This method guarantees the total number of signatures is unchanged
generateRandMuts <- function(x){
  #consider the following method as a replacement
  full_r <- matrix(nrow = dim(x)[1],ncol = dim(x)[2])
  colnames(full_r) <- colnames(x)
  row.names(full_r) <- row.names(x)
  for (i in 1:ncol(x)){
    samples <- sample(1:nrow(x),size = sum(x[,i]),prob = x[,i]/sum(x[,i]),replace = TRUE)
    r <- unlist(lapply(1:nrow(x),function(p) sum(samples==p)))
    names(r) <- rownames(x)
    full_r[,i] <- r
  }
  return(full_r)
}

## Remove unused Rows and Columns from the Catalogue
## Remove Mutations with small numbers
preprocessCatalgue <- function(d, mut_thr){
  
  ## Remove Mutations
  nmut <- apply(d, 1, sum)
  nmut <- nmut/sum(nmut)
  pos <- which(nmut<=mut_thr)
  if(length(pos)>0){
    d <- d[-pos,]
  }
  return(d)
}


sortCatalogue <- function(cat){
  all_bp <- c("A", "C", "G", "T")
  pyr_muts <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  nm <- c()
  for(m in pyr_muts){
    for(a in all_bp){
      for(b in all_bp){
        nm <- c(nm, paste(a, "[",m, "]", b, sep=""))
      }
    }
  }
  return(cat[nm,,drop=FALSE])
  
}



computeCorrelation <- function(x){
  if (ncol(x)==1){
    #if there is only one column, correlation matrix is 1
    return(matrix(1, 1, 1)) 
  }
  out <- matrix(NA, ncol(x), ncol(x))
  #diagonal is 1
  for(i in 1:ncol(x)){
    out[i,i] <- 1
  }
  colnames(out) <- colnames(x)
  rownames(out) <- colnames(x)
  for(i in 2:ncol(x)){
    for(j in 1:(i-1)){ #up to i-1, diag already set to 1
      #message(i, " ", j)
      out[i,j] <- cos.sim(as.numeric(x[,i]), as.numeric(x[,j]))
      out[j,i] <- out[i,j] #upper triangular is the same
    }
  }
  return(out)
}

computeCorrelation_parallel <- function(x,nparallel=1,parallel=FALSE){
  if (ncol(x)==1){
    #if there is only one column, correlation matrix is 1
    return(matrix(1, 1, 1)) 
  }
  out <- matrix(NA, ncol(x), ncol(x))
  #diagonal is 1
  for(i in 1:ncol(x)){
    out[i,i] <- 1
  }
  colnames(out) <- colnames(x)
  rownames(out) <- colnames(x)
  #correlation matrix is symmetric
  if (parallel){
    library(foreach)
    library(doParallel)
    library(doMC)
    registerDoMC(nparallel)
    par_res <- foreach(i=2:ncol(x)) %dopar% {
      current_res <- c()
      for(j in 1:(i-1)){ #up to i-1, diag already set to 1
        #message(i, " ", j)
        current_res <- c(current_res,cos.sim(as.numeric(x[,i]), as.numeric(x[,j])))
      }
      current_res
    }
    for(i in 2:ncol(x)){
      out[i,1:(i-1)] <- par_res[[i-1]]
      out[1:(i-1),i] <- par_res[[i-1]]
    }
  }else{
    for(i in 2:ncol(x)){
      for(j in 1:(i-1)){ #up to i-1, diag already set to 1
        #message(i, " ", j)
        out[i,j] <- cos.sim(as.numeric(x[,i]), as.numeric(x[,j]))
        out[j,i] <- out[i,j] #upper triangular is the same
      }
    }
  }
  return(out)
}

#compute the cosine similarity
cos.sim <- function(a, b){
  return( sum(a*b)/sqrt(sum(a^2)*sum(b^2)) )
} 

#function to compute the average within cluster cosine similarity,
#uses as input the distance matrix constructed as 1 - cosSimMatrix
withinClusterCosSim <- function(clustering,distMatrix,parallel){
  clusters <- unique(clustering)
  clusters <- clusters[order(clusters)]
  res <- c()
  if (parallel==TRUE){
    res <- foreach(i=clusters) %dopar%{
      samplesInCluster <- names(clustering[clustering==i])
      nSamples <- length(samplesInCluster)
      meanCosSim <- 0
      if(nSamples==1){
        meanCosSim <- 1
      }else{
        combinations <- nSamples*(nSamples - 1)/2
        #for efficiency, instead of summing lots of (1 - dist)
        #start from 1*combinations and subtract the dists
        sumCosSim <- combinations 
        for(j in 1:(nSamples-1)){
          #for(w in (j+1):(nSamples)){
            sumCosSim <- sumCosSim - sum(distMatrix[samplesInCluster[j],samplesInCluster[1:nSamples > j]])
          #}
        }
        meanCosSim <- sumCosSim/combinations
      }
      meanCosSim
    }
    res <- unlist(res)
  }else{
    for(i in clusters){
      samplesInCluster <- names(clustering[clustering==i])
      nSamples <- length(samplesInCluster)
      combinations <- nSamples*(nSamples - 1)/2
      #for efficiency, instead of summing lots of (1 - dist)
      #start from 1*combinations and subtract the dists
      sumCosSim <- combinations 
      for(j in 1:(nSamples-1)){
        #for(w in (j+1):(nSamples)){
          sumCosSim <- sumCosSim - sum(distMatrix[samplesInCluster[j],samplesInCluster[1:nSamples > j]])
        #}
      }
      res <- c(res,sumCosSim/combinations)
    }
  }
  return(res)
}

#function to compute the minimum within cluster cosine similarity,
#uses as input the distance matrix constructed as 1 - cosSimMatrix
minWithinClusterCosSim <- function(clustering,distMatrix,parallel){
  #find max distance and then do minCosSim <- 1 - maxDist
  clusters <- unique(clustering)
  clusters <- clusters[order(clusters)]
  res <- c()
  if (parallel==TRUE){
    res <- foreach(i=clusters) %dopar%{
      samplesInCluster <- names(clustering[clustering==i])
      nSamples <- length(samplesInCluster)
      maxDist <- 0
      if(nSamples==1){
        #maxDist <- 0
      }else{
        for(j in 1:(nSamples-1)){
          maxDist <- max(maxDist,max(distMatrix[samplesInCluster[j],samplesInCluster[1:nSamples > j]]))
        }
      }
      (1 - maxDist)
    }
    res <- unlist(res)
  }else{
    for(i in clusters){
      samplesInCluster <- names(clustering[clustering==i])
      nSamples <- length(samplesInCluster)
      maxDist <- 0
      if(nSamples==1){
        #maxDist <- 0
      }else{
        for(j in 1:(nSamples-1)){
          maxDist <- max(maxDist,max(distMatrix[samplesInCluster[j],samplesInCluster[1:nSamples > j]]))
        }
      }
      res <- c(res,1 - maxDist)
    }
  }
  return(res)
}

#function to compute the maximum between clusters cosine similarity,
#uses as input the distance matrix constructed as 1 - cosSimMatrix
maxBetweenClustersCosSim <- function(clustering,distMatrix,parallel){
  clusters <- unique(clustering)
  clusters <- clusters[order(clusters)]
  res <- matrix(1,nrow = length(clusters),ncol = length(clusters))
  colnames(res) <- clusters
  rownames(res) <- clusters
  if (parallel==TRUE){
    res_list <- foreach(c=1:(length(clusters) - 1)) %dopar%{
      i <- clusters[c]
      samplesInCluster <- names(clustering[clustering==i])
      res_table <- data.frame()
      for (d in (c+1):length(clusters)){
        j <- clusters[d]
        samplesInCluster2 <- names(clustering[clustering==j])
        ms <- 1 - min(distMatrix[samplesInCluster,samplesInCluster2])
        res_table <- rbind(res_table,data.frame(r=i,c=j,maxSim=ms))
      }
      res_table
    }
    
    for (i in 1:length(res_list)){
      currentTable <- res_list[[i]]
      for(j in 1:nrow(currentTable)){
        res[currentTable$r[j],currentTable$c[j]] <- currentTable$maxSim[j]
        res[currentTable$c[j],currentTable$r[j]] <- currentTable$maxSim[j]
      }
    }
  }else{
    for (c in 1:(length(clusters) - 1)){
      i <- clusters[c]
      samplesInCluster <- names(clustering[clustering==i])
      for (d in (c+1):length(clusters)){
        j <- clusters[d]
        samplesInCluster2 <- names(clustering[clustering==j])
        ms <- 1 - min(distMatrix[samplesInCluster,samplesInCluster2])
        res[i,j] <- ms
        res[j,i] <- ms
      }
    }
  }
  return(res)
}

plotHierarchicalCluster <- function(fit_clust,outFilePath,group,ns,nboots){
  output_file <- paste0(outFilePath,"Sigs_Cluster_Average_",group,"_ns",ns,"_nboots",nboots,".jpg")
  jpeg(output_file,width = 12*nboots*ns,height = 500,res = 80)
  plot(fit_clust)
  abline(a=0.1,b=0,col="red")
  dev.off()
}


plotWithinClusterCosSim <- function(cosSimHClust,cosSimPAM,cosSimMC,outFilePath,group,ns,nboots){
  dists <- c(cosSimHClust,cosSimPAM,cosSimMC)
  clustermethod <- c(rep("hclust",length(cosSimHClust)),rep("pam",length(cosSimPAM)),rep("MC",length(cosSimMC)))
  output_file <- paste0(outFilePath,"Sigs_WithinClusterCosSim_",group,"_ns",ns,"_nboots",nboots,".jpg")
  jpeg(output_file,width = 600,height = 500,res = 100)
  boxplot(dists ~ clustermethod, lwd = 2, ylab = 'mean Cosine Similarity',xlab = 'method',
          main = paste0("Within Cluster Cosine Similarity\n",group,", nSig=",ns))
  stripchart( dists ~ clustermethod, vertical = TRUE, 
             method = "jitter", add = TRUE, pch = 20, col = 'blue')
  dev.off()
}

plotWithinClusterSilWidth <- function(sil_hclust,sil_pam,sil_MC,outFilePath,group,ns,nboots){
  dists <- c(sil_hclust$clus.avg.widths,sil_pam$clus.avg.widths,sil_MC$clus.avg.widths)
  clustermethod <- c(rep("hclust",length(sil_hclust$clus.avg.widths)),rep("pam",length(sil_pam$clus.avg.widths)),rep("MC",length(sil_pam$clus.avg.widths)))
  output_file <- paste0(outFilePath,"Sigs_WithinClusterSilWidth_",group,"_ns",ns,"_nboots",nboots,".jpg")
  jpeg(output_file,width = 600,height = 500,res = 100)
  boxplot(dists ~ clustermethod, lwd = 2, ylab = 'mean Silhouette Width',xlab = 'method',
          main = paste0("Within Cluster Silhouette Width\n",group,", nSig=",ns))
  stripchart( dists ~ clustermethod, vertical = TRUE, 
              method = "jitter", add = TRUE, pch = 20, col = 'blue')
  dev.off()
}


plotOverallMetrics <- function(overall_metrics,whattoplot,overall_metrics_file,group,nboots,nmfmethod){
  jpeg(overall_metrics_file,width = 800,height = 500,res = 100)
  par(mar = c(5, 4, 4, 12),mgp = c(2.5,1,0))
  if(nmfmethod=="lee"){
    max_error <- max(overall_metrics$ave.RMSE,overall_metrics$ave.RMSE.orig)
    plot(overall_metrics$nsig,overall_metrics$ave.RMSE/max_error,type="l",
         ylab = "Score",xlab = "Number of extracted signatures",lwd = 2,
         ylim = c(min(min(overall_metrics$ave.RMSE/max_error,overall_metrics$ave.RMSE.orig/max_error),
                      min(overall_metrics[,whattoplot],na.rm = TRUE))*0.95,1),
         main = paste0("Overall Metrics\n(",group,", bootstraps=",nboots,")"))
    lines(overall_metrics$nsig,overall_metrics$ave.RMSE.orig/max_error,type="l",col="black",lwd = 2,lty = 2)
  }else{
    max_error <- max(overall_metrics$ave.KLD,overall_metrics$ave.KLD.orig)
    plot(overall_metrics$nsig,overall_metrics$ave.KLD/max_error,type="l",
         ylab = "Score",xlab = "Number of extracted signatures",lwd = 2,
         ylim = c(min(min(overall_metrics$ave.KLD/max_error,overall_metrics$ave.KLD.orig/max_error),
                      min(overall_metrics[,whattoplot],na.rm = TRUE))*0.95,1),
         main = paste0("Overall Metrics\n(",group,", bootstraps=",nboots,")"))
    lines(overall_metrics$nsig,overall_metrics$ave.KLD.orig/max_error,type="l",col="black",lwd = 2, lty = 2)
  }
  abline(h = 0.9,lty = 2)
  colours_list <- c("red","green","blue","purple","orange","brown","yellow")
  for(i in 1:length(whattoplot)){
    pos <- whattoplot[i]
    lines(overall_metrics$nsig,overall_metrics[,pos],type="l",col=colours_list[i],lwd = 2)
  }
  # lines(overall_metrics$nsig,overall_metrics$proportion.tooSimilar.Signatures,type="l",col="brown",lwd = 2)
  # legend("right", c("norm.Error","ave.CosSim.hclust","ave.CosSim.PAM","ave.SilWid.hclust","ave.SilWid.PAM","cophenetic.corr.hclust","prop.tooSimilar.Sig"), xpd = TRUE, horiz = FALSE, inset = c(-0.45,-0.4),lty = rep(1,7),
  #        bty = "n", col = c("black","red","green","blue","purple","orange","brown"),lwd = 2, cex = 0.9)
  legend("right", c("norm.Error","norm.Error (orig. cat.)",colnames(overall_metrics)[whattoplot]), xpd = TRUE, horiz = FALSE, inset = c(-0.45,-0.4),lty = c(1,2,rep(1,length(whattoplot))),
         bty = "n", col = c("black","black",colours_list[1:length(whattoplot)]),lwd = 2, cex = 0.9)
  dev.off()
}

findMedoidsHclust <- function(distMatrix,cut_res){
  clusters <- unique(cut_res)
  clusters <- clusters[order(clusters)]
  medoids_hclust <- c()
  for (j in clusters){
    current_cluster <- cut_res[cut_res==j]
    current_cluster_names <- names(current_cluster)
    nInCluster <- length(current_cluster)
    if (nInCluster==1){
      medoids_hclust <- c(medoids_hclust,names(current_cluster))
    }else{
      for (p in 1:nInCluster){
        current_cluster[current_cluster_names[p]] <- sum(distMatrix[current_cluster_names[p],current_cluster_names[-p]])/(nInCluster-1)
      }
      medoids_hclust <- c(medoids_hclust,names(which.min(current_cluster)[1]))
    }
    #message(names(which.min(current_cluster)))
    #message(nInCluster)
    #message(medoids_hclust)
  }
  return(medoids_hclust)
}

average_medoids_cosSim <- function(distMatrix,medoids){
  #at least two clusters!
  nClusters <- length(medoids)
  medoidsNumber <- 1:nClusters
  nCombinations <- nClusters*(nClusters - 1)/2
  amcs <- 0
  for (j in (medoidsNumber-1)){
    amcs <- amcs + sum(distMatrix[medoids[j],medoids[medoidsNumber>j]])
  }
  amcs <- (nCombinations - amcs)/nCombinations
  return(amcs)
}

medoids_cosSimMatrix <- function(p_boot,medoids){
  return(computeCorrelation_parallel(p_boot[,medoids]))
}

computePropTooSimilar <- function(distMatrix,saved_nmf_runs,ns){
  countTooSimilar <- 0
  for (i in 1:length(saved_nmf_runs)){
    if (max(distMatrix[((i-1)*ns + 1):(ns*i),((i-1)*ns + 1):(ns*i)]) > 0.9) countTooSimilar <- countTooSimilar + 1
  }
  return(countTooSimilar/saved_nmf_runs)
}

plotGenericSignatures <- function(signature_data_matrix,output_file = NULL,plot_sum = TRUE,overall_title = "",mar=NULL){
  # rearr.colours <- c(rep("blue",16),rep("black",16),rep("red",16),rep("grey",16),rep("green",16),rep("pink",16))
  nplotrows <- ceiling(ncol(signature_data_matrix)/3)
  if(!is.null(output_file)) {
    jpeg(output_file,width = 3*800,height = nplotrows*400,res = 190)
    par(mfrow = c(nplotrows, 3),oma=c(0,0,2,0))
  }
  for (pos in 1:ncol(signature_data_matrix)){
    par(mgp=c(1,1,0))
    if(is.null(mar)){
      par(mar=c(3,3,2,2))
    }else{
      par(mar=mar)
    }
    title <- colnames(signature_data_matrix)[pos]
    if (plot_sum) title <- paste0(title," (",round(sum(signature_data_matrix[,pos]))," substitutions)")
    # xlabels <- rep("",96)
    # xlabels[8] <- "C > A"
    # xlabels[24] <- "C > G"
    # xlabels[40] <- "C > T"
    # xlabels[56] <- "T > A"
    # xlabels[72] <- "T > C"
    # xlabels[88] <- "T > G"
    barplot(signature_data_matrix[,pos],
            main = title,
            #names.arg = row.names(signature_data_matrix),
            names.arg = c(rep("",nrow(signature_data_matrix))),
            col="skyblue",
            beside = TRUE,
            xlab = "channels",
            cex.names = 1)
  }
  title(main = overall_title,outer = TRUE,cex.main = 2)
  if(!is.null(output_file)) dev.off()
}

plotGenericSignatures_withMeanSd <- function(signature_data_matrix,mean_matrix,sd_matrix,output_file = NULL,plot_sum = TRUE,overall_title = "",mar=NULL){
  # rearr.colours <- c(rep("blue",16),rep("black",16),rep("red",16),rep("grey",16),rep("green",16),rep("pink",16))
  nplotrows <- ncol(signature_data_matrix)
  if(!is.null(output_file)) jpeg(output_file,width = 2*800,height = nplotrows*400,res = 190)
  par(mfrow = c(nplotrows, 2),oma=c(0,0,2,0))
  par(mgp=c(1,1,0))
  if(is.null(mar)){
    par(mar=c(3,3,2,2))
  }else{
    par(mar=mar)
  }
  for (pos in 1:ncol(signature_data_matrix)){
    ylimit <- c(0,max(signature_data_matrix[,pos],mean_matrix[,pos]+sd_matrix[,pos]))
    title <- colnames(signature_data_matrix)[pos]
    if (plot_sum) title <- paste0(title," (",round(sum(signature_data_matrix[,pos]))," substitutions)")
    # xlabels <- rep("",96)
    # xlabels[8] <- "C > A"
    # xlabels[24] <- "C > G"
    # xlabels[40] <- "C > T"
    # xlabels[56] <- "T > A"
    # xlabels[72] <- "T > C"
    # xlabels[88] <- "T > G"
    barplot(signature_data_matrix[,pos],
            main = title,
            #names.arg = row.names(signature_data_matrix),
            names.arg = c(rep("",nrow(signature_data_matrix))),
            col="skyblue",
            beside = TRUE,
            xlab = "channels",
            ylim = ylimit,
            cex.names = 1)
    barCenters <- barplot(mean_matrix[,pos],
            main = "mean and sd of cluster",
            #names.arg = row.names(signature_data_matrix),
            names.arg = c(rep("",nrow(signature_data_matrix))),
            col="skyblue",
            beside = TRUE,
            xlab = "channels",
            ylim = ylimit,
            cex.names = 1)
    segments(barCenters, mean_matrix[,pos] - sd_matrix[,pos], barCenters,
             mean_matrix[,pos] + sd_matrix[,pos], lwd = 1.5)
    
    arrows(barCenters, mean_matrix[,pos] - sd_matrix[,pos], barCenters,
           mean_matrix[,pos] + sd_matrix[,pos], lwd = 1.5, angle = 90,
           code = 3, length = 0.05)
  }
  title(main = overall_title,outer = TRUE,cex.main = 2)
  if(!is.null(output_file)) dev.off()
}

plotSubsSignatures <- function(signature_data_matrix,output_file = NULL,plot_sum = TRUE,overall_title = "",add_to_titles = NULL,mar=NULL){
  colnames(signature_data_matrix) <- sapply(colnames(signature_data_matrix),function(x) if (nchar(x)>20) paste0(substr(x,1,17),"...") else x)
  rearr.colours <- c(rep("blue",16),rep("black",16),rep("red",16),rep("grey",16),rep("green",16),rep("pink",16))
  nplotrows <- ceiling(ncol(signature_data_matrix)/3)
  if(!is.null(output_file)) {
    jpeg(output_file,width = 3*800,height = nplotrows*320,res = 220)
    par(mfrow = c(nplotrows, 3),oma=c(0,0,2,0))
  }
  for (pos in 1:ncol(signature_data_matrix)){
    if(is.null(mar)){
      par(mar=c(2,3,2,2))
    }else{
      par(mar=mar)
    }
    title <- colnames(signature_data_matrix)[pos]
    if (plot_sum) title <- paste0(title," (",round(sum(signature_data_matrix[,pos]))," substitutions)")
    if (!is.null(add_to_titles)) title <- paste0(title,"\n",add_to_titles[pos])
    muttypes <- c("C>A","C>G","C>T","T>A","T>C","T>G")
    xlabels <- rep("",96)
    # xlabels[8] <- "C > A"
    # xlabels[24] <- "C > G"
    # xlabels[40] <- "C > T"
    # xlabels[56] <- "T > A"
    # xlabels[72] <- "T > C"
    # xlabels[88] <- "T > G"
    barplot(signature_data_matrix[,pos],
            main = title,
            #names.arg = row.names(signature_data_matrix),
            names.arg = xlabels,
            col=rearr.colours,
            beside = TRUE,
            las=2,
            cex.names = 1,border = NA,space = 0.2)
    par(xpd=TRUE)
    par(usr = c(0, 1, 0, 1))
    recttop <- -0.02
    rectbottom <- -0.16
    start1 <- 0.035
    gap <- 0.155
    rect(start1, rectbottom, start1+gap, recttop,col = "blue",lwd = 0)
    rect(start1+gap, rectbottom, start1+2*gap, recttop,col = "black",lwd = 0)
    rect(start1+2*gap, rectbottom, start1+3*gap, recttop,col = "red",lwd = 0)
    rect(start1+3*gap, rectbottom, start1+4*gap, recttop,col = "grey",lwd = 0)
    rect(start1+4*gap, rectbottom, start1+5*gap, recttop,col = "green",lwd = 0)
    rect(start1+5*gap, rectbottom, start1+6*gap, recttop,col = "pink",lwd = 0)
    textposx <- 0.04+seq(8,88,16)/104
    text(x = textposx[1:3],y = -0.09,labels = muttypes[1:3],col = "white",font = 2)
    text(x = textposx[4:6],y = -0.09,labels = muttypes[4:6],col = "black",font = 2)
    #shadowtext(x = 0.04+seq(8,88,16)/104,y = rep(-0.09,6),labels = muttypes,col = "white",bg = "black",r=0.2)
    par(xpd=FALSE)
  }
  title(main = overall_title,outer = TRUE,cex.main = 2)
  if(!is.null(output_file)) dev.off()
}

plotSubsSignatures_withMeanSd <- function(signature_data_matrix,mean_matrix,sd_matrix,output_file = NULL,plot_sum = TRUE,overall_title = "",add_to_titles = NULL,mar=NULL){
  colnames(signature_data_matrix) <- sapply(colnames(signature_data_matrix),function(x) if (nchar(x)>20) paste0(substr(x,1,17),"...") else x)
  rearr.colours <- c(rep("blue",16),rep("black",16),rep("red",16),rep("grey",16),rep("green",16),rep("pink",16))
  nplotrows <- ncol(signature_data_matrix)
  if(!is.null(output_file)) jpeg(output_file,width = 2*800,height = nplotrows*320,res = 220)
  par(mfrow = c(nplotrows, 2),oma=c(0,0,2,0))
  if(is.null(mar)){
    par(mar=c(2,3,2,2))
  }else{
    par(mar=mar)
  }
  muttypes <- c("C>A","C>G","C>T","T>A","T>C","T>G")
  xlabels <- rep("",96)
  # xlabels[8] <- "C > A"
  # xlabels[24] <- "C > G"
  # xlabels[40] <- "C > T"
  # xlabels[56] <- "T > A"
  # xlabels[72] <- "T > C"
  # xlabels[88] <- "T > G"
  for (pos in 1:ncol(signature_data_matrix)){
    ylimit <- c(0,max(signature_data_matrix[,pos],mean_matrix[,pos]+sd_matrix[,pos]))
    title <- colnames(signature_data_matrix)[pos]
    if (plot_sum) title <- paste0(title," (",round(sum(signature_data_matrix[,pos]))," substitutions)")
    if (!is.null(add_to_titles)) title <- paste0(title,"\n",add_to_titles[pos])
    barplot(signature_data_matrix[,pos],
            main = title,
            #names.arg = row.names(signature_data_matrix),
            names.arg = xlabels,
            col=rearr.colours,
            beside = TRUE,
            ylim = ylimit,
            las=2,
            cex.names = 1)
    par(xpd=TRUE)
    par(usr = c(0, 1, 0, 1))
    recttop <- -0.02
    rectbottom <- -0.16
    start1 <- 0.035
    gap <- 0.155
    rect(start1, rectbottom, start1+gap, recttop,col = "blue",lwd = 0)
    rect(start1+gap, rectbottom, start1+2*gap, recttop,col = "black",lwd = 0)
    rect(start1+2*gap, rectbottom, start1+3*gap, recttop,col = "red",lwd = 0)
    rect(start1+3*gap, rectbottom, start1+4*gap, recttop,col = "grey",lwd = 0)
    rect(start1+4*gap, rectbottom, start1+5*gap, recttop,col = "green",lwd = 0)
    rect(start1+5*gap, rectbottom, start1+6*gap, recttop,col = "pink",lwd = 0)
    textposx <- 0.04+seq(8,88,16)/104
    text(x = textposx[1:3],y = -0.09,labels = muttypes[1:3],col = "white",font = 2)
    text(x = textposx[4:6],y = -0.09,labels = muttypes[4:6],col = "black",font = 2)
    par(xpd=FALSE)
    barCenters <- barplot(mean_matrix[,pos],
                          main = "mean and sd of cluster",
                          #names.arg = row.names(signature_data_matrix),
                          names.arg = xlabels,
                          col=rearr.colours,
                          #border = NA,
                          beside = TRUE,
                          ylim = ylimit,
                          las=2,
                          cex.names = 1,border = NA,space = 0.2)
    # segments(barCenters, mean_matrix[,pos] - sd_matrix[,pos], barCenters,
    #          mean_matrix[,pos] + sd_matrix[,pos], lwd = 1.5)
    segments(barCenters, mean_matrix[,pos], barCenters,
             mean_matrix[,pos] + sd_matrix[,pos], lwd = 1)
    # arrows(barCenters, mean_matrix[,pos] - sd_matrix[,pos], barCenters,
    #        mean_matrix[,pos] + sd_matrix[,pos], lwd = 1.5, angle = 90,
    #        code = 3, length = 0.05)
    par(xpd=TRUE)
    par(usr = c(0, 1, 0, 1))
    recttop <- -0.02
    rectbottom <- -0.16
    start1 <- 0.035
    gap <- 0.155
    rect(start1, rectbottom, start1+gap, recttop,col = "blue",lwd = 0)
    rect(start1+gap, rectbottom, start1+2*gap, recttop,col = "black",lwd = 0)
    rect(start1+2*gap, rectbottom, start1+3*gap, recttop,col = "red",lwd = 0)
    rect(start1+3*gap, rectbottom, start1+4*gap, recttop,col = "grey",lwd = 0)
    rect(start1+4*gap, rectbottom, start1+5*gap, recttop,col = "green",lwd = 0)
    rect(start1+5*gap, rectbottom, start1+6*gap, recttop,col = "pink",lwd = 0)
    textposx <- 0.04+seq(8,88,16)/104
    text(x = textposx[1:3],y = -0.09,labels = muttypes[1:3],col = "white",font = 2)
    text(x = textposx[4:6],y = -0.09,labels = muttypes[4:6],col = "black",font = 2)
    par(xpd=FALSE)
  }
  title(main = overall_title,outer = TRUE,cex.main = 2)
  if(!is.null(output_file)) dev.off()
}



plotRearrSignatures <-function(signature_data_matrix,output_file = NULL,plot_sum = TRUE,overall_title = "",add_to_titles = NULL,mar=NULL){
  #This function plots a set of signatures in a single file, three signatures for each row.
  #signature_data_matrix is a data frame that contains the rearrangement signatures.
  #                      The columns are the signatures, while the rows are the 32 features
  colnames(signature_data_matrix) <- sapply(colnames(signature_data_matrix),function(x) if (nchar(x)>20) paste0(substr(x,1,17),"...") else x)
  del_col = rgb(228,26,28, maxColorValue = 255)
  td_col = rgb(77,175,74, maxColorValue =255)
  inv_col  = rgb(55,126,184, maxColorValue = 255)
  transloc_col = rgb(152,78,163, maxColorValue =255)
  non_clust_col = rgb(240,240,240, maxColorValue =255)
  #rearr.colours <- c(rep("darkblue",16),rep("red",16))
  rearr.colours <- rep(c(rep(del_col,5),rep(td_col,5),rep(inv_col,5),transloc_col),2)
  nplotrows <- ceiling(ncol(signature_data_matrix)/3)
  if(!is.null(output_file)){
    jpeg(output_file,width = 3*800,height = nplotrows*500,res = 220)
    par(mfrow = c(nplotrows, 3),oma=c(0,0,2,0))
  }
  sizes <- c("1-10Kb",
             "10-100Kb",
             "100Kb-1Mb",
             "1Mb-10Mb",
             ">10Mb")
  sizes_names <- c(rep(sizes,3),"",rep(sizes,3),"")
  for (pos in 1:ncol(signature_data_matrix)){
    if(is.null(mar)){
      par(mar=c(8,3,2,2))
    }else{
      par(mar=mar)
    }
    title <- colnames(signature_data_matrix)[pos]
    if (plot_sum) title <- paste0(title," (",sum(signature_data_matrix[,pos])," rearrangements)")
    if (!is.null(add_to_titles)) title <- paste0(title,"\n",add_to_titles[pos])
    pos <- barplot(signature_data_matrix[,pos],
                   main = title,
                   names.arg = NA,
                   #names.arg = sizes_names,
                   col=rearr.colours,
                   beside = FALSE,
                   #las=2,
                   cex.names = 0.8,
                   #mgp=c(3,2,0),
                   border = 0,
                   space = 0.1)
    axis(1,
         las=2,
         #hadj=0.5,
         at=pos,
         lab=sizes_names,
         #mgp=c(3,2,0),
         col = "transparent",
         line = 1,
         cex.axis = 0.8)
    #save old plot coordinates
    op <- par("usr")
    #set new coordinates
    par(usr = c(0, 1, 0, 1))
    #add graphics
    par(xpd=TRUE)
    start1 <- 0.035
    xsep = 0.145
    start1_text <- 0.11
    tr_size <- 0.03
    #rect(0.1, 0.1, 0.2, 0.2,col = "blue",lwd = 0)
    #rect(0.1, 0.1, 0.11, 0.11,col = "red",lwd = 0)
    stop <- start1
    for(i in 1:2){
      start <- stop
      stop <- start + xsep
      rect(start, -0.14, stop, -0.02,col = del_col,lwd = 0)
      text(x = start+0.5*xsep,y = -0.08,"del",col = "white")
      start <- stop
      stop <- start + xsep
      rect(start, -0.14, stop, -0.02,col = td_col,lwd = 0)
      text(x = start+0.5*xsep,y = -0.08,"tds",col = "white")
      start <- stop
      stop <- start + xsep
      rect(start, -0.14, stop, -0.02,col = inv_col,lwd = 0)
      text(x = start+0.5*xsep,y = -0.08,"inv",col = "white")
      start <- stop
      stop <- start + tr_size
      rect(start, -0.14, stop, -0.02,col = transloc_col,lwd = 0)
      text(x = start+0.5*tr_size,y = -0.08,"tr",col = "white")
    }
    xsep2 <- 3*xsep+tr_size
    rect(start1, -0.26, start1+xsep2, -0.14,col = "black",lwd = 0)
    text(x = start1+0.5*xsep2,y = -0.2,"clustered",col = "white")
    rect(start1+xsep2, -0.26, start1+2*xsep2, -0.14,col = non_clust_col,lwd = 0)
    text(x = start1+1.5*xsep2,y = -0.2,"non-clustered",col = "black")
    
    #restore old coordinates
    par(usr = op)
  }
  title(main = overall_title,outer = TRUE,cex.main = 2)
  if(!is.null(output_file)) dev.off()
}


plotRearrSignatures_withMeanSd <-function(signature_data_matrix,mean_matrix,sd_matrix,output_file = NULL,plot_sum = TRUE,overall_title = "",add_to_titles = NULL,mar=NULL){
  #This function plots a set of signatures in a single file, three signatures for each row.
  #signature_data_matrix is a data frame that contains the rearrangement signatures.
  #                      The columns are the signatures, while the rows are the 32 features
  colnames(signature_data_matrix) <- sapply(colnames(signature_data_matrix),function(x) if (nchar(x)>20) paste0(substr(x,1,17),"...") else x)
  del_col = rgb(228,26,28, maxColorValue = 255)
  td_col = rgb(77,175,74, maxColorValue =255)
  inv_col  = rgb(55,126,184, maxColorValue = 255)
  transloc_col = rgb(152,78,163, maxColorValue =255)
  non_clust_col = rgb(240,240,240, maxColorValue =255)
  #rearr.colours <- c(rep("darkblue",16),rep("red",16))
  rearr.colours <- rep(c(rep(del_col,5),rep(td_col,5),rep(inv_col,5),transloc_col),2)
  nplotrows <- ncol(signature_data_matrix)
  if(!is.null(output_file)) jpeg(output_file,width = 2*800,height = nplotrows*500,res = 220)
  par(mfrow = c(nplotrows, 2),oma=c(0,0,2,0))
  sizes <- c("1-10Kb",
             "10-100Kb",
             "100Kb-1Mb",
             "1Mb-10Mb",
             ">10Mb")
  sizes_names <- c(rep(sizes,3),"",rep(sizes,3),"")
  
  rearrAxis <- function(barCenters,sizes_names){
    axis(1,
         las=2,
         #hadj=0.5,
         at=barCenters,
         lab=sizes_names,
         #mgp=c(3,2,0),
         col = "transparent",
         line = 1,
         cex.axis = 0.8)
    #save old plot coordinates
    op <- par("usr")
    #set new coordinates
    par(usr = c(0, 1, 0, 1))
    #add graphics
    par(xpd=TRUE)
    start1 <- 0.035
    xsep = 0.145
    start1_text <- 0.11
    tr_size <- 0.03
    #rect(0.1, 0.1, 0.2, 0.2,col = "blue",lwd = 0)
    #rect(0.1, 0.1, 0.11, 0.11,col = "red",lwd = 0)
    stop <- start1
    for(i in 1:2){
      start <- stop
      stop <- start + xsep
      rect(start, -0.14, stop, -0.02,col = del_col,lwd = 0)
      text(x = start+0.5*xsep,y = -0.08,"del",col = "white")
      start <- stop
      stop <- start + xsep
      rect(start, -0.14, stop, -0.02,col = td_col,lwd = 0)
      text(x = start+0.5*xsep,y = -0.08,"tds",col = "white")
      start <- stop
      stop <- start + xsep
      rect(start, -0.14, stop, -0.02,col = inv_col,lwd = 0)
      text(x = start+0.5*xsep,y = -0.08,"inv",col = "white")
      start <- stop
      stop <- start + tr_size
      rect(start, -0.14, stop, -0.02,col = transloc_col,lwd = 0)
      text(x = start+0.5*tr_size,y = -0.08,"tr",col = "white")
    }
    xsep2 <- 3*xsep+tr_size
    rect(start1, -0.26, start1+xsep2, -0.14,col = "black",lwd = 0)
    text(x = start1+0.5*xsep2,y = -0.2,"clustered",col = "white")
    rect(start1+xsep2, -0.26, start1+2*xsep2, -0.14,col = non_clust_col,lwd = 0)
    text(x = start1+1.5*xsep2,y = -0.2,"non-clustered",col = "black")
    
    #restore old coordinates
    par(usr = op)
  }
  
  for (pos in 1:ncol(signature_data_matrix)){
    ylimit <- c(0,max(signature_data_matrix[,pos],mean_matrix[,pos]+sd_matrix[,pos]))
    if(is.null(mar)){
      par(mar=c(8,3,2,2))
    }else{
      par(mar=mar)
    }
    title <- colnames(signature_data_matrix)[pos]
    if (plot_sum) title <- paste0(title," (",sum(signature_data_matrix[,pos])," rearrangements)")
    if (!is.null(add_to_titles)) title <- paste0(title,"\n",add_to_titles[pos])
    barCenters <- barplot(signature_data_matrix[,pos],
                   main = title,
                   names.arg = NA,
                   #names.arg = sizes_names,
                   col=rearr.colours,
                   beside = FALSE,
                   #las=2,
                   cex.names = 0.8,
                   #mgp=c(3,2,0),
                   border = 0,
                   ylim = ylimit,
                   space = 0.1)
    rearrAxis(barCenters,sizes_names)
    barCenters <- barplot(mean_matrix[,pos],
                          main = "mean and sd of cluster",
                          #names.arg = row.names(signature_data_matrix),
                          names.arg = NA,
                          col=rearr.colours,
                          #border = NA,
                          beside = FALSE,
                          ylim = ylimit,
                          las=2,
                          border = 0,
                          space = 0.1)
    # segments(barCenters, mean_matrix[,pos] - sd_matrix[,pos], barCenters,
    #          mean_matrix[,pos] + sd_matrix[,pos], lwd = 1.5)
    segments(barCenters, mean_matrix[,pos], barCenters,
             mean_matrix[,pos] + sd_matrix[,pos], lwd = 1.5)
    # arrows(barCenters, mean_matrix[,pos] - sd_matrix[,pos], barCenters,
    #        mean_matrix[,pos] + sd_matrix[,pos], lwd = 1.5, angle = 90,
    #        code = 3, length = 0.05)
    rearrAxis(barCenters,sizes_names)
  }
  title(main = overall_title,outer = TRUE,cex.main = 2)
  if(!is.null(output_file)) dev.off()
}


plot.CosSimMatrix <- function(CosSimMatrix,output_file,dpi=300,xlabel = "",ylabel = "",thresholdMark = 0.9,extraWidth = 500,extraHeight = 500){
  library("ggplot2")
  
  # Set up the vectors                           
  signatures.names <- colnames(CosSimMatrix)
  sample.names <- row.names(CosSimMatrix)
  
  # Create the data frame
  df <- expand.grid(sample.names,signatures.names)
  df$value <- unlist(CosSimMatrix)   
  df$labels <- sprintf("%.2f", df$value)
  df$labels[df$value==0] <- ""
  
  #Plot the Data (500+150*nsamples)x1200
  g <- ggplot(df, aes(Var1, Var2)) + geom_point(aes(size = value, colour = value>thresholdMark)) + theme_bw() + xlab(xlabel) + ylab(ylabel)
  g <- g + scale_size_continuous(range=c(0,10)) + geom_text(aes(label = labels))
  g + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=14),
            axis.text.y = element_text(vjust = 1, size=14)) + theme(legend.position="none")
  w <- (extraWidth+150*length(sample.names))/dpi
  h <- (extraHeight+150*length(signatures.names))/dpi
  ggsave(filename = output_file,dpi = dpi,height = h,width = w)
}


plot.CosSimSignatures <- function(sig1,sig2,output_file,dpi=300,xlabel = "",ylabel = ""){
  cos.sim <- function(a, b){
    return( sum(a*b)/sqrt(sum(a^2)*sum(b^2)) )
  }  
  cos_sim_df <- data.frame()
  for (s in colnames(sig1)){
    for(a in colnames(sig2)){
      cos_sim_df[s,a] <- cos.sim(sig1[,s],sig2[,a])
    }
  }
  plot.CosSimMatrix(cos_sim_df,output_file,dpi=dpi,xlabel = xlabel,ylabel = ylabel)
}

#returns the list of signatures identified. For example,
#c("C1","C3","N1","C13","N2")
#means that Cosmic (C) signatures 1, 3 and 13 were found, while
#signatures N1 and N2 are unknown signatures (N for not found), based on the fact
#that no similarity >=threshold was found with the Cosmic 30
findClosestCOSMIC30 <- function(sigs,threshold){
  #load COSMIC30
  cosmic30 <- read.table("../data/COSMIC_signatures.txt", sep="\t", header=T, as.is=T, check.names = FALSE)
  #compute cos sim matrix
  cos_sim_df <- data.frame()
  for (s in colnames(sigs)){
    for(a in colnames(cosmic30)){
      cos_sim_df[s,a] <- cos.sim(sigs[,s],cosmic30[,a])
    }
  }
  max.sim <- apply(cos_sim_df,1,max)
  closestCosmic <- apply(cos_sim_df,1,which.max)
  notFound <- max.sim < threshold
  closestCosmic[notFound] <- paste0("N",1:sum(notFound))
  closestCosmic[!notFound] <- paste0("C",closestCosmic[!notFound])
  return(closestCosmic)
}

#automatically detect similarity with sum of two COSMIC30
findClosestCOSMIC30andCombinations <- function(sigs,threshold){
  #load COSMIC30
  cosmic30 <- read.table("../data/COSMIC_signatures.txt", sep="\t", header=T, as.is=T, check.names = FALSE)
  ncols30 <- length(cosmic30)
  colnames(cosmic30) <- 1:ncols30
  for(i in 1:(ncols30-1)){
    for(j in (i+1):ncols30){
      #message(paste0(i,"+",j))
      cosmic30[,paste0(i,"+",j)] <- (cosmic30[,i]+cosmic30[,j])/sum(cosmic30[,i]+cosmic30[,j])
    }
  }
  #compute cos sim matrix
  cos_sim_df <- data.frame()
  for (s in colnames(sigs)){
    for(a in colnames(cosmic30)){
      cos_sim_df[s,a] <- cos.sim(sigs[,s],cosmic30[,a])
    }
  }
  max.sim <- apply(cos_sim_df,1,max)
  closestCosmic <- colnames(cosmic30)[apply(cos_sim_df,1,which.max)]
  notFound <- max.sim < threshold
  closestCosmic[notFound] <- paste0("N",1:sum(notFound))
  closestCosmic[!notFound] <- paste0("C",closestCosmic[!notFound])
  return(closestCosmic)
}

#returns the list of signatures identified and the corresponding similarity. For example,
#res$cosmic = c("C1","C3","C13")
#res$cos.sim = c(0.94,0.85,0.7)
#means that Cosmic (C) signatures 1, 3 and 13 were found, while
#the corrsponding similarities to those signatures are 0.94, 0.85 and 0.7
findClosestCOSMIC30_withSimilarity <- function(sigs){
  #load COSMIC30
  cosmic30 <- read.table("../data/COSMIC_signatures.txt", sep="\t", header=T, as.is=T, check.names = FALSE)
  #compute cos sim matrix
  cos_sim_df <- data.frame()
  for (s in colnames(sigs)){
    for(a in colnames(cosmic30)){
      cos_sim_df[s,a] <- cos.sim(sigs[,s],cosmic30[,a])
    }
  }
  max.sim <- apply(cos_sim_df,1,max)
  closestCosmic <- apply(cos_sim_df,1,which.max)
  closestCosmic <- paste0("C",closestCosmic)
  res <- list()
  res[["cosmic"]] <- closestCosmic
  res[["cos.sim"]] <- max.sim
  return(res)
}

#automatically detect similarity with sum of two COSMIC30
findClosestCOSMIC30andCombinations_withSimilarity <- function(sigs){
  #load COSMIC30
  cosmic30 <- read.table("../data/COSMIC_signatures.txt", sep="\t", header=T, as.is=T, check.names = FALSE)
  ncols30 <- length(cosmic30)
  colnames(cosmic30) <- 1:ncols30
  for(i in 1:(ncols30-1)){
    for(j in (i+1):ncols30){
      #message(paste0(i,"+",j))
      cosmic30[,paste0(i,"+",j)] <- (cosmic30[,i]+cosmic30[,j])/sum(cosmic30[,i]+cosmic30[,j])
    }
  }
  #compute cos sim matrix
  cos_sim_df <- data.frame()
  for (s in colnames(sigs)){
    for(a in colnames(cosmic30)){
      cos_sim_df[s,a] <- cos.sim(sigs[,s],cosmic30[,a])
    }
  }
  max.sim <- apply(cos_sim_df,1,max)
  closestCosmic <- colnames(cosmic30)[apply(cos_sim_df,1,which.max)]
  closestCosmic <- paste0("C",closestCosmic)
  res <- list()
  res[["cosmic"]] <- closestCosmic
  res[["cos.sim"]] <- max.sim
  return(res)
}

#returns the list of signatures identified. For example,
#c("R1","R3","N1","R5","N2")
#means that Rearrangement (R) signatures 1, 3 and 5 were found, while
#signatures N1 and N2 are unknown signatures (N for not found), based on the fact
#that no similarity >=threshold was found with the Rearr Sigs from Breast 560 study
findClosestRearrSigsBreast560 <- function(sigs,threshold){
  #load RS.Breast560
  RS.Breast560 <- read.table("../data/rearrangement.signatures.txt", sep="\t", header=T, as.is=T, check.names = FALSE)
  #compute cos sim matrix
  cos_sim_df <- data.frame()
  for (s in colnames(sigs)){
    for(a in colnames(RS.Breast560)){
      cos_sim_df[s,a] <- cos.sim(sigs[,s],RS.Breast560[,a])
    }
  }
  max.sim <- apply(cos_sim_df,1,max)
  closestRS.Breast560 <- apply(cos_sim_df,1,which.max)
  notFound <- max.sim < threshold
  closestRS.Breast560[notFound] <- paste0("N",1:sum(notFound))
  closestRS.Breast560[!notFound] <- paste0("R",closestRS.Breast560[!notFound])
  return(closestRS.Breast560)
}

#returns the list of signatures identified and the corresponding similarity. For example,
#res$RS.Breast560 = c("R1","R3","R5")
#res$cos.sim = c(0.94,0.85,0.7)
#means that Rearrangement (R) signatures 1, 3 and 5 were found, while
#the corrsponding similarities to those signatures are 0.94, 0.85 and 0.7
findClosestRearrSigsBreast560_withSimilarity <- function(sigs){
  #load RS.Breast560
  RS.Breast560 <- read.table("../data/rearrangement.signatures.txt", sep="\t", header=T, as.is=T, check.names = FALSE)
  #compute cos sim matrix
  cos_sim_df <- data.frame()
  for (s in colnames(sigs)){
    for(a in colnames(RS.Breast560)){
      cos_sim_df[s,a] <- cos.sim(sigs[,s],RS.Breast560[,a])
    }
  }
  max.sim <- apply(cos_sim_df,1,max)
  closestRS.Breast560 <- apply(cos_sim_df,1,which.max)
  closestRS.Breast560 <- paste0("R",closestRS.Breast560)
  res <- list()
  res[["RS.Breast560"]] <- closestRS.Breast560
  res[["cos.sim"]] <- max.sim
  return(res)
}

KLD <- function(m1,m2){
  # print(sessionInfo())
  # print(m1)
  # print(m2)
  # m1 <- as.vector(as.matrix(cat))
  # m2 <- as.vector(as.matrix(m2))
  m1[m1==0] <- .Machine$double.eps
  m2[m2==0] <- .Machine$double.eps
  return(sum(m1*(log(m1)-log(m2)) - m1 + m2))
}

#samples/sigantures are ararnged by columns
computeCorrelationOfTwoSetsOfSigs <- function(sigs1,sigs2){
  cos_sim_df <- data.frame()
  for (s in colnames(sigs1)){
    for(a in colnames(sigs2)){
      cos_sim_df[s,a] <- cos.sim(sigs1[,s],sigs2[,a])
    }
  }
  return(cos_sim_df)
}

removeSimilarCatalogueSamples <- function(cat,cosSimThreshold = 0.99){
  catCosSimCorr <- computeCorrelation(cat) - diag(ncol(cat))
  newcat <- data.frame(cat[,1],row.names = row.names(cat))
  colnames(newcat) = colnames(cat)[1]
  #as you are building the new catalogue, add samples only if they have cos sim
  #less than cosSimThreshold w.r.t. samples in the new catalgue
  for (j in 2:ncol(cat)){
    cossimres <- catCosSimCorr[j,colnames(newcat)]
    if(!any(cossimres>cosSimThreshold)){
      newcat <- cbind(newcat,cat[,j])
      colnames(newcat)[ncol(newcat)] <- colnames(cat)[j]
    }
  }
  return(newcat)
}

normaliseSamples <- function(cat){
  return(cat/matrix(rep(apply(cat,2,sum),nrow(cat)),nrow = nrow(cat),byrow = TRUE))
}

writeTable <- function(t,file){
  write.table(t,file = file,sep = "\t",quote = FALSE,row.names = TRUE,col.names = TRUE)
}
readTable <- function(file){
  read.table(file = file,sep = "\t",check.names = FALSE,header = TRUE,as.is = TRUE)
}

#############################################

# test below
test_library <- function(){
  setwd("~/sandbox/git/signature-tools/lib")
  dateOfAnalysis <- "2018_05_15"
  n_row <- 32
  n_col <- 50
  rnd_matrix <- round(matrix(runif(n_row*n_col,min = 0,max = 50),nrow = n_row,ncol = n_col))
  colnames(rnd_matrix) <- paste0("C",1:n_col)
  row.names(rnd_matrix) <- paste0("R",1:n_row)
  SignatureExtraction(cat = rnd_matrix,
                      outFilePath = paste0("../results/",dateOfAnalysis,"_tests/"),
                      nrepeats = 10,
                      nboots = 2,
                      nparallel = 2,
                      nsig = 2:3,
                      mut_thr = 0,
                      type_of_extraction = "rearr",
                      project = "test",
                      parallel = TRUE,
                      nmfmethod = "brunet")
}
