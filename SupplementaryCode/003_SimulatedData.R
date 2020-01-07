runLocally <- TRUE

if (runLocally){
  setwd("~/sandbox/git/signature-tools/tests")
  resultdir <- "../results/testsSimulatedData/"
  ngenomes <- 30
}else{
  args = commandArgs(trailingOnly=TRUE)
  ngenomes <- as.numeric(args[1])
  resultdir <- "~/scratch116/testsSimulatedData/"
}

source("../lib/SignatureExtractionLib.R")

n <- 5
nsig_range <- 8:12
nparallel <- 10
nmf_algorithms <- c("brunet","lee","nsNMF")
nmf_algorithms <- nmf_algorithms[1]

nrepeats_withFilter <- 50
nboots_withFilter <- 20
nrepeats_withoutFilter <- 1
nboots_withoutFilter <- 1000

#local options run for debug
if (runLocally){
  nmf_algorithms <- nmf_algorithms[1]
  nrepeats_withFilter <- 20
  nboots_withFilter <- 10
  nrepeats_withoutFilter <- 1
  nboots_withoutFilter <- 80
}

#Cosmic Signatures
cosmic30 <- read.table("../data/COSMIC_signatures.txt",sep = "\t",header = TRUE,as.is = TRUE,check.names = FALSE)
cosmic30 <- sortCatalogue(cosmic30)
signatures <- c(1,2,3,5,6,8,12,13,17,18)
P <- cosmic30[,signatures]

#parameters
nsignatures_per_sample <- 5
min_mutations <- 1000
max_mutations <- 50000

#change directory
resultdir <- paste0(resultdir,"sigPerSample",nsignatures_per_sample,"_samples",ngenomes,"/")

generateCatalogue <- function(P,ngenomes,nsignatures_per_sample,min_mutations,max_mutations){
  nmutations <- exp(runif(ngenomes,log(min_mutations),log(max_mutations)))
  nsignatures <- ncol(P)
  E <- matrix(runif(ngenomes*nsignatures),nrow = nsignatures)
  #make more sparse
  if (nsignatures>nsignatures_per_sample){
    for (i in 1:ncol(E)){
      selected_rows <- sample(1:nsignatures,nsignatures - nsignatures_per_sample)
      E[selected_rows,i] <- 0
    }
  }
  E <- E/matrix(rep(apply(E,2,sum),nsignatures),nrow = nsignatures,byrow = TRUE)
  E <- E*matrix(rep(nmutations,nsignatures),nrow = nsignatures,byrow = TRUE)
  
  #compose simulated catalogue
  catalogue <- as.matrix(P) %*% E
  catalogue <- round(catalogue)
  row.names(catalogue) <- row.names(P)
  colnames(catalogue) <- paste0("Sample_",1:ngenomes)
  #add poisson noise
  for (i in nrow(catalogue)){
    catalogue[i,] <- rpois(ngenomes,catalogue[i,])
  }
  
  return(list(catalogue=catalogue,E=E))
}

for (algo in nmf_algorithms){ 
  for (i in 1:n){
  
    #check if catalogue has already been generated
    cat_file <- paste0(resultdir,"simulated_catalogue_",i,".tsv")
    E_file <- paste0(resultdir,"simulated_E_",i,".tsv")
    P_file <- paste0(resultdir,"simulated_P_",i,".tsv")
    if (file.exists(cat_file)){
      catalogue <- read.table(file = cat_file,sep = "\t",header = TRUE,check.names = FALSE,as.is = TRUE)
      E <- read.table(file = E_file,sep = "\t",header = TRUE,check.names = FALSE,as.is = TRUE)
      P <- read.table(file = P_file,sep = "\t",header = TRUE,check.names = FALSE,as.is = TRUE)
      message("catalogue ",i," loaded from file")
    } else {
      dir.create(resultdir,recursive = TRUE,showWarnings = FALSE)
      res <- generateCatalogue(P,ngenomes,nsignatures_per_sample,min_mutations,max_mutations)
      catalogue <- res$catalogue
      E <- res$E
      write.table(catalogue,file = cat_file,sep = "\t",col.names = TRUE,row.names = TRUE,quote = FALSE)
      write.table(E,file = E_file,sep = "\t",col.names = TRUE,row.names = TRUE,quote = FALSE)
      write.table(P,file = P_file,sep = "\t",col.names = TRUE,row.names = TRUE,quote = FALSE)
    }

    message("Extraction on catalogue ",i," with ",algo," and filtering best solutions")
    SignatureExtraction(cat=catalogue,
                        outFilePath=paste0(resultdir,"with_filterBest_",algo,"_",i,"/"),
                        nrepeats=nrepeats_withFilter,
                        nboots=nboots_withFilter,
                        clusteringMethod = "MC",
                        nparallel=nparallel,
                        nsig=nsig_range,
                        mut_thr=0,
                        type_of_extraction = "subs",
                        project = "SimulatedData",
                        parallel=TRUE,
                        nmfmethod = algo,
                        normaliseCatalogue = FALSE,
                        removeDuplicatesInCatalogue = FALSE,
                        plotCatalogue = TRUE)

    message("Extraction on catalogue ",i," with ",algo," without filtering best solutions")    
    SignatureExtraction(cat=catalogue,
                        outFilePath=paste0(resultdir,"without_filterBest_",algo,"_",i,"/"),
                        nrepeats=nrepeats_withoutFilter,
                        nboots=nboots_withoutFilter,
                        clusteringMethod = "MC",
                        filterBestOfEachBootstrap=FALSE,
                        nparallel=nparallel,
                        nsig=nsig_range,
                        mut_thr=0,
                        type_of_extraction = "subs",
                        project = "SimulatedData",
                        parallel=TRUE,
                        nmfmethod = algo,
                        normaliseCatalogue = FALSE,
                        removeDuplicatesInCatalogue = FALSE,
                        plotCatalogue = TRUE)
  }
}


