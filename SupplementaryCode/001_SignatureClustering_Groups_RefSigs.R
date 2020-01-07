source("../lib/SignatureExtractionLib.R")
#library(signature.tools.lib)
library(NMF)

#functions----

subs_sigs_names <- c("1","18","17","MMR1","2",
                     "MMR2","5","8","13","3",
                     "30","MIXED1","4","N1","PLATINUM",
                     "MIXED2","33","22","10","36",
                     "7","MIXED3","16","19","N2",
                     "9","N3","52","11","38",
                     "51","N4","N5","N6","N7",
                     "N8","N9","N10","N11","24",
                     "N12")

rearr_sigs_names <- c("R2","R4","R6a","R1","R7",
                     "R5","R6b","R8","R3","R9",
                     "R10","R11","R12","R13","R14",
                     "R15","R16","R17","R18","R19",
                     "R20")

#change signature names
change_numbers_to_letters <- function(x){
  splitted_text <- base::strsplit(x,split = "_")
  result <- c()
  for (i in 1:length(splitted_text)) {
    current_text <- splitted_text[[i]]
    n <- length(current_text)
    organ <- paste(current_text[1:(n-1)],collapse = "_")
    sigletter <- LETTERS[as.numeric(current_text[n])]
    #hardcode for Skin and Colorectal to remove A/B extraction 
    if(current_text[1]=="SkinA" | current_text[1]=="SkinB") organ <- "Skin"
    if(current_text[1]=="ColorectalA" | current_text[1]=="ColorectalB") organ <- "Colorectal"
    if(current_text[1]=="SkinB") sigletter <- "J"
    if(current_text[1]=="ColorectalB" ) sigletter <- "I"
    result <- c(result,paste(organ,sigletter,sep = "."))
  }
  return(result)
}

findClosestCOSMIC_PCAWG_withSimilarity <- function(sigs){
  #load COSMIC_PCAWG
  cosmic30 <- read.table("../data/sigProfiler_SBS_signatures.csv", sep=",", header=T, as.is=T, check.names = FALSE)
  cosmic30 <- cosmic30[,3:ncol(cosmic30)]
  colnames(cosmic30) <- paste0("C",substr(colnames(cosmic30),4,6))
  #compute cos sim matrix
  cos_sim_df <- data.frame()
  for (s in colnames(sigs)){
    for(a in colnames(cosmic30)){
      cos_sim_df[s,a] <- cos.sim(sigs[,s],cosmic30[,a])
    }
  }
  max.sim <- apply(cos_sim_df,1,max)
  closestCosmic <- colnames(cosmic30)[apply(cos_sim_df,1,which.max)]
  res <- list()
  res[["cosmic"]] <- closestCosmic
  res[["cos.sim"]] <- max.sim
  return(res)
}

results_dir <- paste0("../results/")
out_dir <- paste0(results_dir,"pancan_extraction_clustering/")
dir.create(out_dir,showWarnings = FALSE,recursive = TRUE)
dateOfAnalysis <- format(Sys.time(), "%Y_%m_%d")
out_dir <- paste0(out_dir,dateOfAnalysis)
nboots <- 20

plot_clustering_of_pancan <- function(out_dir,partitions,hc_fit,t,where_to_cut=NULL,tag){
  #find out how many signatures are in each partition
  table_n_part <- table(partitions)
  #reorder
  new_partitions_order <- order(table_n_part,decreasing = TRUE)
  partitions_copy <- partitions
  for (i in 1:length(partitions_copy)) partitions[partitions_copy==new_partitions_order[i]] <- i
  
  if(t=="subs"){
    closestCosmic <- findClosestCOSMIC30andCombinations_withSimilarity(all_signatures[["subs"]])
    closestCosmic_label <- paste0("Closest Cosmic: ",closestCosmic$cosmic," (",sprintf("%.2f",closestCosmic$cos.sim),")")
    names(closestCosmic_label) <- colnames(all_signatures[["subs"]])
  }else if (t=="rearr"){
    closestCosmic <- findClosestRearrSigsBreast560_withSimilarity(all_signatures[["rearr"]])
    closestCosmic_label <- paste0("Closest Breast560: ",closestCosmic$RS.Breast560," (",sprintf("%.2f",closestCosmic$cos.sim),")")
    names(closestCosmic_label) <- colnames(all_signatures[["rearr"]])
  }

  # partitions <- cutree(hc_fit, h = where_to_cut)
  partition_list <- unique(partitions)
  boxes_order <- partitions[hc_fit$order]
  #plot
  output_file <- paste0(out_dir,"_",t,"_Cluster_Average_",tag,".jpg")
  jpeg(output_file,width = 6000,height = 1000,res = 160)
  plot(hc_fit,hang = -1,xlab = "",ylab = "1 - cosine similarity, average linkage",sub = "")
  if (!is.null(where_to_cut)) abline(a=where_to_cut,b=0,col="red")
  start_draw <- 0.5
  if(t=="subs"){
    box_bottom <- -0.48
  }else if (t=="rearr"){
    box_bottom <- -0.55
  }
  for (i in unique(boxes_order)){
    end_draw <- start_draw + sum(boxes_order==i)
    par(xpd=TRUE)
    rect(start_draw,box_bottom,end_draw,0,border = "red")
    text(i,x = start_draw + (end_draw - start_draw)/2,y = box_bottom - 0.05,cex = 0.8)
    par(xpd=FALSE)
    start_draw <- end_draw
  }
  dev.off()
  
  #plot all partitions together
  out_dir_par <- paste0(out_dir,"_",t,"_partitions_",tag,"/")
  dir.create(out_dir_par,showWarnings = FALSE,recursive = TRUE)
  for (i in partition_list){
    if(t=="subs"){
      plotSubsSignatures(all_signatures[[t]][,names(partitions[partitions==i]),drop = FALSE],
                         output_file = paste0(out_dir_par,"partition_",i,".jpg"),
                         plot_sum = FALSE,overall_title = paste0("Pan Can Cluster ",i),
                         #add_to_titles = closestCosmic_label[names(partitions[partitions==i])])
                         add_to_titles = NULL)
    }else if (t=="rearr"){
      plotRearrSignatures(all_signatures[[t]][,names(partitions[partitions==i]),drop = FALSE],
                          output_file = paste0(out_dir_par,"partition_",i,".jpg"),
                          plot_sum = FALSE,overall_title = paste0("Pan Can Cluster ",i),
                          #add_to_titles = closestCosmic_label[names(partitions[partitions==i])])
                          add_to_titles = NULL)
    }
    writeTable(all_signatures[[t]][,names(partitions[partitions==i]),drop = FALSE],
               file = paste0(out_dir_par,"partition_",i,".tsv"))
  }
  
  #average signatures (Reference Signatures)
  mean_sig_matrix <- matrix(NA,nrow = nrow(all_signatures[[t]]),ncol=max(partitions),
                            dimnames = list(rownames(all_signatures[[t]]),1:max(partitions)))
  sd_sig_matrix <- matrix(NA,nrow = nrow(all_signatures[[t]]),ncol=max(partitions),
                          dimnames = list(rownames(all_signatures[[t]]),1:max(partitions)))
  partition_size <- c()
  for (i in 1:max(partitions)){
    current_sigs <- names(partitions)[partitions==i]
    partition_size <- c(partition_size,length(current_sigs))
    mean_sig_matrix[,i] <- apply(all_signatures[[t]][,current_sigs,drop=FALSE],1,mean)
    if(length(current_sigs)>1){
      sd_sig_matrix[,i] <- apply(all_signatures[[t]][,current_sigs,drop=FALSE],1,sd)
    }else{
      sd_sig_matrix[,i] <- 0
    }
  }
  if(t=="subs"){
    #closestCosmic <- findClosestCOSMIC30_withSimilarity(mean_sig_matrix)
    closestCosmic <- findClosestCOSMIC_PCAWG_withSimilarity(mean_sig_matrix)
    #closestCosmic_label <- paste0("Closest COSMIC: ",closestCosmic$cosmic," (",sprintf("%.2f",closestCosmic$cos.sim),")")
    closestCosmic_label <- paste0(" (n=",partition_size,"): RefSig",subs_sigs_names,"\nClosest COSMIC: ",closestCosmic$cosmic," (",sprintf("%.2f",closestCosmic$cos.sim),")")
    names(closestCosmic_label) <- colnames(mean_sig_matrix)
    
    # out_dir_par <- paste0(out_dir,"_subs_partitions_cuth",where_to_cut,"/")
    plotSubsSignatures_MeanSd(mean_matrix = mean_sig_matrix,sd_matrix = sd_sig_matrix,
                              output_file = paste0(out_dir_par,"partitions_mean_sd.jpg"),
                              add_to_titles = paste0(closestCosmic_label),
                              plot_sum = FALSE)
  }else if (t=="rearr"){
    closestCosmic <- findClosestRearrSigsBreast560_withSimilarity(mean_sig_matrix)
    #closestCosmic_label <- paste0("Closest Breast560: ",closestCosmic$RS.Breast560," (",sprintf("%.2f",closestCosmic$cos.sim),")")
    closestCosmic_label <- paste0(" (n=",partition_size,"): RefSig",rearr_sigs_names,"\nClosest Breast560: ",closestCosmic$RS.Breast560," (",sprintf("%.2f",closestCosmic$cos.sim),")")
    names(closestCosmic_label) <- colnames(mean_sig_matrix)
    
    # out_dir_par <- paste0(out_dir,"_rearr_partitions_cuth",where_to_cut,"/")
    plotRearrSignatures_MeanSd(mean_matrix = mean_sig_matrix,sd_matrix = sd_sig_matrix,
                              output_file = paste0(out_dir_par,"partitions_mean_sd.jpg"),
                              add_to_titles = paste0(closestCosmic_label),
                              plot_sum = FALSE)
  }
  mean_cosmic_similar <- data.frame(closestCosmic,row.names = paste0("partition.",1:length(closestCosmic$cos.sim)))
  writeTable(mean_sig_matrix,file = paste0(out_dir_par,"partitions_mean.tsv"))
  writeTable(sd_sig_matrix,file = paste0(out_dir_par,"partitions_sd.tsv"))
  writeTable(partition_size,file = paste0(out_dir_par,"partitions_n.tsv"))
  writeTable(mean_cosmic_similar,file = paste0(out_dir_par,"partitions_CosSim.tsv"))
}


plotSubsSignatures_MeanSd <- function(mean_matrix,sd_matrix,output_file = NULL,plot_sum = TRUE,overall_title = "",add_to_titles = NULL,mar=NULL){
  colnames(mean_matrix) <- sapply(colnames(mean_matrix),function(x) if (nchar(x)>20) paste0(substr(x,1,17),"...") else x)
  rearr.colours <- c(rep("blue",16),rep("black",16),rep("red",16),rep("grey",16),rep("green",16),rep("pink",16))
  nplotrows <- ceiling(ncol(mean_matrix)/3)
  if(!is.null(output_file)) jpeg(output_file,width = 3*800,height = nplotrows*400,res = 220)
  par(mfrow = c(nplotrows, 3),oma=c(0,0,2,0))
  if(is.null(mar)){
    par(mar=c(5,3,2,2))
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
  for (pos in 1:ncol(mean_matrix)){
    ylimit <- c(0,max(mean_matrix[,pos]+sd_matrix[,pos])+0.1*max(mean_matrix[,pos]+sd_matrix[,pos]))
    title <- paste0("Group ",pos)
    if (plot_sum) title <- paste0(title," (",round(sum(mean_matrix[,pos]))," substitutions)")
    if (!is.null(add_to_titles)) title <- paste0(title,"",add_to_titles[pos])
    barCenters <- barplot(mean_matrix[,pos],
                          main = title,
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



plotRearrSignatures_MeanSd <-function(mean_matrix,sd_matrix,output_file = NULL,plot_sum = TRUE,overall_title = "",add_to_titles = NULL,mar=NULL){
  colnames(mean_matrix) <- sapply(colnames(mean_matrix),function(x) if (nchar(x)>20) paste0(substr(x,1,17),"...") else x)
  del_col = rgb(228,26,28, maxColorValue = 255)
  td_col = rgb(77,175,74, maxColorValue =255)
  inv_col  = rgb(55,126,184, maxColorValue = 255)
  transloc_col = rgb(152,78,163, maxColorValue =255)
  non_clust_col = rgb(240,240,240, maxColorValue =255)
  #rearr.colours <- c(rep("darkblue",16),rep("red",16))
  rearr.colours <- rep(c(rep(del_col,5),rep(td_col,5),rep(inv_col,5),transloc_col),2)
  nplotrows <- ceiling(ncol(mean_matrix)/3)
  if(!is.null(output_file)) jpeg(output_file,width = 3*800,height = nplotrows*500,res = 220)
  par(mfrow = c(nplotrows, 3),oma=c(0,0,2,0))
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
  
  for (pos in 1:ncol(mean_matrix)){
    ylimit <- c(0,max(mean_matrix[,pos]+sd_matrix[,pos])+0.1*max(mean_matrix[,pos]+sd_matrix[,pos]))
    if(is.null(mar)){
      par(mar=c(8,3,2,2))
    }else{
      par(mar=mar)
    }
    title <- paste0("Group ",pos)
    if (plot_sum) title <- paste0(title," (",sum(mean_matrix[,pos])," rearrangements)")
    if (!is.null(add_to_titles)) title <- paste0(title,"",add_to_titles[pos])
    barCenters <- barplot(mean_matrix[,pos],
                          main = title,
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

#functions end----







#Which results?
choice_table <- list()
choice_table[["subs"]] <- read.table("../data/2018_02_18_fullExtraction_pancan_signatures_subs.txt",
                                header = TRUE,as.is = TRUE,check.names = FALSE,sep = "\t")
choice_table[["rearr"]] <- read.table("../data/2018_02_18_fullExtraction_pancan_signatures_rearr.txt",
                                     header = TRUE,as.is = TRUE,check.names = FALSE,sep = "\t")

#load COSMIC signatures
signatures_subs_COSMIC <- sortCatalogue(read.table(paste0("../data/COSMIC_signatures.txt"),
                                sep = "\t",header = TRUE,as.is = TRUE,check.names = FALSE))
colnames(signatures_subs_COSMIC) <- paste0(rep("COSMIC-",ncol(signatures_subs_COSMIC)),seq(1,ncol(signatures_subs_COSMIC)))

signatures_rearr_breast560 <- read.table(paste0("../data/rearrangement.signatures.txt"),
                                       sep = "\t",header = TRUE,as.is = TRUE,check.names = FALSE)
colnames(signatures_rearr_breast560) <- paste0(rep("Breast560-",ncol(signatures_rearr_breast560)),seq(1,ncol(signatures_rearr_breast560)))

all_signatures <- list()
all_signatures[["subs"]] <- data.frame(row.names = row.names(signatures_subs_COSMIC))
#all_signatures[["subs"]] <- cbind(all_signatures[["subs"]],signatures_subs_COSMIC)
all_signatures[["rearr"]] <- data.frame(row.names = row.names(signatures_rearr_breast560))
# all_signatures[["rearr"]] <- cbind(all_signatures[["rearr"]],signatures_rearr_breast560)
for(t in c("subs","rearr")){
  for (group in choice_table[[t]]$Organ){

    ns <- choice_table[[t]]$nsig[choice_table[[t]]$Organ==group]
    tag <- ""
    if (choice_table[[t]]$scale[choice_table[[t]]$Organ==group]) tag <- paste0(tag,"_scale")
    signatures <- read.table(paste0("../results/",t,"/extraction/brunet/extraction_",group,"_",t,"_brunet",tag,"/round_1/sig_",ns,"/Sigs_plot_",group,"_ns",ns,"_nboots",nboots,".tsv"),
                                    sep = "\t",header = TRUE,as.is = TRUE,check.names = FALSE)
    if (t=="subs") signatures <- sortCatalogue(signatures)
    if (group=="Breast741"){
      colnames(signatures) <- paste0(rep("Breast",ncol(signatures)),"_",seq(1,ncol(signatures)))
    }else{
      colnames(signatures) <- paste0(rep(group,ncol(signatures)),"_",seq(1,ncol(signatures)))
    }
    
    all_signatures[[t]] <- cbind(all_signatures[[t]],signatures)
    
  }
}

colnames(all_signatures[["subs"]]) <- change_numbers_to_letters(colnames(all_signatures[["subs"]]))
colnames(all_signatures[["rearr"]]) <- change_numbers_to_letters(colnames(all_signatures[["rearr"]]))

#Hierachical clustering subs
#compute the distance matrix based on the cosine similarity matrix
# distMatrix_subs <- 1 - computeCorrelation(all_subs_signatures)
# #fit
# fit_subs <- hclust(as.dist(distMatrix_subs), method="ward.D") 
# #plot
# output_file <- paste0(results_dir,dateOfAnalysis,"_Subs_Cluster_Ward.jpg")
# jpeg(output_file,width = 3000,height = 500,res = 80)
# plot(fit_subs)
# dev.off()
#fit
# closestCosmic <- findClosestCOSMIC30andCombinations_withSimilarity(all_signatures[["subs"]])
# closestCosmic_label <- paste0("Closest Cosmic: ",closestCosmic$cosmic," (",sprintf("%.2f",closestCosmic$cos.sim),")")
# names(closestCosmic_label) <- colnames(all_signatures[["subs"]])
# 
# closestBreast560 <- findClosestRearrSigsBreast560_withSimilarity(all_signatures[["rearr"]])
# closestBreast560_label <- paste0("Closest Breast560: ",closestBreast560$RS.Breast560," (",sprintf("%.2f",closestBreast560$cos.sim),")")
# names(closestBreast560_label) <- colnames(all_signatures[["rearr"]])


distMatrix_subs <- 1 - computeCorrelation(all_signatures[["subs"]])
fit_subs <- hclust(as.dist(distMatrix_subs), method="average") 

distMatrix_rearr <- 1 - computeCorrelation(all_signatures[["rearr"]])
fit_rearr <- hclust(as.dist(distMatrix_rearr), method="average") 



where_to_cut <- 0.25
#plot dendrogram and partitions
# plot_clustering_of_pancan(out_dir = out_dir,partitions = cutree(fit_subs, h = where_to_cut),
#                           hc_fit = fit_subs,t = "subs",where_to_cut = where_to_cut,tag = paste0("cuth",where_to_cut))
# plot_clustering_of_pancan(out_dir = out_dir,partitions = cutree(fit_rearr, h = where_to_cut),
#                           hc_fit = fit_rearr,t = "rearr",where_to_cut = where_to_cut,tag = paste0("cuth",where_to_cut))


#now some manual curation of the subs clustering following dendrogram
#I just need to change the partitions object
partitions_subs <- cutree(fit_subs, h = where_to_cut)
sigs_names <- names(partitions_subs)
partitions_subs <- partitions_subs[fit_subs$order]
partitions_subs_copy <- partitions_subs

partitions_subs[partitions_subs_copy==26] <- 1
partitions_subs[partitions_subs_copy==31] <- 2
partitions_subs[partitions_subs_copy==30] <- 3
partitions_subs[partitions_subs_copy==1] <- 4
partitions_subs[partitions_subs_copy==35] <- 5
partitions_subs[partitions_subs_copy==21][1:2] <- 6
partitions_subs[partitions_subs_copy==21][3:4] <- 7
partitions_subs[partitions_subs_copy==10][1:8] <- 8
partitions_subs[partitions_subs_copy==10][9:14] <- 9
partitions_subs[partitions_subs_copy==15] <- 10
partitions_subs[partitions_subs_copy==37] <- 11
partitions_subs[partitions_subs_copy==3][1:3] <- 12
partitions_subs[partitions_subs_copy==3][4:14] <- 13
partitions_subs[partitions_subs_copy==22] <- 14
partitions_subs[partitions_subs_copy==38] <- 15
partitions_subs[partitions_subs_copy==34] <- 16
partitions_subs[partitions_subs_copy==2][1] <- 17
partitions_subs[partitions_subs_copy==2][2:23] <- 18
partitions_subs[partitions_subs_copy==39] <- 19
partitions_subs[partitions_subs_copy==16] <- 20
partitions_subs[partitions_subs_copy==32] <- 21
partitions_subs[partitions_subs_copy==25] <- 22
partitions_subs[partitions_subs_copy==5] <- 23
partitions_subs[partitions_subs_copy==33] <- 24
partitions_subs[partitions_subs_copy==13] <- 24
partitions_subs[partitions_subs_copy==12] <- 25
partitions_subs[partitions_subs_copy==11] <- 25
partitions_subs[partitions_subs_copy==29] <- 26
partitions_subs[partitions_subs_copy==24] <- 27
partitions_subs[partitions_subs_copy==36][1] <- 28
partitions_subs[partitions_subs_copy==36][2] <- 29
partitions_subs[partitions_subs_copy==9][1] <- 30
partitions_subs[partitions_subs_copy==9][2:3] <- 31
partitions_subs[partitions_subs_copy==4] <- 31
partitions_subs[partitions_subs_copy==23] <- 32
partitions_subs[partitions_subs_copy==17] <- 33
partitions_subs[partitions_subs_copy==8] <- 34
partitions_subs[partitions_subs_copy==6] <- 35
partitions_subs[partitions_subs_copy==28] <- 36
partitions_subs[partitions_subs_copy==7] <- 37
partitions_subs[partitions_subs_copy==19] <- 38
partitions_subs[partitions_subs_copy==18] <- 39
partitions_subs[partitions_subs_copy==14] <- 39
partitions_subs[partitions_subs_copy==27] <- 40
partitions_subs[partitions_subs_copy==20] <- 41

recover_order <- array(NA,length(partitions_subs))
for (i in 1:length(partitions_subs)){
  recover_order[fit_subs$order[i]] <- partitions_subs[i]
}
names(recover_order) <- sigs_names
partitions_subs <- recover_order

#plot dendrogram and partitions
plot_clustering_of_pancan(out_dir = out_dir,partitions = partitions_subs,
                          hc_fit = fit_subs,t = "subs",where_to_cut = NULL,tag = "manually_updated")


#-----------------


#now some manual curation of the rearr clustering following dendrogram
#I just need to change the partitions object
partitions_rearr <- cutree(fit_rearr, h = where_to_cut)
sigs_names <- names(partitions_rearr)
partitions_rearr <- partitions_rearr[fit_rearr$order]
partitions_rearr_copy <- partitions_rearr

partitions_rearr[partitions_rearr_copy==4][1] <- 1
partitions_rearr[partitions_rearr_copy==4][2:7] <- 2
partitions_rearr[partitions_rearr_copy==9] <- 3
partitions_rearr[partitions_rearr_copy==3][1:5] <- 4
partitions_rearr[partitions_rearr_copy==3][6:14] <- 5
partitions_rearr[partitions_rearr_copy==17] <- 6
partitions_rearr[partitions_rearr_copy==11] <- 7
partitions_rearr[partitions_rearr_copy==12] <- 8
partitions_rearr[partitions_rearr_copy==2] <- 9
partitions_rearr[partitions_rearr_copy==18] <- 10
partitions_rearr[partitions_rearr_copy==10] <- 11
partitions_rearr[partitions_rearr_copy==8] <- 12
partitions_rearr[partitions_rearr_copy==15] <- 13
partitions_rearr[partitions_rearr_copy==7] <- 14
partitions_rearr[partitions_rearr_copy==16] <- 15
partitions_rearr[partitions_rearr_copy==19] <- 16
partitions_rearr[partitions_rearr_copy==1] <- 17
partitions_rearr[partitions_rearr_copy==6] <- 18
partitions_rearr[partitions_rearr_copy==13] <- 19
partitions_rearr[partitions_rearr_copy==14] <- 20
partitions_rearr[partitions_rearr_copy==5] <- 21

recover_order <- array(NA,length(partitions_rearr))
for (i in 1:length(partitions_rearr)){
  recover_order[fit_rearr$order[i]] <- partitions_rearr[i]
}
names(recover_order) <- sigs_names
partitions_rearr <- recover_order

#plot dendrogram and partitions
plot_clustering_of_pancan(out_dir = out_dir,partitions = partitions_rearr,
                          hc_fit = fit_rearr,t = "rearr",where_to_cut = NULL,tag = "manually_updated")


#-----------------
