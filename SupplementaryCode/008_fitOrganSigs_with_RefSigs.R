#requires signature.tools.lib package
setwd("~/sandbox/git/signature-tools/bin/")
source("../lib/sigantureAssignmentLib.R")

type_mut <- "subs"
#type_mut <- "rearr"

replot <- FALSE

max_combinations <- 3

if(type_mut=="subs"){
  pancan_sigs <- read.table("../results/2018_05_21_pancan_signatures_subs.tsv",sep = "\t",stringsAsFactors = FALSE,header = TRUE,check.names = FALSE)
  mean_sigs <- read.table("../results/pancan_extraction_clustering/2018_09_26_subs_partitions_manually_updated/partitions_mean.tsv",sep = "\t",stringsAsFactors = FALSE,header = TRUE,check.names = FALSE)
  mean_sigs <- mean_sigs[,- c(12,16,22)] #remove redundant
}else if(type_mut=="rearr"){
  pancan_sigs <- read.table("../results/2018_05_21_pancan_signatures_rearr.tsv",sep = "\t",stringsAsFactors = FALSE,header = TRUE,check.names = FALSE)
  mean_sigs <- read.table("../results/pancan_extraction_clustering/2018_09_26_rearr_partitions_manually_updated/partitions_mean.tsv",sep = "\t",stringsAsFactors = FALSE,header = TRUE,check.names = FALSE)
}
out_dir <- "../results/FitClusterAveragesToPanCanSigs/"
dir.create(paste0(out_dir,"/",type_mut),showWarnings = FALSE,recursive = TRUE)

signatures_to_fit <- pancan_sigs
signatures_to_use <- mean_sigs

full_results_file <- paste0(out_dir,type_mut,"/2018_10_12_full_results_Fits.rData")
if(!file.exists(full_results_file)){
  final_res <- fitSignaturesWithSignatures(signatures_to_fit,
                                           signatures_to_use,
                                           max_combinations = max_combinations,
                                           nparallel = 8)
  save(file = full_results_file,final_res,signatures_to_fit,signatures_to_use,max_combinations)
}else{
  load(file = full_results_file)
}

#plotting!
sigs_names <- names(final_res$fit_tables)
if(replot){
  for(current_s in sigs_names){
    
    current_fit_table <- final_res$fit_tables[[current_s]]
    current_cossim <- final_res$cossim_list[[current_s]]
    current_comb_class <- final_res$combinations_class_list[[current_s]]
    
    jpeg(filename = paste0(out_dir,type_mut,"/Average_cluster_fit_",current_s,".jpg"),
         width = 640*(nrow(current_fit_table)),
         height = 480*5,
         res = 150)
    
    par(mfrow=c(5,nrow(current_fit_table)))
    for(i in 1:nrow(current_fit_table)){
      if (type_mut=="subs"){
        plotSubsSignatures(signatures_to_fit[,current_s,drop=FALSE],plot_sum = FALSE,add_to_titles = "original")
      }else if (type_mut=="rearr"){
        plotRearrSignatures(signatures_to_fit[,current_s,drop=FALSE],plot_sum = FALSE,add_to_titles = "original")
      }
    }
    if (type_mut=="subs"){
      plotSubsSignatures(as.matrix(signatures_to_use) %*% t(current_fit_table),plot_sum = FALSE,add_to_titles = paste0("model (cos sim=",sprintf("%.2f",current_cossim),")"))
    }else if (type_mut=="rearr"){
      plotRearrSignatures(as.matrix(signatures_to_use) %*% t(current_fit_table),plot_sum = FALSE,add_to_titles = paste0("model (cos sim=",sprintf("%.2f",current_cossim),")"))
    }
    track_plots <- rep(0,nrow(current_fit_table))
    for(c in 1:max_combinations){
      for(i in 1:nrow(current_fit_table)){
        track_plots[i] <- track_plots[i]+1
        pos <- which(current_fit_table[i,]>0)[track_plots[i]]
        if(!is.na(pos)){
          if (type_mut=="subs"){
            plotSubsSignatures(signatures_to_use[,pos,drop=FALSE],plot_sum = FALSE)
          }else if (type_mut=="rearr"){
            plotRearrSignatures(signatures_to_use[,pos,drop=FALSE],plot_sum = FALSE)
          }
        }else{
          #fake plot to fill the space
          par(mar=c(0,0,0,0))
          plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
        }
      }
    }
    dev.off()
    
  }
}
  
#------------------------------------------------------------------------------ 
#------ now let's build the conversion matrices based on the manual selection
#------------------------------------------------------------------------------ 

manual_selection <- read.table(paste0("../results/FitClusterAveragesToPanCanSigs/",type_mut,"/2018_10_13_SignatureAssignment_",type_mut,".txt"),sep = '\t',
                               header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
manual_selection_list <- list()

conversion_table <- NULL
for (i in 1:nrow(manual_selection)){
  sig <- manual_selection$signature[i]
  group_names <- unlist(strsplit(manual_selection$groups[i],","))
  
  if(length(group_names)>1){
    current_fit_table <- final_res$fit_tables[[sig]]
    current_cossim <- final_res$cossim_list[[sig]]
    current_comb_class <- final_res$combinations_class_list[[sig]]
    
    rowselection <- apply(current_fit_table[,group_names,drop=FALSE] > 0,1,all) & current_comb_class==length(group_names)
    newrow <- current_fit_table[rowselection,,drop=FALSE]
    rownames(newrow) <- sig
  }else{
    newrow <- data.frame(matrix(0,ncol = ncol(signatures_to_use),nrow = 1,dimnames = list(sig,colnames(signatures_to_use))),check.names = FALSE)
    newrow[1,group_names] <- 1
  }
  conversion_table <- rbind(conversion_table,newrow)
}

write.table(conversion_table,file = paste0(out_dir,type_mut,"/Conversion_table_",type_mut,".tsv"),
            sep="\t",col.names = TRUE,row.names = TRUE,quote = FALSE)


#------------------------------------------------------------------------------ 
#------ now plot conversion matrices per organ
#------------------------------------------------------------------------------

draw_legend <- function(col,xl,xr,yb,yt){
  par(xpd=TRUE)
  rect(xl,yb,xr,yt)
  rect(
    xl,
    head(seq(yb,yt,(yt-yb)/length(col)),-1),
    xr,
    tail(seq(yb,yt,(yt-yb)/length(col)),-1),
    col=col,border = NA
  )
  text(x = xl-0.05*xl, y = yt,labels = "1")
  text(x = xl-0.05*xl, y = (yt-yb)/2,labels = "0.5")
  text(x = xl-0.05*xl, y = yb,labels = "0")
  par(xpd=FALSE)
}

out_dir_plot_matrices <- paste0(out_dir,type_mut,"/conversion_matrices/")
dir.create(out_dir_plot_matrices,showWarnings = FALSE,recursive = TRUE)

conversion_table <- read.table(file = paste0(out_dir,type_mut,"/Conversion_table_",type_mut,".tsv"),
                                sep="\t",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)

choice_table <- read.table(paste0("../data/2018_02_18_fullExtraction_pancan_signatures_",type_mut,".txt"),
                           header = TRUE,as.is = TRUE,check.names = FALSE,sep = "\t")
organs <- choice_table$Organ

for (i in 1:nrow(choice_table)){
  organ <- choice_table$Organ[i]
  nsigs <- choice_table$nsig[i]
  sigs_names <- paste0(organ,"_",1:nsigs)
  tmp_matrix <- conversion_table[sigs_names,,drop=FALSE]
  tmp_matrix <- tmp_matrix[,apply(tmp_matrix,2,sum)>0,drop=FALSE]
  colnames(tmp_matrix) <- paste0("Ref_",type_mut,"_",colnames(tmp_matrix))
  
  #plot it!
  jpeg(filename = paste0(out_dir_plot_matrices,"conversionMatrix_",organ,".jpg"),
       width = 800,
       height = 600,
       res = 140)
  current_table <- as.matrix(t(tmp_matrix))
  par(mar=c(8,9,3,6))
  par(xpd=FALSE)
  col<- colorRampPalette(c("white", "blue"))(51)
  image(current_table,col = col,zlim = c(0,1), axes=F,main="Signatures conversion matrix",cex.main = 0.9)
  extrabit_v <- 1/(ncol(current_table)-1)/2
  if (ncol(current_table)>1) abline(h=seq(0-extrabit_v,1+extrabit_v,length.out = ncol(current_table)+1),col="grey",lty=2)
  extrabit_h <- 1/(nrow(current_table)-1)/2
  if (nrow(current_table)>1) abline(v=seq(0-extrabit_h,1+extrabit_h,length.out = nrow(current_table)+1),col="grey",lty=2)
  axis(2,at = seq(0,1,length.out = ncol(current_table)),labels = colnames(current_table),las=1,cex.lab=0.8,tick = FALSE)
  axis(1,at = seq(0,1,length.out = nrow(current_table)),labels = rownames(current_table),las=2,cex.lab=0.8,tick = FALSE)
  par(xpd=TRUE)
  #text(labels = "combinations",x = 0.5,y = - 0.1)
  draw_legend(col,1.25,1.3,0,1)
  dev.off()
}

