library(signature.tools.lib)

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

thresholds_perc <- 5
nboots_fit <- 100
p.value <- 0.05

#collect all the res.cor_triangular matrices
res.cor_triangular_list <- list()

for (k in 1:n){
  #load data
  cat_file <- paste0(resultdir,"simulated_catalogue_",k,".tsv")
  E_file <- paste0(resultdir,"simulated_E_",k,".tsv")
  if (file.exists(cat_file)){
    catalogue <- read.table(file = cat_file,sep = "\t",header = TRUE,check.names = FALSE,as.is = TRUE)
    E <- read.table(file = E_file,sep = "\t",header = TRUE,check.names = FALSE,as.is = TRUE)
    message("catalogue ",k," loaded from file")
  } else {
    message("catalogue",k,"not found")
  }
  
  #Fit analysis
  outdir <- paste0(resultdir,"fitAnalysis_part2/rep_",k,"/")
  dir.create(outdir,showWarnings = FALSE,recursive = TRUE)
  # fit_res_file <- paste0(resultdir,"with_filterBest_brunet_",k,"/round_1/sig_10/",
  #        "SigFit_withBootstrap_Summary_mKLD_bfmCosSim_alpha-1_tr",
  #        thresholds_perc,"_p",p.value,"_MC.rData")
  # file.copy(fit_res_file,paste0(outdir,"SigFit_withBootstrap_Summary_mKLD_bfmCosSim_alpha-1_tr",
  #                               thresholds_perc,"_p",p.value,".rData"))
  res <- SignatureFit_withBootstrap_Analysis(outdir = outdir,
                                             cat = catalogue,
                                             signature_data_matrix = P,
                                             nboot = nboots_fit,
                                             type_of_mutations = "subs",
                                             threshold_percent = thresholds_perc,
                                             threshold_p.value = p.value,
                                             nparallel = nparallel)
  unexpl_tab <- unexplainedSamples(fileout = paste0(outdir,"unexplained.jpg"),
                     catalogue = catalogue,
                     sigs = P,
                     exposures = res$E_median_filtered,
                     pvalue_threshold = 0.01)
  
  
  #plots copied from SignatureFitLib and adapted 
  cat <- catalogue
  signature_data_matrix <- P
  nboot <- nboots_fit
  type_of_mutations <- "subs"
  threshold_percent <- thresholds_perc
  threshold_p.value <- p.value
  
  #function to draw a legend for the heatmap of the correlation matrix
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
    text(x = 1.2, y = yt,labels = "1")
    text(x = 1.2, y = (yt-yb)/2,labels = "0")
    text(x = 1.2, y = yb,labels = "-1")
  }
  
  reconstructed_with_median <- round(as.matrix(signature_data_matrix) %*% res$E_median_filtered)
  #provide a series of plots for each sample
  #plot_nrows <- ncol(cat)
  rows_ordered_from_best <- order(res$KLD_samples)
  plot_nrows <- 2
  plot_ncol <- 4
  nplots <- plot_nrows*plot_ncol
  which_to_plot <- unexpl_tab$index
  howmanyplots <- ncol(cat)
  plotsdir <- paste0(outdir)
  dir.create(plotsdir,recursive = TRUE,showWarnings = FALSE)
  
  #collect the res.cor_triangular
  for (i in 1:howmanyplots){
    if(sum(cat[,i,drop=FALSE])>0){
      if(ncol(signature_data_matrix)>1){
        #5 top correlated signatures
        res.cor <- cor(t(res$samples_list[[i]]),method = "spearman")
        res.cor_triangular <- res.cor
        res.cor_triangular[row(res.cor)+(ncol(res.cor)-col(res.cor))>=ncol(res.cor)] <- 0
        res.cor_triangular_list[[i+howmanyplots*(k-1)]] <- res.cor_triangular
      }
    }
  }
  
  for(p in which_to_plot){
    # if (plot_nrows+(p-1)*plot_nrows>=length(rows_ordered_from_best)){
    #   current_samples <- rows_ordered_from_best[1:plot_nrows+(p-1)*plot_nrows]
    # }else{
    #   current_samples <- rows_ordered_from_best[(1+(p-1)*plot_nrows):length(rows_ordered_from_best)]
    # }
    current_samples <- p
    jpeg(filename = paste0(plotsdir,"sigfit_bootstrap_",p,"of",howmanyplots,".jpg"),
         width = 800*(plot_ncol),
         height = 600*plot_nrows,
         res = 220)
    par(mfrow=c(plot_nrows,plot_ncol),oma=c(0,1,0,0))
    for(i in current_samples){
      unassigned_mut <- sprintf("%.2f",(sum(cat[,i,drop=FALSE]) - sum(reconstructed_with_median[,i,drop=FALSE]))/sum(cat[,i,drop=FALSE])*100)
      percentdiff <- sprintf("%.2f",sum(abs(cat[,i,drop=FALSE] - reconstructed_with_median[,i,drop=FALSE]))/sum(cat[,i,drop=FALSE])*100) 
      if(type_of_mutations=="subs"){
        #1 original
        plotSubsSignatures(signature_data_matrix = cat[,i,drop=FALSE],add_to_titles = "Catalogue",mar=c(6,3,5,2))
        if(sum(cat[,i,drop=FALSE])>0){
          #2 reconstructed
          plotSubsSignatures(signature_data_matrix = reconstructed_with_median[,i,drop=FALSE],add_to_titles = "Model",mar=c(6,3,5,2))
          #3 difference
          plotSubsSignatures(signature_data_matrix = cat[,i,drop=FALSE] - reconstructed_with_median[,i,drop=FALSE],add_to_titles = paste0("Unassigned, ",unassigned_mut,"%"),mar=c(6,3,5,2),plot_sum = FALSE)
        }
      }else if(type_of_mutations=="rearr"){
        #1 original
        plotRearrSignatures(signature_data_matrix = cat[,i,drop=FALSE],add_to_titles = "Catalogue",mar=c(12,3,5,2))
        if(sum(cat[,i,drop=FALSE])>0){
          #2 reconstructed
          plotRearrSignatures(signature_data_matrix = reconstructed_with_median[,i,drop=FALSE],add_to_titles = "Model",mar=c(12,3,5,2))
          #3 difference
          plotRearrSignatures(signature_data_matrix = cat[,i,drop=FALSE] - reconstructed_with_median[,i,drop=FALSE],add_to_titles = paste0("Unassigned, ",unassigned_mut,"%"),mar=c(12,3,5,2),plot_sum = FALSE)
        }
      }else if(type_of_mutations=="generic"){
        #1 original
        plotGenericSignatures(signature_data_matrix = cat[,i,drop=FALSE],add_to_titles = "Catalogue",mar=c(6,3,5,2))
        if(sum(cat[,i,drop=FALSE])>0){
          #2 reconstructed
          plotGenericSignatures(signature_data_matrix = reconstructed_with_median[,i,drop=FALSE],add_to_titles = "Model",mar=c(6,3,5,2))
          #3 difference
          plotGenericSignatures(signature_data_matrix = cat[,i,drop=FALSE] - reconstructed_with_median[,i,drop=FALSE],add_to_titles = paste0("Unassigned, ",unassigned_mut,"%"),mar=c(6,3,5,2),plot_sum = FALSE)
        }
      }
      if(sum(cat[,i,drop=FALSE])>0){
        #4 bootstraps
        par(mar=c(6,4,5,2))
        boxplot(t(res$samples_list[[i]]),las=3,cex.axes=0.9,
                ylab="n mutations",
                ylim=c(0,max(res$samples_list[[i]])),
                main=paste0("Exposures, of ",colnames(res$E_median_filtered)[i],"\nthreshold=",threshold_percent,"%, p-value=",threshold_p.value,", n=",nboot))
        points(1:length(res$E_median_filtered[,i,drop=FALSE]),res$E_median_filtered[,i,drop=FALSE],col="red")
        points(1:length(E[,i,drop=TRUE]),E[,i,drop=TRUE],col="#0000FF", cex = 1.8,pch = 4)
        abline(h=threshold_percent/100*sum(cat[,i,drop=FALSE]),col="green")
        ylegend <- max(res$samples_list[[i]])/720*875
        legend(x=-0.5,y=ylegend,legend = c("consensus exposures"),col = "red",pch = 1,cex = 1,bty = "n",inset = c(0,-0.16),xpd = TRUE)
        legend(x=8,y=ylegend,legend = c("threshold"),col = "green",lty = 1,cex = 1,bty = "n",inset = c(0,-0.16),xpd = TRUE)
        legend(x=5.3,y=ylegend,legend = c("original"),col = "#0000FF",pch = 4,cex = 1,bty = "n",xpd = TRUE)
        if(ncol(signature_data_matrix)>1){
          #5 top correlated signatures
          res.cor <- cor(t(res$samples_list[[i]]),method = "spearman")
          res.cor_triangular <- res.cor
          res.cor_triangular[row(res.cor)+(ncol(res.cor)-col(res.cor))>=ncol(res.cor)] <- 0
          res.cor_triangular_label <- matrix(sprintf("%0.2f",res.cor_triangular),nrow = nrow(res.cor_triangular))
          res.cor_triangular_label[row(res.cor)+(ncol(res.cor)-col(res.cor))>=ncol(res.cor)] <- ""
          # heatmap(res.cor_triangular,
          #           Rowv = NA,
          #           Colv = NA,
          #           scale = "none",
          #           col = col,
          #           symm = TRUE,
          #           breaks=seq(-1,1,length.out = 52))
          par(mar=c(6,8,5,6))
          par(xpd=FALSE)
          col<- colorRampPalette(c("blue", "white", "red"))(51)
          image(res.cor_triangular,col = col,zlim = c(-1,1), axes=F,main="Exposures Correlation (spearman)")
          extrabit <- 1/(ncol(signature_data_matrix)-1)/2
          abline(h=seq(0-extrabit,1+extrabit,length.out = ncol(signature_data_matrix)+1),col="grey",lty=2)
          abline(v=seq(0-extrabit,1+extrabit,length.out = ncol(signature_data_matrix)+1),col="grey",lty=2)
          axis(2,at = seq(0,1,length.out = ncol(signature_data_matrix)),labels = colnames(signature_data_matrix),las=1,cex.lab=0.8)
          axis(1,at = seq(0,1,length.out = ncol(signature_data_matrix)),labels = colnames(signature_data_matrix),las=2,cex.lab=0.8)
          draw_legend(col,1.25,1.3,0,1)
          
          #6 some correlation plots
          #pos <- which(max(abs(res.cor_triangular))==abs(res.cor_triangular),arr.ind = TRUE)
          vals <- res.cor_triangular[order(abs(res.cor_triangular),decreasing = TRUE)]
          for (j in 1:(nplots-5)){
            pos <- which(vals[j]==res.cor_triangular,arr.ind = TRUE)
            mainpar <- paste0("Exposures across bootstraps, n=",nboot,"\nspearman correlation ",sprintf("%.2f",vals[j]))
            plot(res$samples_list[[i]][pos[1],],res$samples_list[[i]][pos[2],],
                 xlab = colnames(signature_data_matrix)[pos[1]],
                 ylab = colnames(signature_data_matrix)[pos[2]],
                 # ylim = c(0,max(res$samples_list[[i]][pos[2],])),
                 # xlim = c(0,max(res$samples_list[[i]][pos[1],]))
                 main=mainpar,col="blue",pch = 16)
            
          }
          #sig.pca <- prcomp(t(res$samples_list[[i]]),center = TRUE,scale. = TRUE)
        }
      }
    }
    dev.off()
  }
}

#compute average correlations
average_cor <- matrix(0,nrow = 10,ncol = 10)
for (b in 1:length(res.cor_triangular_list)){
  average_cor <- average_cor + res.cor_triangular_list[[b]]
}
average_cor <- average_cor/length(res.cor_triangular_list)

jpeg(filename = paste0(resultdir,"fitAnalysis_part2/sigfit_bootstrap_expCor_average_all.jpg"),
     width = 800,
     height = 600,
     res = 150)
par(mar=c(6,8,5,6))
par(xpd=FALSE)
col<- colorRampPalette(c("blue", "white", "red"))(51)
image(average_cor,col = col,zlim = c(-1,1), axes=F,main=paste0("Exposures Correlation (spearman)\nAverage across ",n*ngenomes," simulated samples"))
extrabit <- 1/(ncol(signature_data_matrix)-1)/2
abline(h=seq(0-extrabit,1+extrabit,length.out = ncol(signature_data_matrix)+1),col="grey",lty=2)
abline(v=seq(0-extrabit,1+extrabit,length.out = ncol(signature_data_matrix)+1),col="grey",lty=2)
axis(2,at = seq(0,1,length.out = ncol(signature_data_matrix)),labels = colnames(signature_data_matrix),las=1,cex.lab=0.8)
axis(1,at = seq(0,1,length.out = ncol(signature_data_matrix)),labels = colnames(signature_data_matrix),las=2,cex.lab=0.8)
draw_legend(col,1.25,1.3,0,1)
dev.off()

