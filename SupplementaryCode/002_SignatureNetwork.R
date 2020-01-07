source("../lib/SignatureExtractionLib.R")
source("../lib/SignatureFitLib.R")

redo_fit <- FALSE
redo_filter <- FALSE
replot <- FALSE

type_of_sigs <- "rearr"
minexp <- 0.15
cos.sim.thr <- 0.89
min_increase_cossim <- 0.02

#version
version <- "_v2" #second version,

#how to order signature fits
order_fits_by <- "cos.sim"
#order_fits_by <- "KLD"

how_many_sigs_in_combination <- 1:2
nboots <- 20
nparallel <- 6


#functions----
#function to draw a legend for the heatmap of the correlation matrix
draw_legend <- function(col,xl,xr,yb,yt){
  par(xpd=TRUE)
  # par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),mar = c(0, 0, 0, 0), new = TRUE)
  # plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n",xlim = c(0,1),ylim = c(0,1))
  #plot(NULL,xlim = c(0,1),ylim = c(0,1),add=TRUE)
  rect(xl,yb,xr,yt)
  rect(
    xl,
    head(seq(yb,yt,(yt-yb)/length(col)),-1),
    xr,
    tail(seq(yb,yt,(yt-yb)/length(col)),-1),
    col=col,border = NA
  )
  text(x = 1.2, y = yt,labels = "1")
  text(x = 1.2, y = (yt-yb)/2,labels = ".5")
  text(x = 1.2, y = yb,labels = "0")
}
#functions end----



# args = commandArgs(trailingOnly=TRUE)
# group <- args[1]
# type_of_sigs <- args[2]

choice_table <- list()
choice_table[["subs"]] <- read.table("../data/2018_02_18_fullExtraction_pancan_signatures_subs.txt",
                                     header = TRUE,as.is = TRUE,check.names = FALSE,sep = "\t")
choice_table[["rearr"]] <- read.table("../data/2018_02_18_fullExtraction_pancan_signatures_rearr.txt",
                                      header = TRUE,as.is = TRUE,check.names = FALSE,sep = "\t")

#i <- which(choice_table[[type_of_sigs]]$Organ==group)


plotsdir <- paste0("../results/signatures_fitDiscovery",version,"/max_nsig_",max(how_many_sigs_in_combination),"_minCS",cos.sim.thr,"_minexp",minexp,"_minCSincr",min_increase_cossim,"_order",order_fits_by,"/",type_of_sigs,"/")
dir.create(plotsdir,showWarnings = FALSE,recursive = TRUE)


full_graph_table <- NULL
closestCOSMIC <- NULL

for(i in 1:nrow(choice_table[[type_of_sigs]])){
  group <- choice_table[[type_of_sigs]]$Organ[i]
  message("Group Organ to Fit: ",group)

  ns <- choice_table[[type_of_sigs]]$nsig[choice_table[[type_of_sigs]]$Organ==group]
  tag <- ""
  if (choice_table[[type_of_sigs]]$scale[choice_table[[type_of_sigs]]$Organ==group]) tag <- paste0(tag,"_scale")
  signatures <- read.table(paste0("../results/",type_of_sigs,"/extraction/brunet/extraction_",group,"_",type_of_sigs,"_brunet",tag,"/round_1/sig_",ns,"/Sigs_plot_",group,"_ns",ns,"_nboots",nboots,".tsv"),
                           sep = "\t",header = TRUE,as.is = TRUE,check.names = FALSE)
  if (type_of_sigs=="subs") signatures <- sortCatalogue(signatures)
  colnames(signatures) <- paste0(rep(group,ncol(signatures)),"_",seq(1,ncol(signatures)))
  if(type_of_sigs=="subs"){
    cc <- findClosestCOSMIC30_withSimilarity(signatures)
    for (c in 1:length(cc$cosmic)){
      closestCOSMIC <- rbind(closestCOSMIC,data.frame(sig=colnames(signatures)[c],
                                                      cosmic=cc$cosmic[c],
                                                      cos.sim=cc$cos.sim[c]))
    }
  }else if(type_of_sigs=="rearr"){
    cc <- findClosestRearrSigsBreast560_withSimilarity(signatures)
    for (c in 1:length(cc$RS.Breast560)){
      closestCOSMIC <- rbind(closestCOSMIC,data.frame(sig=colnames(signatures)[c],
                                                      cosmic=cc$RS.Breast560[c],
                                                      cos.sim=cc$cos.sim[c]))
    }
  }
  
  r_file <- paste0(plotsdir,"sigFitDiscovery_",group,"_cstr",cos.sim.thr,".Rdata")
  r_file_filtered <- paste0(plotsdir,"sigFitDiscovery_",group,"_cstr",cos.sim.thr,"_minexp",minexp,"_filtered.Rdata")
  
  if(file.exists(r_file) & !redo_fit){
    load(r_file)
    message("Results loaded from file")
  }else{
    all_results <- list()
    count <- 0
    for(h in how_many_sigs_in_combination){
      all_results[[h]] <- list()
      for(j in 1:nrow(choice_table[[type_of_sigs]])){
        #use sigs from other organs j to fit these sigs in organ i
        if(i!=j & choice_table[[type_of_sigs]][j,"nsig"]>=h){ 
          group2 <- choice_table[[type_of_sigs]][j,"Organ"]
          message("Organ: ",group2,", using ",h," signature(s)")
          ns2 <- choice_table[[type_of_sigs]]$nsig[choice_table[[type_of_sigs]]$Organ==group2]
          tag2 <- ""
          if (choice_table[[type_of_sigs]]$scale[choice_table[[type_of_sigs]]$Organ==group2]) tag2 <- paste0(tag2,"_scale")
          signatures2 <- read.table(paste0("../results/",type_of_sigs,"/extraction/brunet/extraction_",group2,"_",type_of_sigs,"_brunet",tag2,"/round_1/sig_",ns2,"/Sigs_plot_",group2,"_ns",ns2,"_nboots",nboots,".tsv"),
                                   sep = "\t",header = TRUE,as.is = TRUE,check.names = FALSE)
          if (type_of_sigs=="subs") signatures2 <- sortCatalogue(signatures2)
          colnames(signatures2) <- paste0(rep(group2,ncol(signatures2)),"_",seq(1,ncol(signatures2)))
          
          comb_of_3 <- combn(1:ns2,h)
          all_res <- list()
          cosSims <- list()
          KLD_list <- list()
          for(z in 1:ncol(signatures)) cosSims[[z]] <- array(NA,dim = ncol(comb_of_3))
          for(z in 1:ncol(signatures)) KLD_list[[z]] <- array(NA,dim = ncol(comb_of_3))
          for(z in 1:ncol(signatures)) all_res[[z]] <- array(NA,dim = c(h,ncol(comb_of_3)))
          for (k in 1:ncol(comb_of_3)){
            #message("k=",k," of ",ncol(comb_of_3))
            sigs_to_use <- signatures2[,comb_of_3[,k]]
            res <- SignatureFit(signatures,sigs_to_use,doRound = FALSE,verbose = FALSE)
            reconstructed <- as.matrix(sigs_to_use) %*% res
            for (z in 1:ncol(signatures)){
              cosSims[[z]][k] <- cos.sim(signatures[,z],reconstructed[,z])
              KLD_list[[z]][k] <- KLD(signatures[,z],reconstructed[,z])
              all_res[[z]][,k] <- res[,z]
            }
          }
          
          #get the best for each signature
          report_tables <- list()
          for(z in 1:ncol(signatures)) {
            select_comb <- cosSims[[z]]>=cos.sim.thr
            report_tables[[z]] <- matrix(0,nrow = sum(select_comb),ncol = ns2)
            if(sum(select_comb)>0){
              if(order_fits_by=="cos.sim"){
                select_pos_ordered <- order(cosSims[[z]][select_comb],decreasing = TRUE)
              }else if (order_fits_by=="KLD"){
                select_pos_ordered <- order(KLD_list[[z]][select_comb],decreasing = FALSE)
              }else{
                stop("I don't know how to order the signature fits, set order_fits_by to cos.sim or KLD")
              }
              colnames(report_tables[[z]]) = colnames(signatures2)
              combs <- comb_of_3[,select_comb,drop=FALSE]
              exps <- all_res[[z]][,select_comb,drop=FALSE]
              for(k in select_pos_ordered){
                report_tables[[z]][k,combs[,k]] <- exps[,k]
              }
              report_tables[[z]] <- cbind(report_tables[[z]][select_pos_ordered,,drop=FALSE],cosSims[[z]][select_comb][select_pos_ordered],KLD_list[[z]][select_comb][select_pos_ordered])
              colnames(report_tables[[z]])[ns2+1] <- "cos.sim"
              colnames(report_tables[[z]])[ns2+2] <- "KLD"
              #count the tables
              count <- count + 1
            }
          }
          all_results[[h]][[group2]] <- report_tables
        }
      }
    }
    save(file = r_file,all_results,count)
  }
  
  if(file.exists(r_file_filtered) & !redo_filter){
    load(r_file_filtered)
    message("Filtered Results loaded from file")
  }else{
    #---processing to select best/representative that I can use to plot graphs---
    filtered_results <- list()
    if(count>0){
      count_filtered <- 0
      for(z in 1:ncol(signatures)){
        filtered_results[[z]] <- list()
        for(j in 1:nrow(choice_table[[type_of_sigs]])){
          group2 <- choice_table[[type_of_sigs]][j,"Organ"]
          ns2 <- choice_table[[type_of_sigs]]$nsig[choice_table[[type_of_sigs]]$Organ==group2]
          filtered_results[[z]][[group2]] <- matrix(0,nrow = 0,ncol = ns2+2)
          for(h in how_many_sigs_in_combination){
            #use sigs from other organs j to fit these sigs in organ i
            if(i!=j & choice_table[[type_of_sigs]][j,"nsig"]>=h){
              current_table <- all_results[[h]][[group2]][[z]]
              #apply filter on exposures, need at least minexp
              if(nrow(current_table)>0){
                current_table <- current_table[apply(current_table[,1:ns2,drop=FALSE],1,function(x) all(x[which(x > 0)]>=minexp)),,drop=FALSE]
              }
              if(nrow(current_table)>0){
                #we can use either best cos.sim or best KLD
                pos_of_best <- 1 #already ordered by cos sim
                if(nrow(filtered_results[[z]][[group2]])>0){
                  #there is a competitor that used less signatures! Compare cosine similarities
                  if(current_table[pos_of_best,"cos.sim"] >= filtered_results[[z]][[group2]][1,"cos.sim"] + min_increase_cossim){
                    #improvement of cosine similarity of at least 1%
                    filtered_results[[z]][[group2]] <- current_table[pos_of_best,,drop=FALSE]
                  }
                }else{
                  #no better combiantion was found with less signatures, so just put the best here:
                  filtered_results[[z]][[group2]] <- current_table[pos_of_best,,drop=FALSE]
                  count_filtered <- count_filtered + 1
                }
              }
            }
          }
        }
      }
    }
    save(file = r_file_filtered,filtered_results,count_filtered)
  }
  
  #---plotting---
  
  if(count>0 & replot){
    jpeg(filename = paste0(plotsdir,"sigFitDiscovery_",group,"_cstr",cos.sim.thr,"_minexp",minexp,".jpg"),
         width = 640*(2+max(how_many_sigs_in_combination)),
         height = 480*count,
         res = 150)
    par(mfrow=c(count,2+max(how_many_sigs_in_combination)))
    for(z in 1:ncol(signatures)){
      for(h in how_many_sigs_in_combination){
        for(j in 1:nrow(choice_table[[type_of_sigs]])){
          #use sigs from other organs j to fit these sigs in organ i
          if(i!=j & choice_table[[type_of_sigs]][j,"nsig"]>=h){ 
            group2 <- choice_table[[type_of_sigs]][j,"Organ"]
            
            ns2 <- choice_table[[type_of_sigs]]$nsig[choice_table[[type_of_sigs]]$Organ==group2]
            tag2 <- ""
            if (choice_table[[type_of_sigs]]$scale[choice_table[[type_of_sigs]]$Organ==group2]) tag2 <- paste0(tag2,"_scale")
            signatures2 <- read.table(paste0("../results/",type_of_sigs,"/extraction/brunet/extraction_",group2,"_",type_of_sigs,"_brunet",tag2,"/round_1/sig_",ns2,"/Sigs_plot_",group2,"_ns",ns2,"_nboots",nboots,".tsv"),
                                      sep = "\t",header = TRUE,as.is = TRUE,check.names = FALSE)
            if (type_of_sigs=="subs") signatures2 <- sortCatalogue(signatures2)
            colnames(signatures2) <- paste0(rep(group2,ncol(signatures2)),"_",seq(1,ncol(signatures2)))
            
            current_table <- all_results[[h]][[group2]][[z]]
            #apply filter on exposures, need at least minexp
            if(nrow(current_table)>0){
              current_table <- current_table[apply(current_table[,1:ns2,drop=FALSE],1,function(x) all(x[which(x > 0)]>=minexp)),,drop=FALSE]
            }
            #1. Draw heatmap of the table
            if(nrow(current_table)>0){
              #scale KLD
              current_table[,"KLD"] <- current_table[,"KLD"]/max(current_table[,"KLD"])
              #we have a row to plot!
              par(mar=c(3,8,3,6))
              par(xpd=FALSE)
              col<- colorRampPalette(c("white", "blue"))(51)
              image(current_table,col = col,zlim = c(0,1), axes=F,main="Proportion of Signatures")
              extrabit_v <- 1/(ncol(current_table)-1)/2
              abline(h=seq(0-extrabit_v,1+extrabit_v,length.out = ncol(current_table)+1),col="grey",lty=2)
              #extrabit_h <- 1/(nrow(current_table)-1)/2
              #abline(v=seq(0-extrabit_h,1+extrabit_h,length.out = nrow(current_table)+1),col="grey",lty=2)
              axis(2,at = seq(0,1,length.out = ncol(current_table)),labels = colnames(current_table),las=1,cex.lab=0.8)
              #axis(1,at = seq(0,1,length.out = nrow(current_table)),labels = rep("",nrow(current_table)),las=2,cex.lab=0.8,tick = FALSE)
              par(xpd=TRUE)
              text(labels = "combinations",x = 0.5,y = - 0.1)
              draw_legend(col,1.25,1.3,0,1)
              
              #find the three most relevant signatures
              top3 <- order(apply(current_table[,1:ns2,drop=FALSE],2,sum),decreasing = TRUE)[1:h]
              
              #now draw the interesting signatures
              if(type_of_sigs=="subs"){
                #2 original
                plotSubsSignatures(signature_data_matrix = signatures[,z,drop=FALSE],mar=c(6,3,5,2),plot_sum = FALSE)
                #top3
                for(q in top3) plotSubsSignatures(signature_data_matrix = signatures2[,q,drop=FALSE],mar=c(6,3,5,2),plot_sum = FALSE)
                if(h<max(how_many_sigs_in_combination)){
                  for (q in 1:(max(how_many_sigs_in_combination)-h)){
                    #fake plot to fill the space
                    par(mar=c(0,0,0,0))
                    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
                  }
                }
              }else if(type_of_sigs=="rearr"){
                #2 original
                plotRearrSignatures(signature_data_matrix = signatures[,z,drop=FALSE],mar=c(12,3,5,2),plot_sum = FALSE)
                #top3
                for(q in top3) plotRearrSignatures(signature_data_matrix = signatures2[,q,drop=FALSE],mar=c(12,3,5,2),plot_sum = FALSE)
                if(h<max(how_many_sigs_in_combination)){
                  for (q in 1:(max(how_many_sigs_in_combination)-h)){
                    #fake plot to fill the space
                    par(mar=c(0,0,0,0))
                    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
                  }
                }
                
              }else if(type_of_sigs=="generic"){
                #2 original
                plotGenericSignatures(signature_data_matrix = signatures[,z,drop=FALSE],mar=c(6,3,5,2),plot_sum = FALSE)
                #top3
                for(q in top3) plotGenericSignatures(signature_data_matrix = signatures2[,q,drop=FALSE],mar=c(6,3,5,2),plot_sum = FALSE)
                if(h<max(how_many_sigs_in_combination)){
                  for (q in 1:(max(how_many_sigs_in_combination)-h)){
                    #fake plot to fill the space
                    par(mar=c(0,0,0,0))
                    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
                  }
                }
                  
                
              }
              
            }
            
          }
        }
      }
    }
    dev.off()
  }else{
    if(count==0) count_filtered <- 0
  }
  
  
  #other plot go!
  
  if(count_filtered>0 & replot){
    jpeg(filename = paste0(plotsdir,"sigFitDiscovery_",group,"_cstr",cos.sim.thr,"_minexp",minexp,"_filtered.jpg"),
         width = 640*(2+max(how_many_sigs_in_combination)),
         height = 480*count_filtered,
         res = 150)
    par(mfrow=c(count_filtered,2+max(how_many_sigs_in_combination)))
    
    for(z in 1:ncol(signatures)){
      for(j in 1:nrow(choice_table[[type_of_sigs]])){
        group2 <- choice_table[[type_of_sigs]][j,"Organ"]
        
        ns2 <- choice_table[[type_of_sigs]]$nsig[choice_table[[type_of_sigs]]$Organ==group2]
        tag2 <- ""
        if (choice_table[[type_of_sigs]]$scale[choice_table[[type_of_sigs]]$Organ==group2]) tag2 <- paste0(tag2,"_scale")
        signatures2 <- read.table(paste0("../results/",type_of_sigs,"/extraction/brunet/extraction_",group2,"_",type_of_sigs,"_brunet",tag2,"/round_1/sig_",ns2,"/Sigs_plot_",group2,"_ns",ns2,"_nboots",nboots,".tsv"),
                                  sep = "\t",header = TRUE,as.is = TRUE,check.names = FALSE)
        if (type_of_sigs=="subs") signatures2 <- sortCatalogue(signatures2)
        colnames(signatures2) <- paste0(rep(group2,ncol(signatures2)),"_",seq(1,ncol(signatures2)))
        
        current_table <- filtered_results[[z]][[group2]]
        if (nrow(current_table)>0){
          sigs_pos <- which(current_table[1,1:ns2]>0)
          
          # for(pp in sigs_pos){
          #   full_graph_table <- rbind(full_graph_table,data.frame(source=colnames(signatures)[z],
          #                                                         target=colnames(current_table)[pp],
          #                                                         value=current_table[1,pp]))
          # }
          
          #plot:
          #1. original signature
          #2. modelled signature
          #3. signatures used, up to three sigs
          
          reconstructed_sig <- as.matrix(signatures2) %*% current_table[,1:ns2]
          #now draw the interesting signatures
          if(type_of_sigs=="subs"){
            #1 original
            plotSubsSignatures(signature_data_matrix = signatures[,z,drop=FALSE],mar=c(6,3,5,2),plot_sum = FALSE,add_to_titles = paste0("Original"))
            #2 model
            plotSubsSignatures(signature_data_matrix = reconstructed_sig,mar=c(6,3,5,2),plot_sum = FALSE,add_to_titles = paste0("Model, cosSim=",sprintf("%.2f",current_table[,"cos.sim"])))
            #top3
            for(q in sigs_pos) plotSubsSignatures(signature_data_matrix = signatures2[,q,drop=FALSE],mar=c(6,3,5,2),plot_sum = FALSE,add_to_titles = sprintf("contrib.=%.2f",current_table[1,q]))
            if(length(sigs_pos)<max(how_many_sigs_in_combination)){
              for (q in 1:(max(how_many_sigs_in_combination)-length(sigs_pos))){
                #fake plot to fill the space
                par(mar=c(0,0,0,0))
                plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
              }
            }
          }else if(type_of_sigs=="rearr"){
            #1 original
            plotRearrSignatures(signature_data_matrix = signatures[,z,drop=FALSE],mar=c(12,3,5,2),plot_sum = FALSE,add_to_titles = paste0("Original"))
            #2 model
            plotRearrSignatures(signature_data_matrix = reconstructed_sig,mar=c(12,3,5,2),plot_sum = FALSE,add_to_titles = paste0("Model, cosSim=",sprintf("%.2f",current_table[,"cos.sim"])))
            #top3
            for(q in sigs_pos) plotRearrSignatures(signature_data_matrix = signatures2[,q,drop=FALSE],mar=c(12,3,5,2),plot_sum = FALSE,add_to_titles = sprintf("contrib.=%.2f",current_table[1,q]))
            if(length(sigs_pos)<max(how_many_sigs_in_combination)){
              for (q in 1:(max(how_many_sigs_in_combination)-length(sigs_pos))){
                #fake plot to fill the space
                par(mar=c(0,0,0,0))
                plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
              }
            }
          }else if(type_of_sigs=="generic"){
            #1 original
            plotGenericSignatures(signature_data_matrix = signatures[,z,drop=FALSE],mar=c(6,3,5,2),plot_sum = FALSE,add_to_titles = paste0("Original"))
            #2 model
            plotGenericSignatures(signature_data_matrix = reconstructed_sig,mar=c(6,3,5,2),plot_sum = FALSE,add_to_titles = paste0("Model, cosSim=",sprintf("%.2f",current_table[,"cos.sim"])))
            #top3
            for(q in sigs_pos) plotGenericSignatures(signature_data_matrix = signatures2[,q,drop=FALSE],mar=c(6,3,5,2),plot_sum = FALSE,add_to_titles = sprintf("contrib.=%.2f",current_table[1,q]))
            if(length(sigs_pos)<max(how_many_sigs_in_combination)){
              for (q in 1:(max(how_many_sigs_in_combination)-length(sigs_pos))){
                #fake plot to fill the space
                par(mar=c(0,0,0,0))
                plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
              }
            }
          }
          
        }
      }
    }
    dev.off()
  }
  
  #get graph
  if(count_filtered>0){
    for(z in 1:ncol(signatures)){
      for(j in 1:nrow(choice_table[[type_of_sigs]])){
        group2 <- choice_table[[type_of_sigs]][j,"Organ"]
        
        ns2 <- choice_table[[type_of_sigs]]$nsig[choice_table[[type_of_sigs]]$Organ==group2]
        # tag2 <- ""
        # if (choice_table[[type_of_sigs]]$scale[choice_table[[type_of_sigs]]$Organ==group2]) tag2 <- paste0(tag2,"_scale")
        # signatures2 <- read.table(paste0("../results/",type_of_sigs,"/extraction/brunet/extraction_",group2,"_",type_of_sigs,"_brunet",tag2,"/round_1/sig_",ns2,"/Sigs_plot_",group2,"_ns",ns2,"_nboots",nboots,".tsv"),
        #                           sep = "\t",header = TRUE,as.is = TRUE,check.names = FALSE)
        # if (type_of_sigs=="subs") signatures2 <- sortCatalogue(signatures2)
        # colnames(signatures2) <- paste0(rep(group2,ncol(signatures2)),"_",seq(1,ncol(signatures2)))
        
        current_table <- filtered_results[[z]][[group2]]
        if (nrow(current_table)>0){
          sigs_pos <- which(current_table[1,1:ns2]>0)
          
          for(pp in sigs_pos){
            full_graph_table <- rbind(full_graph_table,data.frame(source=colnames(current_table)[pp],
                                                                  target=colnames(signatures)[z],
                                                                  value=current_table[1,pp]))
          }
        }
      }
    }
    
  }
}

full_graph_table$source <- gsub("\\.","\\_",full_graph_table$source)

write.table(full_graph_table,
            file=paste0(plotsdir,"sigFitDiscovery_graphD3js_tr",cos.sim.thr,"_minexp",minexp,".csv"),
            sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE)

nodes <- unique(union(full_graph_table$source,full_graph_table$target))
nodes_table <- closestCOSMIC[closestCOSMIC$sig %in% nodes,]
# nodes_table$colour <- NA
# #nodes_table$simpleid <- gsub("\\.", "\\\\\\\\\\\\\\\\\\\\.", nodes_table$sig)
# unique_cosmic <- unique(nodes_table$cosmic)
# for (c in 1:length(unique_cosmic)){
#   nodes_table[nodes_table$cosmic==unique_cosmic[c],"colour"] <- paste0("#",kelly_colors[c])
# }
write.table(nodes_table,
            file=paste0(plotsdir,"sigFitDiscovery_graphD3js_tr",cos.sim.thr,"_minexp",minexp,"_cosmic.csv"),
            sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE)

#colours
kelly_colors <- c('F2F3F4', '222222', 'F3C300', '875692', 'F38400', 'A1CAF1', 'BE0032', 
                  'C2B280', '848482', '008856', 'E68FAC', '0067A5', 'F99379', '604E97', 
                  'F6A600', 'B3446C', 'DCD300', '882D17', '8DB600', '654522', 'E25822', '2B3D26','CCCCCC','CCCCCC','CCCCCC')



#----look at main sigs
threshold_value <- 0.6

tmp_table <- full_graph_table[full_graph_table$value>=threshold_value,]
write.table(tmp_table,
            file=paste0(plotsdir,"sigFitDiscovery_graphD3js_tr",cos.sim.thr,"_minexp",minexp,"_value",threshold_value,".csv"),
            sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE)

nodes <- unique(union(tmp_table$source,tmp_table$target))
nodes_table <- closestCOSMIC[closestCOSMIC$sig %in% nodes,]
nodes_table$colour <- NA
#nodes_table$simpleid <- gsub("\\.", "\\\\\\\\\\\\\\\\\\\\.", nodes_table$sig)
unique_cosmic <- unique(nodes_table$cosmic)
for (c in 1:length(unique_cosmic)){
  nodes_table[nodes_table$cosmic==unique_cosmic[c],"colour"] <- paste0("#",kelly_colors[c])
}
write.table(nodes_table,
            file=paste0(plotsdir,"sigFitDiscovery_graphD3js_tr",cos.sim.thr,"_minexp",minexp,"_value",threshold_value,"_cosmic.csv"),
            sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE)

colour_table <- unique(nodes_table[,c("cosmic","colour")])
write.table(colour_table,
            file=paste0(plotsdir,"sigFitDiscovery_graphD3js_tr",cos.sim.thr,"_minexp",minexp,"_value",threshold_value,"_cosmic_colours.csv"),
            sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE)


#----look at low contributions
tmp_table <- full_graph_table[full_graph_table$value<threshold_value,]
write.table(tmp_table,
            file=paste0(plotsdir,"sigFitDiscovery_graphD3js_tr",cos.sim.thr,"_minexp",minexp,"_value",threshold_value,"less.csv"),
            sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE)

nodes <- unique(union(tmp_table$source,tmp_table$target))
nodes_table <- closestCOSMIC[closestCOSMIC$sig %in% nodes,]
nodes_table$colour <- NA
#nodes_table$simpleid <- gsub("\\.", "\\\\\\\\\\\\\\\\\\\\.", nodes_table$sig)
unique_cosmic <- unique(nodes_table$cosmic)
for (c in 1:length(unique_cosmic)){
  nodes_table[nodes_table$cosmic==unique_cosmic[c],"colour"] <- paste0("#",kelly_colors[c])
}
write.table(nodes_table,
            file=paste0(plotsdir,"sigFitDiscovery_graphD3js_tr",cos.sim.thr,"_minexp",minexp,"_value",threshold_value,"less_cosmic.csv"),
            sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE)

colour_table <- unique(nodes_table[,c("cosmic","colour")])
write.table(colour_table,
            file=paste0(plotsdir,"sigFitDiscovery_graphD3js_tr",cos.sim.thr,"_minexp",minexp,"_value",threshold_value,"less_cosmic_colours.csv"),
            sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE)

# s <- 1
# for (g in names(all_results[[1]])) {print(g); print(all_results[[1]][[g]][s])}
# for (g in names(all_results[[2]])) {print(g); print(all_results[[2]][[g]][s])}
# for (g in names(all_results[[3]])) {print(g); print(all_results[[3]][[g]][s])}