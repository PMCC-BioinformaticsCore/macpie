
library(RUVSeq)
library(reshape2)
library(gridExtra)
library(rlang)
library(edgeR)
library(RColorBrewer)
library(ggrepel)
# suppressMessages(library(factoextra))
# suppressMessages(library(rlang))

####################### limma voom ######################

limma_voom<-function(coef=1,seurat_obj,min_cpm=10, design.matrix, contrast.matrix, voom_plot=T, treatment_conc=treatment_conc){
  counts<-seurat_obj@assays$RNA$counts
  colnames(counts)<-seurat_obj$Well_ID
  #factor
  columns<-as.factor(seurat_obj$Column)
  rows<-as.factor(seurat_obj$Row)
  
  counts <- as.matrix(counts)
  keep <- rowSums(cpm(counts) > min_cpm )>= min(unique(table(treatment_conc)))
  #over-write keep
  # keep<-filterByExpr(counts,design.matrix)
  raw_counts <- counts[keep, ]
  
  
  set<-DGEList(raw_counts)
  set$samples$treatment_conc<-treatment_conc
  set$samples$columns<-columns
  set$samples$rows<-rows
  
  set<-calcNormFactors(set, method = "TMM")
  
  # plotMDS(cpm(raw_counts, log=TRUE), labels = set$samples$treatment_conc)
  set_voom<-voom(set, design.matrix, plot=voom_plot)
  
  fit<-lmFit(set_voom, design = design.matrix)
  
  fit2<-contrasts.fit(fit, contrast.matrix)
  fit3<-eBayes(fit2, robust = TRUE)
  # plotSA(fit3, main="Final:limma-voom mean variance trend on residual variances")
  
  res<-topTable(fit3, coef =coef, sort.by = "P", n=Inf, adjust.method = "BH")
  
  return(res)
}


####################### PCA ######################

PCAplot<-function(count_matrix, feature, treatment, center=TRUE ,scale=FALSE, extraPC=FALSE){
  if (any(count_matrix<0)){

    sv<-svd(as.matrix(count_matrix))
    
    treatment <- as.factor(treatment)
    feature<-as.factor(feature)
    
    p<-ggplot(data = as.data.frame(sv$u), aes(x = sv$u[,1], y = sv$u[,2], colour = feature, 
                                              shape = treatment))+scale_shape_manual(values = c(0:nlevels(as.factor(treatment)))) +
      geom_point() + xlab(paste0("PC1: ", round(sv$d[1]^2/sum(sv$d^2)*100, digits = 3), "%")) +
      ylab(paste0("PC2: ", round(as.numeric(sv$d[2]^2/sum(sv$d^2)*100), digits = 3), "%")) +
      scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(nlevels(as.factor(feature)))) +
      theme_bw() + 
      labs(colour = "Feature", shape = "Treatment")
    p
    
    
  }else{
    count_matrix<-apply(log(count_matrix + 1), 1, function(y) scale(y, center = TRUE, scale = FALSE))
    sv<-svd(as.matrix(count_matrix))
    
    treatment <- as.factor(treatment)
    feature<-as.factor(feature)
    
    p<-ggplot(data = as.data.frame(sv$u), aes(x = sv$u[,1], y = sv$u[,2], colour = feature, 
                                               shape = treatment))+scale_shape_manual(values = c(0:20,35:38)) +
      geom_point() + xlab(paste0("PC1: ", round(sv$d[1]^2/sum(sv$d^2)*100, digits = 3), "%")) +
      ylab(paste0("PC2: ", round(as.numeric(sv$d[2]^2/sum(sv$d^2)*100), digits = 3), "%")) +
      scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(nlevels(as.factor(feature)))) +
      # scale_color_gradient(low="blue", high="red")+ 
      # for showing gradient for read counts?
      theme_bw()
    
    

    p1<-ggplot(data = as.data.frame(sv$u), aes(x = sv$u[,1], y = sv$u[,2], colour = feature, 
                                              shape = treatment))+scale_shape_manual(values = c(0:20,35:38)) +
      geom_point() + xlab(paste0("PC1: ", round(sv$d[1]^2/sum(sv$d^2)*100, digits = 3), "%")) +
      ylab(paste0("PC2: ", round(as.numeric(sv$d[2]^2/sum(sv$d^2)*100), digits = 3), "%")) +
      scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(nlevels(as.factor(feature)))) +
      # scale_color_gradient(low="blue", high="red")+ 
      # for showing gradient for read counts?
      theme_bw()+NoLegend()
    
    p2<-ggplot(data = as.data.frame(sv$u), aes(x = sv$u[,3], y = sv$u[,4], colour = feature, 
                                              shape = treatment))+scale_shape_manual(values = c(0:20,35:38)) +
      geom_point() + xlab(paste0("PC3: ", round(sv$d[3]^2/sum(sv$d^2)*100, digits = 3), "%")) +
      ylab(paste0("PC4: ", round(as.numeric(sv$d[4]^2/sum(sv$d^2)*100), digits = 3), "%")) +
      scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(nlevels(as.factor(feature)))) +
      # scale_color_gradient(low="blue", high="red")+ 
      # for showing gradient for read counts?
      theme_bw() + 
      labs(colour = "Feature", shape = "Treatment")
    
    if (extraPC==TRUE){
      p1+p2
    }else{
      p
    }
 
    
  }
  
}

####################### RLE ######################


RLEplot<-function(count_matrix,ID,feature){
  count_matrix<-log2(count_matrix+1)
  count_matrix<-as.matrix(count_matrix)
  colnames(count_matrix)<-as.factor(ID)
  rle<-count_matrix-Biobase::rowMedians(count_matrix)
  feature<-as.factor(feature)
  rle<-rle[, sort.list(feature)]
  rledf <- t(apply(rle, 2, function(x) grDevices::boxplot.stats(x)$stats))
  rledf <- as.data.frame(rledf)
  colnames(rledf) <- c("ymin", "lower", "middle", "upper", "ymax")
  rledf$feature<-as.factor(feature[sort.list(feature)])
  rledf$sample<-rownames(rledf)
  rledf$sample<-reorder(as.factor(rledf$sample), unique(rledf$sample))#updated
  
  p<-ggplot(rledf, aes(x=sample, fill=feature))+geom_boxplot(aes(ymin=ymin, lower=lower, middle=middle, upper=upper, ymax=ymax), stat="identity")+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_x_discrete(limits=rledf$sample)+
    scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(nlevels(as.factor(feature))))+
    geom_hline(yintercept = 0, linetype="dotted", col="red", size=1)+ylab("RLE")
  p
}




#################### Volcano plots ####################


plot_volcano<-function(top, FDR_cutoff=0.05){
  top$diffexpressed<-"NO"
  top[top$logFC>=1 & top$FDR<=FDR_cutoff,]$diffexpressed<-"UP"
  top[top$logFC<=-1 & top$FDR<=FDR_cutoff,]$diffexpressed<-"DOWN"

  top$labelgenes<-""
  top[top$diffexpressed!="NO",]$labelgenes<-rownames(top[top$diffexpressed!="NO",])

  ggplot(top, aes(x=logFC, y=-log10(FDR), col=diffexpressed, label=labelgenes))+geom_point()+theme_classic()+
    geom_text_repel(min.segment.length=5)+scale_color_manual(values=c("blue", "black", "red"))+
    geom_vline(xintercept=c(-1, 1), col="red") + geom_hline(yintercept=-log10(0.05), col="red")
}





plot_volcano_lv_new<-function(top, FDR_cutoff=0.05, title="", min.segment.length=5, size=5){
  top$diffexpressed<-"NO"
  if(dim(top[top$logFC>=1 & top$adj.P.Val<=FDR_cutoff,])[1]>0){
    top[top$logFC>=1 & top$adj.P.Val<=FDR_cutoff,]$diffexpressed<-"UP"
  }
  if(dim(top[top$logFC<=-1 & top$adj.P.Val<=FDR_cutoff,])[1]>0){
    top[top$logFC<=-1 & top$adj.P.Val<=FDR_cutoff,]$diffexpressed<-"DOWN"
  }
  
  top$labelgenes<-""
  top[top$diffexpressed!="NO",]$labelgenes<-rownames(top[top$diffexpressed!="NO",])  
  
  
  
  if(length(unique(top$diffexpressed))==3){
    plot(ggplot(top, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=labelgenes))+geom_point()+theme_classic()+
           geom_text_repel(min.segment.length=min.segment.length, show.legend = F, size=size)+scale_color_manual(values=c("blue", "black", "red"))+
           geom_vline(xintercept=c(-1, 1), col="red") + geom_hline(yintercept=-log10(FDR_cutoff), col="red")+ggtitle(title))
  }
  
  else if(length(unique(top$diffexpressed))==2 & length(intersect(unique(top$diffexpressed), c("UP","NO")))==2){
    plot(ggplot(top, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=labelgenes))+geom_point()+theme_classic()+
           geom_text_repel(min.segment.length=min.segment.length, show.legend = F, size=size)+scale_color_manual(values=c("black", "red"))+
           geom_vline(xintercept=c(-1, 1), col="red") + geom_hline(yintercept=-log10(FDR_cutoff), col="red")+ggtitle(title))
    
  }
  
  else if(length(unique(top$diffexpressed))==2 & length(intersect(unique(top$diffexpressed), c("DOWN","NO")))==2){
    plot(ggplot(top, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=labelgenes))+geom_point()+theme_classic()+
           geom_text_repel(min.segment.length=min.segment.length, show.legend = F, size=size)+scale_color_manual(values=c("blue","black"))+
           geom_vline(xintercept=c(-1, 1), col="red") + geom_hline(yintercept=-log10(FDR_cutoff), col="red")+ggtitle(title))
  }
  
  else if(length(unique(top$diffexpressed))==2 & length(intersect(unique(top$diffexpressed), c("DOWN","UP")))==2){
    plot(ggplot(top, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=labelgenes))+geom_point()+theme_classic()+
           geom_text_repel(min.segment.length=min.segment.length, show.legend = F, size=size)+scale_color_manual(values=c("blue","red"))+
           geom_vline(xintercept=c(-1, 1), col="red") + geom_hline(yintercept=-log10(FDR_cutoff), col="red")+ggtitle(title))
  }  
  
  
  else{
    plot(ggplot(top, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=labelgenes))+geom_point()+theme_classic()+
           geom_text_repel(min.segment.length=min.segment.length, show.legend = F, size=size)+scale_color_manual(values=c("black"))+
           geom_vline(xintercept=c(-1, 1), col="red") + geom_hline(yintercept=-log10(FDR_cutoff), col="red")+ggtitle(title))
  }
  
  
  
}  




####################### RUVg ######################
RUVgDetails <- function(kfactor, mac, empirical, treatment1, treatment2) {
  treatment1 <- grep(treatment1, paste0(mac$Treatment_1, ".", mac$Concentration_1))
  treatment2 <- grep(treatment2, mac$Treatment_1)
  treatment <- make.names(as.factor(mac$Treatment_1[c(treatment1, treatment2)]))
  
  counts <- mac@assays$RNA$counts[, c(treatment1, treatment2)]
  colnames(counts) <- mac$Well_ID[c(treatment1, treatment2)]
  
  counts <- as.matrix(counts)
  keep <- rowSums(cpm(counts) > 10) >= 2
  raw_counts <- counts[keep, ]
  set <- newSeqExpressionSet(as.matrix(raw_counts), phenoData = data.frame(treatment,
                                                                           row.names = colnames(counts)))
  
  RUVset <- RUVg(set, empirical, kfactor)
  
  colors <- brewer.pal(8, "Set2")
  col <- mac$Column[c(treatment1, treatment2)]
  
  plotRLE(RUVset, outline = FALSE, ylim = c(-4, 4), col = colors[col]) + title(paste0("Adjusting RUV; k=",
                                                                                      kfactor))
  
  object <- RUVset@assayData$normalizedCounts
  Y <- apply(log(object + 1), 1, function(y) scale(y, center = TRUE, scale = FALSE))
  rownames(Y) <- c(colnames(object))
  prcomp_pca <- prcomp(Y, center = FALSE, scale. = FALSE)
  pca_summary <- summary(prcomp_pca)
  plot(prcomp_pca$x, type = "n", xlab = paste0("PC1: ", pca_summary$importance[2, 1] * 100, "%"), ylab = paste0("PC2: ", pca_summary$importance[2, 2] * 100,
                                                                                                                "%")) + title(paste0("Adjusting RUV; k=", kfactor))
  text(prcomp_pca$x, labels = rownames(prcomp_pca$x))
  
  if (kfactor > nrow(prcomp_pca$x)) {
    stop()
  } else {
    kmeans_res <- eclust(prcomp_pca$x, "kmeans", hc_metric = "eucliden", k = 2)
    # kmeans_res$cluster
    
    # within group identity automate number of samples each group in the
    # future
    if (kmeans_res$cluster[1] == kmeans_res$cluster[2]) {
      Group1 <- "T"
    } else {
      Group1 <- "F"
    }
    
    if (kmeans_res$cluster[3] == kmeans_res$cluster[4]) {
      Group2 <- "T"
    } else {
      Group2 <- "F"
    }
    
    k_groups <- c(kfactor, Group1, Group2)
    names(k_groups) <- c(paste0("k="), paste0(names(kmeans_res$cluster[1]), "-",
                                              names(kmeans_res$cluster[2])), paste0(names(kmeans_res$cluster[3]), "-",
                                                                                    names(kmeans_res$cluster[4])))
    # return(k_groups)
    return(list(set = RUVset, k_details = k_groups))
  }
  
}


initial_DE<-function(mac,treatment1, treatment2, topgenes_keep=5000){
  treatment1<-grep(treatment1, paste0(mac$Treatment_1,".", mac$Concentration_1))
  treatment2<-grep(treatment2, mac$Treatment_1)
  treatment <- make.names(as.factor(mac$Treatment_1[c(treatment1, treatment2)]))
  
  counts <- mac@assays$RNA$counts[, c(treatment1, treatment2)]
  colnames(counts) <- mac$Well_ID[c(treatment1, treatment2)]
  
  counts <- as.matrix(counts)
  keep <- rowSums(cpm(counts) > 10) >= 2
  raw_counts <- counts[keep, ]
  set <- newSeqExpressionSet(as.matrix(raw_counts), phenoData = data.frame(treatment,
                                                                           row.names = colnames(counts)))
  design <- model.matrix(~treatment, data = pData(set))
  
  y <- DGEList(counts = counts(set), group = treatment)
  y <- calcNormFactors(y, method = "upperquartile")
  # glm
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  plotBCV(y)+title("before RUV")
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit)
  top <- topTags(lrt, n = nrow(set))$table
  
  plot(ggplot(top, aes(x = PValue)) + geom_histogram()+ggtitle(label = "before RUV")+theme_bw())
  colors <- c(brewer.pal(8, "Set2"),brewer.pal(12, "Set3"), brewer.pal(9, "Set1"))
  col <- mac$Column[c(treatment1, treatment2)]
  plotRLE(set, outline = FALSE, ylim = c(-4, 4), col = colors[col])+title("before RUV")
  plotPCA(set, col = colors[col])+title("before RUV")
  
  #plot PCA in kmeans
  object <- set@assayData$counts
  Y <- apply(log(object + 1), 1, function(y) scale(y, center = TRUE, scale = FALSE))
  rownames(Y) <- c(colnames(object))
  prcomp_pca <- prcomp(Y, center = FALSE, scale. = FALSE)
  
  
  #Kmeans
  #if the original PCA looks fine, don't need RUV -> direct to DEGs
  #else RUV 
  kmeans_res <- eclust(prcomp_pca$x, "kmeans", hc_metric = "eucliden", k = 2)
  
  
  # kmeans_res$cluster
  
  # within group identity automate number of samples each group in the
  # future
  if (kmeans_res$cluster[1] == kmeans_res$cluster[2]) {
    Group1 <- "T"
  } else {
    Group1 <- "F"
  }
  
  if (kmeans_res$cluster[3] == kmeans_res$cluster[4]) {
    Group2 <- "T"
  } else {
    Group2 <- "F"
  }
  
  k_groups <- c("original", Group1, Group2)
  names(k_groups) <- c(paste0("original"), paste0(names(kmeans_res$cluster[1]), "-",
                                                  names(kmeans_res$cluster[2])), paste0(names(kmeans_res$cluster[3]), "-",
                                                                                        names(kmeans_res$cluster[4])))
  
  k_groups<-as.data.frame(k_groups)
  
  if(k_groups[2,]=="T" & k_groups[3,]=="T"){
    print("Original p value distribution, PCA and Kmeans groups make sense.")
  }else{
    print("Considering RUV or other approaches.")
  }
  
  empirical_genes <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:topgenes_keep]))]
  return(list(set=set, empirical_genes=empirical_genes, k_details=k_groups))
  
  
}  


k_summary<-function(list_set_k){
  klist_df<-bind_rows(lapply(list_set_k, function(x) {
    x$k_details
  }))
  colnames(klist_df)<-c("k","g1","g2")
  k_opt<-klist_df%>%filter(klist_df[,2]=="T" & klist_df[,3]=="T")%>%filter(k==min(k))
  return(k_opt)
}





DE_RUVg<-function(mac, treatment1, treatment2, design_matrix){
  treatment1<-grep(treatment1, paste0(mac$Treatment_1,".", mac$Concentration_1))
  treatment2<-grep(treatment2, mac$Treatment_1)
  treatment <- make.names(as.factor(mac$Treatment_1[c(treatment1, treatment2)]))
  
  counts <- mac@assays$RNA$counts[, c(treatment1, treatment2)]
  colnames(counts) <- mac$Well_ID[c(treatment1, treatment2)]
  
  counts <- as.matrix(counts)
  keep <- rowSums(cpm(counts) > 10) >= 2
  raw_counts <- counts[keep, ]
  set <- newSeqExpressionSet(as.matrix(raw_counts), phenoData = data.frame(treatment,
                                                                           row.names = colnames(counts)))
  
  y <- DGEList(counts = counts(set), group = treatment)
  y <- calcNormFactors(y, method = "upperquartile")
  # glm
  y <- estimateGLMCommonDisp(y, design_matrix)
  y <- estimateGLMTagwiseDisp(y, design_matrix)
  plotBCV(y)+title("after RUVg")
  fit <- glmFit(y, design_matrix)
  lrt <- glmLRT(fit)
  top <- topTags(lrt, n = nrow(set))$table
  
  plot(ggplot(top, aes(x = PValue)) + geom_histogram()+ggtitle(label = "RUVg then DE, Pvalue")+theme_bw())
  plotMD(lrt)
  abline(h=c(-1,-1), col="blue")
  summary(decideTests.DGELRT(lrt))
  
  return(top)
  
}




