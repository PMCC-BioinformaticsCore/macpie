suppressMessages(library(ggplot2))
library(SeuratObject)
library(gridExtra)

library(RColorBrewer)
library(ggrepel)
####################### Seurat's clustering ######################
clustering <- function(mac) {
  sct <- mac %>%
    PercentageFeatureSet(., pattern = "^MT-", col.name = "percent.mt") %>%
    PercentageFeatureSet(., pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",
                         col.name = "percent.ribo") %>%
    SCTransform(do.scale = T, return.only.var.genes = F, vars.to.regress = "percent.mt",
                verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(verbose = FALSE) %>%
    FindClusters(method = "igraph", verbose = FALSE)
  return(sct)
}

plot_clusters <- function(mac) {
  p1 <- FeaturePlot(mac, features = "nCount_RNA")
  p2 <- FeaturePlot(mac, features = "percent.mt")
  p3 <- FeaturePlot(mac, features = "percent.ribo")
  p4 <- DimPlot(mac, group.by = "Row")
  p5 <- DimPlot(mac, group.by = "Column")
  p6 <- DimPlot(mac, group.by = "Concentration_1")
  p7 <- DimPlot(mac, group.by = "seurat_clusters")
  p8 <- DimPlot(mac, group.by = "Cell_type")
  p10 <- DimPlot(mac, group.by = "Sample_type")
  p9 <- DimPlot(mac, group.by = "edge_wells", split.by = "Cell_type")
  plot(CombinePlots(plots = list(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10), nrow = 3))
}




################ plot QC #################
### not using this. refer to plot_grid_platelayout function at the bottom ###

#function for plotting QC
# plot_QC<-function(data,feature_type,treatment){
#   data<-data %>%
#     select(Well_ID,Row,Column,{{feature_type}},{{treatment}}) %>% 
#     mutate(Col=as.character(Column)) %>%
#     mutate(Col=factor(Col,levels=mixedsort(unique(Col)))) %>%
#     mutate(median_value=median(!!rlang::sym(feature_type))) %>%
#     mutate(max_value=max(!!rlang::sym(feature_type))) %>%
#     mutate(min_value=min(!!rlang::sym(feature_type)))
#   p<-ggplot(data, aes(Col,
#                       fct_rev(as_factor(Row)),
#                       fill=!!rlang::sym(feature_type),
#   )
#   ) + 
#     geom_tile() +
#     scale_x_discrete(position = "top") +
#     #scale_fill_continuous(trans="log10") +
#     scale_fill_gradient2(low = "blue",
#                          mid = "white",
#                          high = "red",
#                          midpoint = unique(data$median_value),
#                          name = "RNA read counts") +
#     theme(
#       panel.grid = element_blank(),
#       axis.title.y = element_blank(),
#       axis.title.x = element_blank()
#     ) +
#     geom_text(aes(label=!!rlang::sym(treatment)),size=2)
#   p
# }

#function for plotting QC with changing legend title 
# plot_QC_name<-function(data,feature_type,treatment, attribute_name){
#   data<-data %>%
#     select(Well_ID,Row,Column,{{feature_type}},{{treatment}}) %>% 
#     mutate(Col=as.character(Column)) %>%
#     mutate(Col=factor(Col,levels=mixedsort(unique(Col)))) %>%
#     mutate(median_value=median(!!rlang::sym(feature_type))) %>%
#     mutate(max_value=max(!!rlang::sym(feature_type))) %>%
#     mutate(min_value=min(!!rlang::sym(feature_type)))
#   p<-ggplot(data, aes(Col,
#                       fct_rev(as_factor(Row)),
#                       fill=!!rlang::sym(feature_type),
#   )
#   ) + 
#     geom_tile() +
#     scale_x_discrete(position = "top") +
#     #scale_fill_continuous(trans="log10") +
#     scale_fill_gradient2(low = "blue",
#                          mid = "white",
#                          high = "red",
#                          midpoint = unique(data$median_value),
#                          name = paste0(attribute_name)) +
#     theme(
#       panel.grid = element_blank(),
#       axis.title.y = element_blank(),
#       axis.title.x = element_blank()
#     ) +
#     geom_text(aes(label=!!rlang::sym(treatment)),size=2)
#   p
# }



################ plot RLE #################
RLEplot<-function(count_matrix,ID,feature){
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




################platetools QC plate layout#################
################log scale option#################
plot_grid_platelayout<-function(data, feature_type, legend_label, title_label,n_wells=384, log_scale=FALSE){
  data<-data%>%select(Well_ID, {{feature_type}})
  if (log_scale==TRUE){
    data[{{feature_type}}]<-log10(data[{{feature_type}}])
    min_value<-min(data[{{feature_type}}][,1])
    median_value<-median(data[{{feature_type}}][,1])
    max_value<-max(data[{{feature_type}}][,1])
    p<-raw_map(data=data[{{feature_type}}][,1], well=data$Well_ID, plate=n_wells)+ggtitle(paste0(params$project_code, ": ", title_label)) +
      scale_fill_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = median_value,
        limits = c(floor(min_value), ceiling(max_value)),
        name = paste0(legend_label)
      )+theme(
        plot.title = element_text(size = 15),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)
      )
    p
    
  }else{
    min_value<-min(data[{{feature_type}}][,1])
    median_value<-median(data[{{feature_type}}][,1])
    max_value<-max(data[{{feature_type}}][,1])
    p<-raw_map(data=data[{{feature_type}}][,1], well=data$Well_ID, plate=n_wells)+ggtitle(paste0(params$project_code, ": ", title_label)) +
      scale_fill_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = median_value,
        limits = c(floor(min_value), ceiling(max_value)),
        name = paste0(legend_label)
      )+  theme(
        plot.title = element_text(size = 15),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)
      )
    p
  }
}


