
#' Connect to enrichR
#' This function connect to enrichR database
#' here, we only focus on main enrichR, not fish/fly/yeastenrichR
#' enrichR contains both human and mouse
#' @param dbs name of database in encirhR such as main enrichR or fish/fly/yeastenrichR
#' @return link for enrichR database
#' @export

#' @import dplyr httr2 ggplot2 DT
enrichr_url<-function(dbs="enrichr"){
  # Check if the input matches "enrichr"
  if (!dbs %in% c("enrichr")) {
    stop("Error: The specified database is not supported. Only 'enrichr' is allowed.")
  }
  dbs<-match.arg(dbs)
  switch(dbs,
         "enrichr"="https://maayanlab.cloud/Enrichr/")
}




#' 2.check enrichr connection and get the stats of datasets from enrichR
#' This function returns a data frame with all available data sets from enrichR
#' from enrichR script, the listEnrichrSites function extracts the big table of all available datasets
#' from "datasetStatistics"
#' @param url link of the enrichR database
#' @return data frame of all available data sets from the enrichR link user provided
#' @export
GS_summary<-function(url){

  resp<-request(paste0(url,"datasetStatistics"))|> req_perform()
  if(!resp_status_desc(resp)=="OK"){
    stop(resp_status_desc(resp))
  }
  data<-resp_body_json(resp)$statistics
  stats<-sapply(data, function(x){
    unlist(x)
    as.data.frame(x)
  })
  stats<-as.data.frame(t(stats))
  stats[c(164:215,217:226),]<-stats[c(164:215,217:226), c(2:5,1,6,7)]
  stats<-stats%>%select(geneCoverage, genesPerTerm, libraryName, link, numTerms)
  stats<-sapply(stats, function(x){
    unlist(x)
  })#need to change type
}





#' Download dataset of interest
#' @param url link of enrichR database
#' @param geneset name of data sets from the enrichR database such as KEGG
#' @export
#' @return list of pathways contains sub-lists of genes that annotated in the corresponding pathways
geneset_download<-function(url, geneset){
  resp_lib<-request(paste0(url,"geneSetLibrary?mode=text&libraryName=",geneset))|> req_perform()
  resp_lib_split<-strsplit(resp_body_string(resp_lib), split = "\n")[[1]]
  resp_lib_split_pathway<-lapply(resp_lib_split, function(x){
    x<-strsplit(x, split = "\t")[[1]]
  })
  names(resp_lib_split_pathway)<-unlist(lapply(resp_lib_split_pathway, function(x) x[1]))

  resp_lib_split_pathway<-lapply(resp_lib_split_pathway, function(x){
    g<-x[3:length(x)]#genes start at index 3
    g<-g[g!=""]
    g<-unique(g)
  })
}



#' Perform gene enrichment analysis
#' by doing a hyper-geometric test
#' @param deg list of DEGs
#' @param genesets selected pathway datasets like KEGG from enrichR
#' @param background number of genes that considered as background, options as "human": number of human genes
#' or "geneset" genes in the previously selected pathway dataset
#' @importFrom methods is
#' @importFrom stats p.adjust
#' @export
#' @return a data frame of enrichment analysis results with all the statistics

hyper_enrich_bg<-function(deg, #list of DEGs
                          genesets, #pathway from enrichr
                          background="human" #define the number of genes in your experiment (interger), or "human" to use human genome or
                          #"geneset" to use number genes in the gene set collection
){

  #checking input
  if(!is(deg,"vector")){
    stop("DEGs expected to be a vector of gene symbols\n")
  }
  if(!is(genesets,"list")){
    stop("Genesets expected to be a list of non-null names\n")
  }

  if(!is.numeric(background)){
    #set background/universe gene number
    background<-switch(background,
                       "human"=20000,#as enrichr set
                       "geneset"=length(unique(unlist(genesets))),
    )
    background<-as.numeric(background)

  }else{
    background=background
  }

  deg<-unique(deg)
  genesets<-lapply(genesets, unique)
  deg_genesets<-deg[deg %in% unique(unlist(genesets))]
  n_hits<-sapply(genesets, function(x,y) length(intersect(x,y)), deg_genesets)#q
  n_hits_updatebg<-n_hits!=0#updating background
  genesets<-genesets[n_hits_updatebg]#updating background
  n_hits<-n_hits[n_hits>0]#exlude 0 overlapping terms

  # Check if there are valid genesets left
  if (length(genesets) == 0 || length(n_hits) == 0) {
    stop("No overlapping genes found between DEGs and genesets. Cannot compute enrichment.\n")
  }


  n_genesets<-sapply(genesets,length)#m
  n_minus<-as.numeric(background) - as.numeric(n_genesets) # n
  n_deg<-length(deg)#k

  #hyper-geometric
  pvals<-stats::phyper(q=n_hits-1,
                       m = as.numeric(n_genesets),
                       n = as.numeric(n_minus),
                       k = as.numeric(n_deg),
                       lower.tail = FALSE)



  #put into a dataframe
  res_data<-data.frame(Term=names(genesets),
                       Overlap=paste0(n_hits,"/",n_genesets),
                       P.value=pvals,
                       Adjusted.Pvalue=p.adjust(pvals, "BH"),
                       Genes=sapply(genesets, function(x,y) paste(intersect(x,y),collapse = ";"), deg_genesets))

  res_data<-arrange(res_data, P.value)
  return(results=res_data)

}


#' generate a basic bar plot
#' @param df output data frame from enrichment analysis
#' @importFrom stats reorder
#' @import ggplot2
#' @export
#' @return bar plot showing top 10 enriched pathways
plot_enrich<-function(df, numTerms=10){

  map<-as.numeric(sub("/.*","", df$Overlap))
  gsn<-as.numeric(sub(".*/","", df$Overlap))
  generatio<-map/gsn

  plot_df<-data.frame(Term=df$Term,
                      Hits=map,
                      GeneRatio=generatio,
                      P.value=df$P.value,
                      AdjP.value=df$Adjusted.Pvalue)

  p<-ggplot(plot_df[1:numTerms,], aes(x=Hits,y=reorder(Term, -AdjP.value),
                                      fill=AdjP.value))+geom_bar(stat="identity")+theme_bw()+scale_fill_continuous(
                                        low="red", high="blue")+
    guides(fill=guide_colorbar(title = "Adjusted P.value"))+ylab("Term")

  return(p)

}

utils::globalVariables(c("Hits", "Term", "AdjP.value","P.value","geneCoverage", "genesPerTerm", "libraryName",
                         "link", "numTerms"))
