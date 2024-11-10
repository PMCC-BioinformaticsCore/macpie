
suppressMessages(library(dplyr))
suppressMessages(library(httr2))
suppressMessages(library(ggplot2))
suppressMessages(library(DT))

### writing functions ###
#get database from enrichR

#1. connect to enrichR
#here, we only focus on main enrichR, not fish/fly/yeastenrichR
#enrichR contains both human and mouse
enrichr_url<-function(dbs="enrichr"){
  dbs<-match.arg(dbs)
  switch(dbs,
         "enrichr"="https://maayanlab.cloud/Enrichr/")
}




#2.check enrichr connection and get the stats of datasets from enrichR
#from enrichR script, the listEnrichrSites function extracts the big table of all available datasets
#from "datasetStatistics"
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





#3. download those datasets
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



# 4. Gene enrichment analysis
# make it as a function for hyper-geometric test
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

  if(typeof(background)!="integer"){
    #set background/universe gene number
    background<-switch(background,
                       "human"=as.character(20000),#as enrichr set
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


  n_genesets<-sapply(genesets,length)#m
  n_minus<-background-n_genesets#n
  n_deg<-length(deg)#k

  #hyper-geometric
  pvals<-stats::phyper(q=n_hits-1,
                       m=n_genesets,
                       n=n_minus,
                       k=n_deg,
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

#5. basic bar plot
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
