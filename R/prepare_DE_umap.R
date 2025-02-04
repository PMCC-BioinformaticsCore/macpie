
prepare_DE_umap <- function(de_list = NULL){

  df<-do.call("rbind",de_list)

  df_wide <- df %>%
    select(gene, combined_id, metric) %>%
    pivot_wider(names_from = combined_id, values_from = metric)

  set.seed(1)
  df_umap<-umap(t(df_wide[,-1]))

  #plot controls on UMAP
  df_umap_data<-df_umap$layout %>%
    as.data.frame() %>%
    rownames_to_column("combined_id") %>%
    left_join(.,mac@meta.data,join_by(combined_id==combined_id)) %>%
    rename("UMAP_1"="V1", "UMAP_2"="V2")
}


