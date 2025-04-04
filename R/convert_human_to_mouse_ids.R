#' For a given gene set of human symbols return mouse symbols
#' 
#' @param gene_list List of human genes to be converted to mouse IDs
#' @returns genes
#' @export
#'
#' @examples
#' \dontrun{
#' convert_human_to_mouse(c("BRCA1", "TRAF1", "MYBL1"))
#' }

convert_human_to_mouse <- function(gene_list) {
    
  validate_inputs <- function(data, genesets, species, direction, p_value_cutoff, n_distinct) {
    if (!inherits(gene_list, "character")) {
      stop("Error: argument 'gene_list' must be of type character.")
    }
  }
  
  mouse_human_genes_file <- "inst/extdata/annotation/mouse_human_genes.csv"
  if(!file.exists(mouse_human_genes_file)){
    mouse_human_genes = read.csv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
    write.csv(mouse_human_genes,file="inst/extdata/annotation/mouse_human_genes.csv",row.names=F)
  } else {
    mouse_human_genes<-read.csv(mouse_human_genes_file)
  }

  output = c()
  for (gene in gene_list) {
    class_key <- mouse_human_genes %>%
      filter(Symbol == gene & Common.Organism.Name == "human") %>%
      pull(DB.Class.Key)
    
    if (length(class_key) > 0) {
      mouse_genes <- mouse_human_genes %>%
        filter(DB.Class.Key == class_key & Common.Organism.Name == "mouse, laboratory") %>%
        pull(Symbol)
      
      # Append while preserving input gene order
      output <- c(output, mouse_genes)
    } else {
      # If no match found, optionally add NA or skip
      output <- c(output, NA)  # or don't include this line if you want to skip
    }
  }
  return (output)
}
