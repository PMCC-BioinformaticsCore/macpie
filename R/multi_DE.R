
multi_DE <- function(data = NULL, treatment_samples = NULL, control_samples = "DMSO", method = "edgeR" ) {

  control.compute=list(save.memory=TRUE)
  de_list <- pmclapply(treatments, function(x) {
    result <- differential_expression(mac, x, control_samples, method = "edgeR")
    result$combined_id <- x
    return(result)
  }, mc.cores = num_cores)

}

