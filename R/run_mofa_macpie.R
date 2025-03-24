#' Run wrapper function for MOFA2
#'
#' @description
#' Wrapper function to train an untrained MOFA object from the package MOFA2. To ensure successful build of the conda environment in the docker container, this wrapper sets use_basilisk = TRUE. Please note that not using the basilisk package will make the function fail.
#'
#' @param object an untrained MOFA object
#' @param outfile output file for the model (.hdf5 format). If NULL, a temporary file is created.
#' @param save_data logical indicating whether to save the training data in the hdf5 file. This is useful for some downstream analysis (mainly functions with the prefix plot_data), but it can take a lot of disk space.
#' @importFrom MOFA2 run_mofa
#' @returns A trained MOFA2 object
#' @export
#'
#' @examples
#' file <- system.file("extdata", "test_data.RData", package = "MOFA2")
#' load(file)
#' MOFAmodel <- create_mofa(dt)
#' MOFAmodel <- prepare_mofa(MOFAmodel)
#' run_mofa_macpie(object = MOFAmodel, outfile = NULL, save_data = TRUE)

run_mofa_macpie <- function(object, outfile = NULL, save_data = TRUE){
  # Helper function to validate input data
  validate_inputs <- function(object, outfile, save_data) {
    if (!inherits(MOFAobject, "MOFA")) {
      stop("Error: argument 'object' must be a MOFA object")
    }
    if (!inherits(outfile, "character") & !inherits(outfile, "NULL")) {
      stop("Error: argument 'outfile' must be a character")
      if (!grepl('\\.hdf5', outfile)){
        stop("Invalid file extention. Argument 'outfile' must end in '.hdf5'")
      }
    }
    if (!inherits(save_data, "logical")) {
      stop("Error: argument 'save_data' must be a logical (TRUE/FALSE)")
    }
  }
  # Validate inputs
  validate_inputs <- validate_inputs(object, outfile, save_data)

  #run original 'run_mofa' function with basilisk for conda
  run_mofa(object = object,
           outfile = outfile,
           save_data = save_data,
           use_basilisk = TRUE)
}
