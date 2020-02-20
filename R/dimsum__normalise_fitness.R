
#' dimsum__normalise_fitness
#'
#' Normalise fitness and error for differences in number of generations.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param input_dt input data.table (required)
#' @param all_reps list of replicates to retain (required)
#' @param fitness_suffix fitness and sigma suffix (default:"")
#'
#' @return data.table with normalised fitness and error
#' @export
#' @import data.table
dimsum__normalise_fitness <- function(
  dimsum_meta,
  input_dt,
  all_reps,
  fitness_suffix=""
  ){

  #Get generations data
  generations_unique <- unique(dimsum_meta[["exp_design"]][dimsum_meta[["exp_design"]][,"selection_id"]==1,c("sample_name", "experiment", "generations")])
  #Mean generations per biological replicate
  generations_mean <- tapply(generations_unique[,"generations"], paste0("e", generations_unique[,"experiment"]), mean)

  #Normalise fitness and sigma by mean generations
  for (E in all_reps) {
    for(i in c("fitness", "sigma")){
      input_dt[,paste0(i, E, fitness_suffix) := log2(exp(.SD/generations_mean[paste0("e", E)])),,.SDcols = paste0(i, E, fitness_suffix)]
    }
  }

  return(input_dt)

}
