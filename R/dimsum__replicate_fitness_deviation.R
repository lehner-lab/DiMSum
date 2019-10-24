
#' dimsum__replicate_fitness_deviation
#'
#' Calculate (square-root of sum of squared) differences between replicate fitness values (for all variants) given normalisation parameters
#'
#' @param p vector of normalisation parameters (required)
#' @param fitness_mat matrix of fitness values (required)
#' @param all_reps list of replicates to retain (required)
#'
#' @return the square-root of sum of squared differences to mean over replicates
#' @export
dimsum__replicate_fitness_deviation <- function(
  p,
  fitness_mat,
  all_reps
  ){
  #Number of input and output replicates
  all_reps_str <- paste0(all_reps, collapse="")
  #Normalised fitness (scale and center/shift)
  F_norm <- (fitness_mat + 
    matrix(p[(nchar(all_reps_str) + 1):(2*nchar(all_reps_str))], 
      nrow = nrow(fitness_mat), ncol = nchar(all_reps_str), byrow = T)) * 
    matrix(p[1:nchar(all_reps_str)], 
      nrow = nrow(fitness_mat), ncol = nchar(all_reps_str), byrow = T)
  #Average fitness (center/shift only)
  F_avg <- rowMeans(fitness_mat + matrix(p[(nchar(all_reps_str) + 1):(2*nchar(all_reps_str))], 
    nrow = nrow(fitness_mat), ncol = nchar(all_reps_str), byrow = T))
  diffF <- sqrt(rowSums((F_norm - F_avg)^2))
  return(sum(diffF))
}
