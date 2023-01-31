
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
  #Normalised fitness (scale and center/shift)
  F_norm <- (fitness_mat + 
    matrix(p[(length(all_reps) + 1):(2*length(all_reps))], 
      nrow = nrow(fitness_mat), ncol = length(all_reps), byrow = T)) * 
    matrix(p[1:length(all_reps)], 
      nrow = nrow(fitness_mat), ncol = length(all_reps), byrow = T)
  #Average fitness (center/shift only)
  F_avg <- rowMeans(fitness_mat + matrix(p[(length(all_reps) + 1):(2*length(all_reps))], 
    nrow = nrow(fitness_mat), ncol = length(all_reps), byrow = T))
  diffF <- sqrt(rowSums((F_norm - F_avg)^2))
  return(sum(diffF))
}
