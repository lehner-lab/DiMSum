
#' dimsum__infer_growth_rates
#'
#' Infer growth rates from fitness.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param input_dt input data.table (required)
#' @param all_reps list of replicates to retain (required)
#'
#' @return data.table with growth rate and growth rate error
#' @export
#' @import data.table
dimsum__infer_growth_rates <- function(
  dimsum_meta,
  input_dt,
  all_reps
  ){

  #If input data.table empty, return empty data.table
  if(nrow(input_dt)==0){
    return(data.table())
  }

  #Get selection_time data
  selection_time_unique <- unique(dimsum_meta[["exp_design"]][dimsum_meta[["exp_design"]][,"selection_id"]==1,c("sample_name", "experiment", "selection_time")])
  #Mean selection_time per biological replicate
  selection_time_mean <- tapply(selection_time_unique[,"selection_time"], paste0("e", selection_time_unique[,"experiment"]), mean)

  #Get cell_density data for output
  cell_density_output_unique <- unique(dimsum_meta[["exp_design"]][dimsum_meta[["exp_design"]][,"selection_id"]==1,c("sample_name", "experiment", "cell_density")])
  #Mean cell_density per output biological replicate
  cell_density_output_mean <- tapply(cell_density_output_unique[,"cell_density"], paste0("e", cell_density_output_unique[,"experiment"]), mean)

  #Get cell_density data for input
  cell_density_input_unique <- unique(dimsum_meta[["exp_design"]][dimsum_meta[["exp_design"]][,"selection_id"]==0,c("sample_name", "experiment", "cell_density")])
  #Mean cell_density per input biological replicate
  cell_density_input_mean <- tapply(cell_density_input_unique[,"cell_density"], paste0("e", cell_density_input_unique[,"experiment"]), mean)

  #Infer growth rates from fitness
  for (E in all_reps) {
    input_dt[error_model==T, paste0("growthrate", E) := log((.SD[[2]]/sum(.SD[[2]], na.rm = T)*cell_density_output_mean[paste0("e", E)]) / (.SD[[1]]/sum(.SD[[1]], na.rm = T)*cell_density_input_mean[paste0("e", E)])) / selection_time_mean[paste0("e", E)],,.SDcols = paste0("count_e", E, "_s", c(0, 1))]

    gr_cor <- input_dt[!is.na(get(paste0("growthrate", E))) & !is.infinite(get(paste0("growthrate", E))) & error_model==T,cor(.SD[[1]], .SD[[2]], use = "pairwise.complete"),,.SDcols = c(paste0("growthrate", E), paste0("fitness", E, "_uncorr"))]
    dimsum__status_message(paste0("Fitness vs. growth rate correlation for replicate ", E, ": ", round(gr_cor, 2), "\n"))

    gr_lm <- lm(y~x, input_dt[!is.na(get(paste0("growthrate", E))) & !is.infinite(get(paste0("growthrate", E))) & error_model==T,.(y = .SD[[1]], x = .SD[[2]]),,.SDcols = c(paste0("growthrate", E), paste0("fitness", E, "_uncorr"))])$coefficients
    input_dt[, paste0("growthrate", E) := .SD[[1]]*gr_lm[2]+gr_lm[1],,.SDcols = paste0("fitness", E, "_uncorr")]
    input_dt[, paste0("growthrate", E, "_sigma") := .SD[[1]]*gr_lm[2],,.SDcols = paste0("sigma", E, "_uncorr")]
  }

  #Merge replicate growth rates and errors
  all_reps_str <- paste0(all_reps, collapse="|")
  fitness_rx <- input_dt[,.SD,.SDcols = grep(paste0("growthrate(", all_reps_str, ")$"),colnames(input_dt))]
  sigma_rx <- input_dt[,.SD,.SDcols = grep(paste0("growthrate(", all_reps_str, ")_sigma"),colnames(input_dt))]
  input_dt[,growthrate := rowSums(fitness_rx/(sigma_rx^2),na.rm=T)/rowSums(1/(sigma_rx^2),na.rm=T)]
  gr_lm <- lm(growthrate~fitness, input_dt[error_model==T])$coefficients
  input_dt[, growthrate := .SD[[1]]*gr_lm[2]+gr_lm[1],,.SDcols = "fitness"]
  input_dt[, growthrate_sigma := .SD[[1]]*gr_lm[2],,.SDcols = "sigma"]

  return(input_dt)

}
