
#create_dimsum_dir
#
# Save experiment metadata.
#
# dimsum_meta: an experiment metadata object (required)
#
# Returns: nothing.
#
save_metadata <- function(dimsum_meta){
  save(dimsum_meta, file=file.path(dimsum_meta[["project_path"]], paste0(dimsum_meta[["project_name"]], '_metadata.RData')))
}
