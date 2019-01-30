
#' dimsum__parse_cutadapt_output_single_end
#'
#' Parse cutadapt stdout (single-end mode).
#'
#' @param file_path path to cutadapt output file (required)
#'
#' @return a list of read trimming count statistics
#' @export
dimsum__parse_cutadapt_output_single_end <- function(
  file_path
  ){
  temp_out <- system(paste0("cat ", file_path), intern=TRUE)
  name_read1 <- basename(rev(unlist(strsplit(temp_out[2], ' ')))[1])
  temp_out <- temp_out[c(9:12, grep("=== Adapter|Sequence: ", temp_out))]
  #Total reads processed
  total_reads <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[1], ' ')))[1]))
  #Total reads trimmed
  total_read1_trimmed <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[2], ' ')))[2]))
  #Adapter position
  adapter_position <- c("a3", "a5")[as.numeric(grepl(" 5';", temp_out[grep("Sequence: ", temp_out)]))+1]
  #Adapter type
  adapter_type <- c("unlinked", "linked")[as.numeric(grepl("Type: linked", temp_out[grep("Sequence: ", temp_out)]))+1]
  #Read number
  trim_counts <- as.numeric(sapply(lapply(strsplit(temp_out[grep("Sequence: ", temp_out)], ' '), rev), '[', 2))
  #Final totals
  total_read1_a3 <- ifelse(sum(adapter_position=="a3"),trim_counts[adapter_position=="a3"],0)
  total_read1_a5 <- ifelse(sum(adapter_position=="a5"),trim_counts[adapter_position=="a5"],0)
  total_read1_both <- total_read1_a3+total_read1_a5-total_read1_trimmed
  #If all adapters linked (mixture of linked and unlinked adapters not supported)
  if(sum(adapter_type!="linked")==0){
    total_read1_a5 <- total_read1_a3
    total_read1_both <- total_read1_a3
  }
  return(list(
    name_read1 = name_read1,
    total_reads = total_reads,
    total_read1_a3 = total_read1_a3,
    total_read1_a5 = total_read1_a5,
    total_read1_both = total_read1_both))
}
