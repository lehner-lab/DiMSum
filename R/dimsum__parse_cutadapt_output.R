
#' dimsum__parse_cutadapt_output
#'
#' Parse cutadapt stdout.
#'
#' @param file_path path to cutadapt output file (required)
#'
#' @return a list of read trimming count statistics
#' @export
dimsum__parse_cutadapt_output <- function(
  file_path
  ){
  temp_out <- system(paste0("cat ", file_path), intern=TRUE)
  name_read1 <- basename(rev(unlist(strsplit(temp_out[2], ' ')))[2])
  name_read2 <- basename(rev(unlist(strsplit(temp_out[2], ' ')))[1])
  temp_out <- temp_out[c(9:13, grep("=== First|=== Second|Sequence: ", temp_out))]
  #Total reads processed
  total_reads <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[1], ' ')))[1]))
  #Total reads trimmed
  total_read1_trimmed <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[2], ' ')))[2]))
  total_read2_trimmed <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[3], ' ')))[2]))
  #Read number
  read_number <- c("read2", "read1")[as.numeric(grepl("First", temp_out[grep("=== ", temp_out)]))+1]
  #Adapter position
  adapter_position <- c("a3", "a5")[as.numeric(grepl(" 5';", temp_out[grep("Sequence: ", temp_out)]))+1]
  #Adapter type
  adapter_type <- c("unlinked", "linked")[as.numeric(grepl("Type: linked", temp_out[grep("Sequence: ", temp_out)]))+1]
  #Read number
  trim_counts <- as.numeric(sapply(lapply(strsplit(temp_out[grep("Sequence: ", temp_out)], ' '), rev), '[', 2))
  #Final totals
  total_read1_a3 <- ifelse(sum(read_number=="read1" & adapter_position=="a3"),trim_counts[read_number=="read1" & adapter_position=="a3"],0)
  total_read1_a5 <- ifelse(sum(read_number=="read1" & adapter_position=="a5"),trim_counts[read_number=="read1" & adapter_position=="a5"],0)
  total_read1_both <- total_read1_a3+total_read1_a5-total_read1_trimmed
  total_read2_a3 <- ifelse(sum(read_number=="read2" & adapter_position=="a3"),trim_counts[read_number=="read2" & adapter_position=="a3"],0)
  total_read2_a5 <- ifelse(sum(read_number=="read2" & adapter_position=="a5"),trim_counts[read_number=="read2" & adapter_position=="a5"],0)
  total_read2_both <- total_read2_a3+total_read2_a5-total_read2_trimmed
  #If all adapters linked (mixture of linked and unlinked adapters not supported)
  if(sum(adapter_type!="linked")==0){
    total_read1_a5 <- total_read1_a3
    total_read1_both <- total_read1_a3
    total_read2_a5 <- total_read2_a3
    total_read2_both <- total_read2_a3
  }
  return(list(
    name_read1 = name_read1,
    name_read2 = name_read2,
    total_reads = total_reads,
    total_read1_a3 = total_read1_a3,
    total_read1_a5 = total_read1_a5,
    total_read1_both = total_read1_both,
    total_read2_a3 = total_read2_a3, 
    total_read2_a5 = total_read2_a5,
    total_read2_both = total_read2_both))
}