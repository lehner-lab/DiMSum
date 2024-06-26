
#' dimsum__parse_cutadapt_output_single_end
#'
#' Parse cutadapt stdout (single-end mode).
#'
#' @param file_path path to cutadapt output file (required)
#' @param ran_cutadapt whether this is a bona fide cutadapt output file (default:T)
#' @param ran_cutadapt_cutonly whether cutadapt was run using "Cut" options only (default:F)
#'
#' @return a list of read trimming count statistics
#' @export
dimsum__parse_cutadapt_output_single_end <- function(
  file_path,
  ran_cutadapt=T,
  ran_cutadapt_cutonly=F
  ){
  temp_out <- readLines(file_path)
  if(!ran_cutadapt){
    name_read1 <- basename(rev(unlist(strsplit(temp_out[1], ' ')))[1])
    total_reads <- as.integer(rev(unlist(strsplit(temp_out[1], ' ')))[2])/4
    total_read1_a3 <- 0
    total_read1_a5 <- 0
    total_read1_both <- 0
  }else if(ran_cutadapt_cutonly){
    name_read1 <- basename(rev(unlist(strsplit(temp_out[2], ' ')))[1])
    summary_index <- grep("=== Summary", temp_out)
    temp_out <- temp_out[c((summary_index+2):(summary_index+5), grep("=== Adapter|Sequence: ", temp_out))]
    #Total reads processed
    total_reads <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[1], ' ')))[1]))
    total_read1_a3 <- 0
    total_read1_a5 <- 0
    total_read1_both <- 0
  }else{
    name_read1 <- basename(rev(unlist(strsplit(temp_out[2], ' ')))[1])
    summary_index <- grep("=== Summary", temp_out)
    temp_out <- temp_out[c((summary_index+2):(summary_index+5), grep("=== Adapter|Sequence: ", temp_out))]
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
    trim_counts_linked_a5 <- as.numeric(sapply(lapply(strsplit(temp_out[grep("Sequence: .*linked", temp_out)], ' '), rev), '[', 6))
    #Final totals
    total_read1_a3 <- ifelse(sum(adapter_position=="a3"),trim_counts[adapter_position=="a3"],0)
    total_read1_a5 <- ifelse(sum(adapter_position=="a5"),trim_counts[adapter_position=="a5"],0)
    total_read1_a5 <- ifelse(sum(adapter_type=="linked"),trim_counts_linked_a5[adapter_type=="linked"],total_read1_a5)
    total_read1_both <- total_read1_a3+total_read1_a5-total_read1_trimmed
  }
  return(list(
    name_read1 = name_read1,
    total_reads = total_reads,
    total_read1_a3 = total_read1_a3,
    total_read1_a5 = total_read1_a5,
    total_read1_both = total_read1_both))
}
