
#' dimsum__parse_cutadapt_output
#'
#' Parse cutadapt stdout.
#'
#' @param file_path path to cutadapt output file (required)
#' @param ran_cutadapt whether this is a bona fide cutadapt output file (default:T)
#' @param ran_cutadapt_cutonly whether cutadapt was run using "Cut" options only (default:F)
#'
#' @return a list of read trimming count statistics
#' @export
dimsum__parse_cutadapt_output <- function(
  file_path,
  ran_cutadapt=T,
  ran_cutadapt_cutonly=F
  ){
  temp_out <- readLines(file_path)
  if(!ran_cutadapt){
    name_read1 <- basename(rev(unlist(strsplit(temp_out[1], ' ')))[1])
    name_read2 <- basename(rev(unlist(strsplit(temp_out[2], ' ')))[1])
    total_reads <- as.integer(rev(unlist(strsplit(temp_out[1], ' ')))[2])/4
    total_read1_a3 <- 0
    total_read1_a5 <- 0
    total_read1_both <- 0
    total_read2_a3 <- 0
    total_read2_a5 <- 0
    total_read2_both <- 0
  }else if(ran_cutadapt_cutonly){
    name_read1 <- basename(rev(unlist(strsplit(temp_out[2], ' ')))[2])
    name_read2 <- basename(rev(unlist(strsplit(temp_out[2], ' ')))[1])
    summary_index <- grep("=== Summary", temp_out)
    temp_out <- temp_out[c((summary_index+2):(summary_index+6), grep("=== First|=== Second|Sequence: ", temp_out))]
    #Total reads processed
    total_reads <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[1], ' ')))[1]))
    total_read1_a3 <- 0
    total_read1_a5 <- 0
    total_read1_both <- 0
    total_read2_a3 <- 0
    total_read2_a5 <- 0
    total_read2_both <- 0  
  }else{
    name_read1 <- basename(rev(unlist(strsplit(temp_out[2], ' ')))[2])
    name_read2 <- basename(rev(unlist(strsplit(temp_out[2], ' ')))[1])
    summary_index <- grep("=== Summary", temp_out)
    temp_out <- temp_out[c((summary_index+2):(summary_index+6), grep("=== First|=== Second|Sequence: ", temp_out))]
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
    #Read1 adapters linked
    if(unique(adapter_type[read_number=="read1"])=="linked"){
      total_read1_a5 <- total_read1_a3
      total_read1_both <- total_read1_a3
    }
    #Read2 adapters linked
    if(unique(adapter_type[read_number=="read2"])=="linked"){
      total_read2_a5 <- total_read2_a3
      total_read2_both <- total_read2_a3
    }
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
