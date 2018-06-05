#!/software/bl/el7.2/R-3.3.0/bin/Rscript
.libPaths(c('/software/bl/el7.2/R/R-3.3.0_packages/', '/software/bl/el7.2/R/R-3.3.0/library/'))

#Software version
version_info <- "DiMSum_v0.1"
#Display software version
print(version_info)

#TODO:
#check for required binaries and scripts before starting pipeline (exit gracefully)
#Gzipped files assumed to have .gz extension

###########################
### PACKAGES
###########################

require(data.table)
library(seqinr)
library(optparse)
library(parallel)
library(reshape2)
library(ggplot2)
library(plyr)

###########################
### COMMAND-LINE OPTIONS
###########################

option_list <- list(
  make_option(opt_str=c("--fastqFileDir", "-i"), help = "Path to directory with input FASTQ files"),
  make_option(opt_str=c("--fastqFileExtension", "-l"), default='.fastq', help = "FASTQ file extension"),
  make_option(opt_str=c("--gzipped", "-g"), type="logical", default=T, help = "Are FASTQ files are gzipped?"),
  make_option(opt_str=c("--stranded"), type="logical", default=T, help = "Is the library design stranded?"),
  make_option(opt_str=c("--barcodeDesignPath", "-b"), help = "Path to barcode design file (tab-separated plain text file with barcode design)"),
  make_option(opt_str=c("--barcodeErrorRate"), type="double", default=0.25, help = "Maximum allowed error rate for the barcode (default:0.25)"),
  make_option(opt_str=c("--experimentDesignPath", "-e"), help = "Path to experimental design file (tab-separated plain text file with replicate structure)"),
  make_option(opt_str=c("--cutadaptCut5First"), type="integer", help = "cutadapt: remove bases from start of first read (before adapter trimming)"),
  make_option(opt_str=c("--cutadaptCut5Second"), type="integer", help = "cutadapt: remove bases from start of second read (before adapter trimming)"),
  make_option(opt_str=c("--cutadaptCut3First"), type="integer", help = "cutadapt: remove bases from end of first read (before adapter trimming)"),
  make_option(opt_str=c("--cutadaptCut3Second"), type="integer", help = "cutadapt: remove bases from end of second read (before adapter trimming)"),
  make_option(opt_str=c("--cutadapt5First"), help = "cutadapt: sequence of an adapter ligated to the 5' end (of the first read)"),
  make_option(opt_str=c("--cutadapt5Second"), help = "cutadapt: sequence of an adapter ligated to the 5' end (of the second read)"),
  make_option(opt_str=c("--cutadapt3First"), help = "cutadapt: sequence of an adapter ligated to the 3' end (of the first read)"),
  make_option(opt_str=c("--cutadapt3Second"), help = "cutadapt: sequence of an adapter ligated to the 3' end (of the second read)"),
  make_option(opt_str=c("--cutadaptMinLength", "-n"), type="integer", default=50, help = "cutadapt: Discard reads shorter than LENGTH (default:50)"),
  make_option(opt_str=c("--cutadaptErrorRate", "-a"), type="double", default=0.2, help = "cutadapt: Maximum allowed error rate (default:0.2)"),
  make_option(opt_str=c("--cutadaptDiscardUntrimmed"), type="logical", default=F, help = "cutadapt: Discard untrimmed read pairs (default:F)"),
  make_option(opt_str=c("--usearchMinQual", "-q"), type="integer", help = "USEARCH: minimum observed base quality to retain read pair"),
  make_option(opt_str=c("--usearchMaxee", "-m"), type="double", help = "USEARCH: maximum number of expected errors to retain read pair"),
  make_option(opt_str=c("--usearchMinovlen"), type="integer", default=16, help = "USEARCH: discard pair if alignment is shorter than given value (default:16)"),
  make_option(opt_str=c("--outputPath", "-o"), help = "Path to directory to use for output files"),
  make_option(opt_str=c("--projectName", "-p"), help = "Project name"),
  make_option(opt_str=c("--wildtypeSequence", "-w"), help = "Wild-type nucleotide sequence"),
  make_option(opt_str=c("--maxAAMutations", "-x"), type="integer", default=2, help = "Maximum considered amino acid sequence mutations"),
  make_option(opt_str=c("--startStage", "-s"), type="integer", default=1, help = "Start at a specified pipeline stage"),
  make_option(opt_str=c("--stopStage", "-t"), type="integer", default=0, help = "Stop at specified pipeline stage (default:0, no stop condition)"),
  make_option(opt_str=c("--numCores", "-c"), type="integer", default=1, help = "Number of available CPU cores")  
)

arg_list <- parse_args(OptionParser(option_list=option_list))
#Display arguments
print(arg_list)

###########################
### FUNCTIONS
###########################

#sum_datatable_columns
#
# Replace a subset of columns with a single column containing the row sums.
#
# dt: Input data.table (required)
# column_patterns: character vector of column patterns to match (required)
# suffix: a character suffix for the appended column (optional)
#
# Returns: a data.table where a subset of columns is replaced with a single column containing the row sums.
#
sum_datatable_columns <- function(dt, column_patterns, suffix=""){
  for(this_column_pattern in column_patterns){
    #Sum columns with given pattern
    temp_data <- apply(dt[,grep(this_column_pattern, colnames(dt)), with=FALSE], 1, sum)
    #Remove summed columns
    dt <- dt[,-grep(this_column_pattern, colnames(dt)), with=FALSE]
    #Append result column 
    dt[,paste0(this_column_pattern, suffix)] <- temp_data
  }
  #Return data.table
  return(dt)
}

#cbind.fill
#
# cbind a list of data.frames of same row names but unequal number.
#
# df_list: a list of data.frames (required)
#
# Returns: a single data.frame where empty rows are filled with NAs.
#
cbind.fill <- function(df_list){
    nm <- lapply(df_list, as.matrix)
    n <- max(sapply(nm, nrow)) 
    temp_rownames <- unique(as.character(unlist(sapply(nm, rownames))))
    temp_df <- as.data.frame(do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(, n-nrow(x), ncol(x))))))
    rownames(temp_df) <- temp_rownames
    temp_df
  }

#dimsum_stage_demultiplex
#
# Run demultiplexing on all fastq files using cutadapt.
#
# dimsum_meta: an experiment metadata object (required)
# demultiplex_outpath: demultiplex output path (required)
# execute: whether or not to execute the system command (default: TRUE)
#
# Returns: an updated experiment metadata object.
#
dimsum_stage_demultiplex <- function(
  dimsum_meta,
  demultiplex_outpath,
  execute = TRUE
  ){
  #Create demultiplex directory (if doesn't already exist)
  demultiplex_outpath <- gsub("/$", "", demultiplex_outpath)
  dir.create(demultiplex_outpath)
  #Demultiplex parameters specified?
  if( !'barcode_design' %in% names(dimsum_meta) ){
    message("Skipping demultiplexing (assuming all fastq files already demultiplexed)")
    return(dimsum_meta)
  }else{
    fastq_pair_list <- unique(dimsum_meta[['barcode_design']][,c('pair1', 'pair2')])
    rownames(fastq_pair_list) = 1:dim(fastq_pair_list)[1]
    #Trim FASTQ file pairs
    message("Demultiplexing FASTQ files with cutadapt:")
    all_fastq <- file.path(dimsum_meta[["exp_design"]]$pair_directory, c(dimsum_meta[['barcode_design']]$pair1, dimsum_meta[['barcode_design']]$pair2))
    print(unique(all_fastq))
    message("Processing...")
    for(pair_name in rownames(fastq_pair_list)){
      #TODO: saber binary path specifiable on commandline?
      print(fastq_pair_list[pair_name,])
      #Check if this system command should be executed
      if(execute){
        #
        temp_design <- dimsum_meta[['barcode_design']][grep(fastq_pair_list[pair_name,'pair1'], dimsum_meta[['barcode_design']][,'pair1']),]
        write(
          x = c(rbind(paste0('>', temp_design$new_pair_prefix), paste0('^', temp_design$barcode))), 
          file = file.path(demultiplex_outpath, paste0('demultiplex_barcode-file_', pair_name, '.fasta')), 
          sep="\n")
        #Demultiplex using cutadapt
        temp_out = system(paste0(
          "cutadapt",
          " -g file:",
          file.path(demultiplex_outpath, paste0('demultiplex_barcode-file_', pair_name, '.fasta')),
          " -G file:",
          file.path(demultiplex_outpath, paste0('demultiplex_barcode-file_', pair_name, '.fasta')),
          " -e ",
          as.character(dimsum_meta[["barcodeErrorRate"]]),
          " --no-indels ",
          " --untrimmed-output ",
          file.path(demultiplex_outpath, paste0(fastq_pair_list[pair_name,][1], ".demultiplex.unknown.fastq")),
          " --untrimmed-paired-output ",
          file.path(demultiplex_outpath, paste0(fastq_pair_list[pair_name,][2], ".demultiplex.unknown.fastq")),
          " -o ",
          file.path(demultiplex_outpath, "{name}1.fastq"),
          " -p ",
          file.path(demultiplex_outpath, "{name}2.fastq"),
          " ",
          file.path(dimsum_meta[["exp_design"]]$pair_directory, fastq_pair_list[pair_name,][1]),
          " ",
          file.path(dimsum_meta[["exp_design"]]$pair_directory, fastq_pair_list[pair_name,][2]),
          " > ",
          file.path(demultiplex_outpath, paste0(fastq_pair_list[pair_name,][1], ".demultiplex.stdout")),
          " 2> ",
          file.path(demultiplex_outpath, paste0(fastq_pair_list[pair_name,][1], ".demultiplex.stderr"))))
      }
    }
    #New experiment metadata
    dimsum_meta_new <- dimsum_meta
    #Update fastq metadata
    dimsum_meta_new[['exp_design']]$pair_directory <- demultiplex_outpath
    return(dimsum_meta_new)
  }
}

#dimsum_stage_fastqc_report
#
# Generate FASTQC summary plots for all fastq files.
#
# dimsum_meta: an experiment metadata object (required)
# report_outpath: FASTQC report output path (required)
#
# Returns: an updated experiment metadata object
#
dimsum_stage_fastqc_report <- function(
  dimsum_meta,
  report_outpath
  ){
  #Create report directory (if doesn't already exist)
  report_outpath <- gsub("/$", "", report_outpath)
  dir.create(report_outpath)
  #Get results for all fastq files
  for(col_name in c('pair1_fastqc', 'pair2_fastqc')){
    fastqc_files <- file.path(dimsum_meta[['exp_design']]$fastqc_directory, dimsum_meta[['exp_design']][,col_name])
    fastqc_list <- list()
    encoding <- ''
    for(f in fastqc_files){
      temp_out <- system(paste0("head -n ", 500, ' ', f), intern=TRUE)
      filename <- gsub('Filename\\t', '', temp_out[4])
      encoding <- gsub('Encoding\\t', '', temp_out[6])
      temp_nlines <- grep('>>END_MODULE', temp_out)[2]
      temp_out <- temp_out[c(13:(temp_nlines-1))]
      temp_out_data <- strsplit(temp_out[2:length(temp_out)], '\\t')
      temp_df <- as.data.frame(do.call('rbind', lapply(lapply(temp_out_data, '[', -1), as.numeric)))
      rownames(temp_df) <- sapply(temp_out_data, '[', 1)
      colnames(temp_df) <- unlist(strsplit(temp_out[1], '\\t'))[-1]
      fastqc_list[[filename]] <- temp_df
    }
    fastqc_df1 <- cbind.fill(lapply(fastqc_list, '[', 'Mean'))
    colnames(fastqc_df1) <- names(fastqc_list)
    fastqc_df1$base_position <- 1:length(rownames(fastqc_df1))
    fastqc_df2 <- cbind.fill(lapply(fastqc_list, '[', '10th Percentile'))
    colnames(fastqc_df2) <- names(fastqc_list)
    fastqc_df2$base_position <- 1:length(rownames(fastqc_df1))
    #Plot
    plot_df1 <- melt(fastqc_df1, id="base_position")
    plot_df1$statistic <- 'Mean'
    plot_df2 <- melt(fastqc_df2, id="base_position")
    plot_df2$statistic <- '10th Percentile'
    plot_df <- rbind(plot_df1, plot_df2)
    #Remove NAs
    plot_df <- plot_df[!is.na(plot_df$value),]
    d <- ggplot(plot_df, aes(base_position, value, color = variable)) +
      geom_line() +
      geom_hline(yintercept=c(20, 28), linetype = 2) +
      theme_bw() +
      coord_cartesian(ylim = c(0, max(plot_df$value))) +
      scale_x_continuous(
      breaks = (1:length(rownames(fastqc_df1)))[seq(1, length(rownames(fastqc_df1)), 5)],
      label = rownames(fastqc_df1)[seq(1, length(rownames(fastqc_df1)), 5)]) +
      labs(x = "Position in read (bp)", y = "Quality score", title = paste0("Quality scores across all bases (", encoding, ")"))
    d <- d + facet_wrap(~statistic, nrow=2, ncol=1)
    ggsave(file.path(report_outpath, paste0('dimsum_stage_fastqc_report_', col_name, '.png')), d, width=12, height=8)
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Update fastq metadata
  return(dimsum_meta_new)
}

#dimsum_stage_fastqc
#
# Run FASTQC on all fastq files.
#
# dimsum_meta: an experiment metadata object (required)
# fastqc_outpath: FASTQC output path (required)
# execute: whether or not to execute the system command (default: TRUE)
#
# Returns: an updated experiment metadata object.
#
dimsum_stage_fastqc <- function(
  dimsum_meta,
  fastqc_outpath,
  execute = TRUE,
  report = TRUE,
  report_outpath = NULL
  ){
  #Create FASTQC directory (if doesn't already exist)
  fastqc_outpath <- gsub("/$", "", fastqc_outpath)
  dir.create(fastqc_outpath)
  #Run FASTQC on all fastq files
  message("Running FASTQC on all files:")
  all_fastq <- file.path(dimsum_meta[['exp_design']]$pair_directory, c(dimsum_meta[['exp_design']]$pair1, dimsum_meta[['exp_design']]$pair2))
  print(all_fastq)
  message("Processing...")
  for(f in all_fastq){
    message(paste0("\t", f))
    #Check if this system command should be executed
    if(execute){
      temp_out = system(paste0(
        "fastqc -o ", 
        fastqc_outpath,
        " --extract ",
        " -t ",
        num_cores,
        " ",
        f,
        " > ", 
        file.path(fastqc_outpath, paste0(basename(f), '.stdout')),
        " 2> ",
        file.path(fastqc_outpath, paste0(basename(f), '.stderr'))))
    }
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Update fastq metadata
  dimsum_meta_new[['exp_design']]$pair1_fastqc <- gsub(dimsum_meta_new[["fastq_file_extension"]], '_fastqc/fastqc_data.txt', gsub('.gz', '', dimsum_meta_new[['exp_design']][,"pair1"]))
  dimsum_meta_new[['exp_design']]$pair2_fastqc <- gsub(dimsum_meta_new[["fastq_file_extension"]], '_fastqc/fastqc_data.txt', gsub('.gz', '', dimsum_meta_new[['exp_design']][,"pair2"]))
  dimsum_meta_new[['exp_design']]$fastqc_directory <- fastqc_outpath
  #Generate FASTQC report
  if(report){
    dimsum_meta_new_report <- dimsum_stage_fastqc_report(dimsum_meta_new, report_outpath)
    return(dimsum_meta_new_report)
  }else{
    return(dimsum_meta_new)
  }
}

#dimsum_stage_unzip
#
# Unzip all fastq files.
#
# dimsum_meta: an experiment metadata object (required)
# fastq_outpath: FASTQ output path (required)
# execute: whether or not to execute the system command (default: TRUE)
#
# Returns: an updated experiment metadata object.
#
dimsum_stage_unzip <- function(
  dimsum_meta,
  fastq_outpath,
  execute = TRUE
  ){
  #All fastq files gzipped?
  if(dimsum_meta[["gzipped"]]){
    #Create unzip directory (if doesn't already exist)
    fastq_outpath <- gsub("/$", "", fastq_outpath)
    dir.create(fastq_outpath)
    message("Unzipping FASTQ files:")
    all_fastq <- file.path(dimsum_meta[["exp_design"]]$pair_directory, c(dimsum_meta[['exp_design']]$pair1, dimsum_meta[['exp_design']]$pair2))
    print(all_fastq)
    message("Processing...")
    for(f in all_fastq){
      message(paste0("\t", f))
      #Check if this system command should be executed
      if(execute){
        temp_out = system(paste0(
          "gunzip -c ", 
          f, 
          " > ", 
          file.path(fastq_outpath, gsub(".gz$", "", basename(f))),
          " 2> ",
          file.path(fastq_outpath, paste0(gsub(".gz$", "", basename(f)), '.stderr'))))
      }
    }
    #New experiment metadata
    dimsum_meta_new <- dimsum_meta
    #Update fastq metadata
    dimsum_meta_new[['exp_design']]$pair1 <- gsub(".gz$", "", dimsum_meta_new[["exp_design"]][,"pair1"])
    dimsum_meta_new[['exp_design']]$pair2 <- gsub(".gz$", "", dimsum_meta_new[["exp_design"]][,"pair2"])
    dimsum_meta_new[['exp_design']]$pair_directory <- fastq_outpath
    return(dimsum_meta_new)
  }else{
    #Copy fastq files
    message("FASTQ files already uncompressed.")
    #New experiment metadata
    dimsum_meta_new <- dimsum_meta
    return(dimsum_meta_new)
  }
}

#dimsum_stage_split
#
# Split all fastq files.
#
# dimsum_meta: an experiment metadata object (required)
# split_outpath: split FASTQ output path (required)
# execute: whether or not to execute the system command (default: TRUE)
#
# Returns: an updated experiment metadata object.
#
dimsum_stage_split <- function(
  dimsum_meta,
  split_outpath,
  execute = TRUE
  ){
  #Create unzip directory (if doesn't already exist)
  split_outpath <- gsub("/$", "", split_outpath)
  dir.create(split_outpath)
  fastq_pair_list <- dimsum_meta[['exp_design']][,c('pair1', 'pair2')]
  rownames(fastq_pair_list) = 1:dim(fastq_pair_list)[1]
  #Split FASTQ files
  message("Splitting FASTQ files:")
  all_fastq <- file.path(dimsum_meta[["exp_design"]]$pair_directory, c(dimsum_meta[['exp_design']]$pair1, dimsum_meta[['exp_design']]$pair2))
  print(all_fastq)
  message("Processing...")
  for(pair_name in rownames(fastq_pair_list)){
    print(fastq_pair_list[pair_name,])
    #Check if this system command should be executed
    if(execute){
      temp_out = system(paste0(
        "fastq_splitter.py -i ", 
        file.path(dimsum_meta[["exp_design"]]$pair_directory, fastq_pair_list[pair_name,][1]), 
        " -o ", 
        file.path(split_outpath, paste0(fastq_pair_list[pair_name,][1], ".split")), 
        " -c 3758096384",
        " > ",
        file.path(split_outpath, paste0(gsub(dimsum_meta[["fastq_file_extension"]], '', fastq_pair_list[pair_name,][1]), ".split.stdout")),
        " 2> ",
        file.path(split_outpath, paste0(gsub(dimsum_meta[["fastq_file_extension"]], '', fastq_pair_list[pair_name,][1]), ".split.stderr"))))
      num_records = as.integer(read.table(file.path(split_outpath, paste0(gsub(dimsum_meta[["fastq_file_extension"]], '', fastq_pair_list[pair_name,][1]), ".split.stdout"))))
      temp_out = system(paste0(
        "fastq_splitter.py -i ", 
        file.path(dimsum_meta[["exp_design"]]$pair_directory, fastq_pair_list[pair_name,][2]), 
        " -o ", 
        file.path(split_outpath, paste0(fastq_pair_list[pair_name,][2], ".split")), 
        " -n ", 
        num_records,
        " >> ",
        file.path(split_outpath, paste0(gsub(dimsum_meta[["fastq_file_extension"]], '', fastq_pair_list[pair_name,][2]), ".split.stdout")),
        " 2>> ",
        file.path(split_outpath, paste0(gsub(dimsum_meta[["fastq_file_extension"]], '', fastq_pair_list[pair_name,][2]), ".split.stderr"))))
    }
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Update fastq metadata
  all_fastq <- file.path(dimsum_meta_new[["exp_design"]]$pair_directory, c(dimsum_meta_new[['exp_design']]$pair1))
  split_list <- list()
  for(f in all_fastq){
    num_splits <- length(list.files(split_outpath, pattern = basename(f)))
    split_list = append(split_list, num_splits)
  }
  dimsum_meta_new[["exp_design"]] = dimsum_meta_new[["exp_design"]][rep(1:length(all_fastq), times = unlist(split_list)),]
  temp_rownames = rownames(dimsum_meta_new[["exp_design"]])
  temp_suffix = rep('.split1', dim(dimsum_meta_new[["exp_design"]])[1])
  temp_suffix[grepl('\\.', temp_rownames)] = paste0('.split', as.integer(sapply(strsplit(temp_rownames[grepl('\\.', temp_rownames)], '\\.'), '[', 2))+1)
  dimsum_meta_new[["exp_design"]]$pair1 = paste0(dimsum_meta_new[["exp_design"]]$pair1, temp_suffix, '.fastq')
  dimsum_meta_new[["exp_design"]]$pair2 = paste0(dimsum_meta_new[["exp_design"]]$pair2, temp_suffix, '.fastq')
  dimsum_meta_new[["exp_design"]]$split = as.integer(gsub(".split", "", temp_suffix))
  dimsum_meta_new[['exp_design']]$pair_directory <- split_outpath
  return(dimsum_meta_new)
}

#dimsum_stage_cutadapt_report
#
# Generate cutadapt summary plots for all fastq files.
#
# dimsum_meta: an experiment metadata object (required)
# report_outpath: cutadapt report output path (required)
#
# Returns: an updated experiment metadata object
#
dimsum_stage_cutadapt_report <- function(
  dimsum_meta,
  report_outpath
  ){
  #Create report directory (if doesn't already exist)
  report_outpath <- gsub("/$", "", report_outpath)
  dir.create(report_outpath)
  #Get cutadapt results for all read pairs
  cutadapt_files <- file.path(dimsum_meta[['exp_design']]$pair_directory, paste0(dimsum_meta[['exp_design']][,'pair1'], '.stdout'))
  cutadapt_read1_list <- list()
  cutadapt_read2_list <- list()
  total_reads_list <- list()
  for(i in 1:length(cutadapt_files)){
    temp_out <- system(paste0("cat ", cutadapt_files[i]), intern=TRUE)
    name_read1 <- basename(rev(unlist(strsplit(temp_out[2], ' ')))[2])
    name_read2 <- basename(rev(unlist(strsplit(temp_out[2], ' ')))[1])
    temp_out <- temp_out[c(9:13, grep("Sequence: ", temp_out))]
    total_reads <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[1], ' ')))[1]))
    total_reads_list[[i]] <- total_reads
    total_read1_trimmed <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[2], ' ')))[2]))
    total_read2_trimmed <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[3], ' ')))[2]))
    total_read1_a3 <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[6], ' ')))[2]))
    total_read1_a5 <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[7], ' ')))[2]))
    total_read1_both <- total_read1_a3+total_read1_a5-total_read1_trimmed
    total_read2_a3 <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[8], ' ')))[2]))
    total_read2_a5 <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[9], ' ')))[2]))
    total_read2_both <- total_read2_a3+total_read2_a5-total_read2_trimmed
    #If linked adapter supplied
    if(grepl('Type: linked', temp_out[6]) & grepl('Type: linked', temp_out[7])){
      total_read1_a3 <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[6], ' ')))[2]))
      total_read1_a5 <- total_read1_a3
      total_read1_both <- total_read1_a3
      total_read2_a3 <- as.integer(gsub(',', '', rev(unlist(strsplit(temp_out[7], ' ')))[2]))
      total_read2_a5 <- total_read2_a3
      total_read2_both <- total_read2_a3
    }
    cutadapt_read1_list[[name_read1]] <- c(
      total_read1_a5-total_read1_both, total_read1_a3-total_read1_both, total_read1_both, total_reads)
    cutadapt_read2_list[[name_read2]] <- c(
      total_read2_a5-total_read2_both, total_read2_a3-total_read2_both, total_read2_both, total_reads)
  }
  #First read
  cutadapt_read1_df <- as.data.frame(do.call('rbind', cutadapt_read1_list))
  colnames(cutadapt_read1_df) <- c('adapter5prime', 'adapter3prime', 'both', 'total_reads')
  cutadapt_read1_df$fastq <- sapply(strsplit(rownames(cutadapt_read1_df), '.split'), '[', 1)
  cutadapt_read1_df_collapse <- ddply(cutadapt_read1_df, "fastq", summarise, 
    adapter5prime = sum(adapter5prime), 
    adapter3prime = sum(adapter3prime), 
    both = sum(both), 
    total_reads = sum(total_reads))
  cutadapt_read1_df_collapse_perc <- cutadapt_read1_df_collapse
  cutadapt_read1_df_collapse_perc[,2:4] <- as.data.frame(t(scale(t(cutadapt_read1_df_collapse_perc[,2:4]), center = F, scale = cutadapt_read1_df_collapse_perc$total_reads)))*100
  cutadapt_read1_df_collapse_perc <- cutadapt_read1_df_collapse_perc[,1:4]
  #Plot
  plot_df <- melt(cutadapt_read1_df_collapse_perc, id="fastq")
  d <- ggplot(plot_df, aes(fastq, value)) +
    geom_col(aes(fill = variable)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "FASTQ files", y = "Reads with adapters (percentage)", title = paste0("Percentage first read trimmed with cutadapt"))
  ggsave(file.path(report_outpath, paste0('dimsum_stage_cutadapt_report_pair1.png')), d, width=12, height=8)
  #Second read
  cutadapt_read2_df <- as.data.frame(do.call('rbind', cutadapt_read2_list))
  colnames(cutadapt_read2_df) <- c('adapter5prime', 'adapter3prime', 'both', 'total_reads')
  cutadapt_read2_df$fastq <- sapply(strsplit(rownames(cutadapt_read2_df), '.split'), '[', 1)
  cutadapt_read2_df_collapse <- ddply(cutadapt_read2_df, "fastq", summarise, 
    adapter5prime = sum(adapter5prime), 
    adapter3prime = sum(adapter3prime), 
    both = sum(both), 
    total_reads = sum(total_reads))
  cutadapt_read2_df_collapse_perc <- cutadapt_read2_df_collapse
  cutadapt_read2_df_collapse_perc[,2:4] <- as.data.frame(t(scale(t(cutadapt_read2_df_collapse_perc[,2:4]), center = F, scale = cutadapt_read2_df_collapse_perc$total_reads)))*100
  cutadapt_read2_df_collapse_perc <- cutadapt_read2_df_collapse_perc[,1:4]
  #Plot
  plot_df <- melt(cutadapt_read2_df_collapse_perc, id="fastq")
  d <- ggplot(plot_df, aes(fastq, value)) +
    geom_col(aes(fill = variable)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "FASTQ files", y = "Reads with adapters (percentage)", title = paste0("Percentage second read trimmed with cutadapt"))
  ggsave(file.path(report_outpath, paste0('dimsum_stage_cutadapt_report_pair2.png')), d, width=12, height=8)
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Update fastq metadata
  dimsum_meta_new[['exp_design']]$total_read_pairs <- as.numeric(unlist(total_reads_list))
  return(dimsum_meta_new)
}

#dimsum_stage_cutadapt
#
# Run cutadapt on all fastq files.
#
# dimsum_meta: an experiment metadata object (required)
# cutadapt_outpath: cutadapt output path (required)
# execute: whether or not to execute the system command (default: TRUE)
#
# Returns: an updated experiment metadata object.
#
dimsum_stage_cutadapt <- function(
  dimsum_meta,
  cutadapt_outpath,
  execute = TRUE,
  report = TRUE,
  report_outpath = NULL
  ){
  #Create cutadapt directory (if doesn't already exist)
  cutadapt_outpath <- gsub("/$", "", cutadapt_outpath)
  dir.create(cutadapt_outpath)
  fastq_pair_list <- dimsum_meta[['exp_design']][,c('pair1', 'pair2')]
  rownames(fastq_pair_list) = 1:dim(fastq_pair_list)[1]
  #Cutadapt parameters specified?
  if( is.null(dimsum_meta[["cutadapt5First"]]) | is.null(dimsum_meta[["cutadapt5Second"]]) ){
    message("Skipping cutadapt adapter removal (all catadapt arguments need to be specified)")
    return(dimsum_meta)
  }else{
    #Options for removing constant regions from beginning or end of either read in pair
    temp_options = paste0(' -g ', dimsum_meta[["cutadapt5First"]], ' -G ', dimsum_meta[["cutadapt5Second"]])
    if( !is.null(dimsum_meta[["cutadapt3First"]]) ){temp_options = paste0(temp_options, " -a ", dimsum_meta[["cutadapt3First"]])}
    if( !is.null(dimsum_meta[["cutadapt3Second"]]) ){temp_options = paste0(temp_options, " -A ", dimsum_meta[["cutadapt3Second"]])}
    if( dimsum_meta[["cutadaptDiscardUntrimmed"]] ){temp_options = paste0(temp_options, " --discard-untrimmed ")}
    #Options for swapping read1 and read2
    temp_options_swap = paste0(' -g forward=', dimsum_meta[["cutadapt5First"]], ' -g reverse=', dimsum_meta[["cutadapt5Second"]])
    #Options for removing a fixed number of bases from beginning or end of either read in pair
    temp_cut_options = ''
    if( !is.null(dimsum_meta[["cutadaptCut5First"]]) ){temp_cut_options = paste0(temp_cut_options, " -u ", dimsum_meta[["cutadaptCut5First"]])}
    if( !is.null(dimsum_meta[["cutadaptCut3First"]]) ){temp_cut_options = paste0(temp_cut_options, " -u ", -dimsum_meta[["cutadaptCut3First"]])}
    if( !is.null(dimsum_meta[["cutadaptCut5Second"]]) ){temp_cut_options = paste0(temp_cut_options, " -U ", dimsum_meta[["cutadaptCut5Second"]])}
    if( !is.null(dimsum_meta[["cutadaptCut3Second"]]) ){temp_cut_options = paste0(temp_cut_options, " -U ", -dimsum_meta[["cutadaptCut3Second"]])}
    #Trim FASTQ file pairs
    message("Trimming FASTQ files with cutadapt:")
    all_fastq <- file.path(dimsum_meta[["exp_design"]]$pair_directory, c(dimsum_meta[['exp_design']]$pair1, dimsum_meta[['exp_design']]$pair2))
    print(all_fastq)
    message("Processing...")
    for(pair_name in rownames(fastq_pair_list)){
      #TODO: cutadapt binary path specifiable on commandline?
      #TEMP: don't trim files
      print(fastq_pair_list[pair_name,])
      #Check if this system command should be executed
      if(execute){
        #Not stranded library
        if( !dimsum_meta[["stranded"]] ){
          #Swap reads in fastq files according to adapter presence
          temp_out = system(paste0(
            "cutadapt",
            temp_options_swap,
            temp_cut_options,
            " --no-trim ",
            " --minimum-length ",
            as.character(dimsum_meta[["cutadaptMinLength"]]),
            " -e ",
            as.character(dimsum_meta[["cutadaptErrorRate"]]),
            " --untrimmed-output ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt.untrimmed.fastq")),
            " --untrimmed-paired-output ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][2], ".cutadapt.untrimmed.fastq")),
            " -o ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt1-{name}.fastq")),
            " -p ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][2], ".cutadapt1-{name}.fastq")),
            " ",
            file.path(dimsum_meta[["exp_design"]]$pair_directory, fastq_pair_list[pair_name,][1]),
            " ",
            file.path(dimsum_meta[["exp_design"]]$pair_directory, fastq_pair_list[pair_name,][2]),
            " > ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt1.stdout")),
            " 2> ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt1.stderr"))))
          #New read1 file
          temp_out = system(paste0(
            "cat ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][2], ".cutadapt1-reverse.fastq")),
            " ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt1-forward.fastq")),
            " ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt.untrimmed.fastq")),
            " > ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt2"))))
          #New read2 file
          temp_out = system(paste0(
            "cat ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt1-reverse.fastq")),
            " ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][2], ".cutadapt1-forward.fastq")),
            " ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][2], ".cutadapt.untrimmed.fastq")),
            " > ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][2], ".cutadapt2"))))
          #Run cutadapt on the swapped  
          temp_out = system(paste0(
            "cutadapt",
            temp_options,
            " --minimum-length ",
            as.character(dimsum_meta[["cutadaptMinLength"]]),
            " -e ",
            as.character(dimsum_meta[["cutadaptErrorRate"]]),
            " -j ",
            num_cores,
            " -o ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt")),
            " -p ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][2], ".cutadapt")),
            " ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt2")),
            " ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][2], ".cutadapt2")),
            " > ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt.stdout")),
            " 2> ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt.stderr")))) 
        #Stranded library   
        }else{
          temp_out = system(paste0(
            "cutadapt",
            temp_options,
            temp_cut_options,
            " --minimum-length ",
            as.character(dimsum_meta[["cutadaptMinLength"]]),
            " -e ",
            as.character(dimsum_meta[["cutadaptErrorRate"]]),
            " -j ",
            num_cores,
            " -o ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt")),
            " -p ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][2], ".cutadapt")),
            " ",
            file.path(dimsum_meta[["exp_design"]]$pair_directory, fastq_pair_list[pair_name,][1]),
            " ",
            file.path(dimsum_meta[["exp_design"]]$pair_directory, fastq_pair_list[pair_name,][2]),
            " > ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt.stdout")),
            " 2> ",
            file.path(cutadapt_outpath, paste0(fastq_pair_list[pair_name,][1], ".cutadapt.stderr"))))
        }
      }
    }
    #New experiment metadata
    dimsum_meta_new <- dimsum_meta
    #Update fastq metadata
    temp_suffix <- ".cutadapt"
    dimsum_meta_new[["exp_design"]]$pair1 = paste0(dimsum_meta_new[["exp_design"]]$pair1, temp_suffix)
    dimsum_meta_new[["exp_design"]]$pair2 = paste0(dimsum_meta_new[["exp_design"]]$pair2, temp_suffix)
    dimsum_meta_new[['exp_design']]$pair_directory <- cutadapt_outpath
    #Generate cutadapt report
    if(report){
      dimsum_meta_new_report <- dimsum_stage_cutadapt_report(dimsum_meta_new, report_outpath)
      return(dimsum_meta_new_report)
    }else{
      return(dimsum_meta_new)
    }
  }
}

#dimsum_stage_usearch_report
#
# Generate USEARCH summary plots for all samples.
#
# dimsum_meta: an experiment metadata object (required)
# report_outpath: USEARCH report output path (required)
#
# Returns: an updated experiment metadata object
#
dimsum_stage_usearch_report <- function(
  dimsum_meta,
  report_outpath
  ){
  #Create report directory (if doesn't already exist)
  report_outpath <- gsub("/$", "", report_outpath)
  dir.create(report_outpath)
  #Get cutadapt results for all read pairs
  usearch_files <- file.path(dimsum_meta[['exp_design']]$aligned_pair_directory, gsub('.usearch$', '.report', dimsum_meta[['exp_design']][,'aligned_pair']))
  usearch_list <- list()
  for(i in 1:length(usearch_files)){
    temp_out <- system(paste0("cat ", usearch_files[i]), intern=TRUE)
    usearch_list[[i]] <- list()
    usearch_list[[i]][['usearch_merge_length_low_quartile']] <- as.integer(rev(unlist(strsplit(temp_out[grep('Low quartile', temp_out)], ' ')))[4])
    usearch_list[[i]][['usearch_merge_length_median']] <- as.integer(rev(unlist(strsplit(temp_out[grep('Median', temp_out)], ' ')))[3])
    usearch_list[[i]][['usearch_merge_length_high_quartile']] <- as.integer(rev(unlist(strsplit(temp_out[grep('High quartile', temp_out)], ' ')))[4])
    usearch_list[[i]][['usearch_total_read_pairs']] <- as.integer(rev(unlist(strsplit(temp_out[grep('Pairs', temp_out)], ' ')))[4])
    usearch_list[[i]][['usearch_merged']] <- as.integer(rev(unlist(strsplit(temp_out[grep('Merged', temp_out)], ' ')))[5])
    usearch_list[[i]][['usearch_too_many_diffs']] <- as.integer(rev(unlist(strsplit(temp_out[grep('Too many diffs', temp_out)], ' ')))[8])
    usearch_list[[i]][['usearch_fwd_too_short']] <- as.integer(rev(unlist(strsplit(temp_out[grep('Fwd too short', temp_out)], ' ')))[11])
    usearch_list[[i]][['usearch_rev_too_short']] <- as.integer(rev(unlist(strsplit(temp_out[grep('Rev too short', temp_out)], ' ')))[11])
    usearch_list[[i]][['usearch_no_alignment_found']] <- as.integer(rev(unlist(strsplit(temp_out[grep('No alignment found', temp_out)], ' ')))[6])
    usearch_list[[i]][['usearch_alignment_too_short']] <- as.integer(rev(unlist(strsplit(temp_out[grep('Alignment too short', temp_out)], ' ')))[8])
    usearch_list[[i]][['usearch_exp_errs_too_high']] <- as.integer(rev(unlist(strsplit(temp_out[grep('Exp.errs. too high', temp_out)], ' ')))[7])
    usearch_list[[i]][['usearch_min_Q_too_low']] <- as.integer(rev(unlist(strsplit(temp_out[grep('Min Q too low', temp_out)], ' ')))[8])
  }
  usearch_df <- as.data.frame(apply(as.data.frame(do.call('rbind', usearch_list)), 2, unlist))
  #Merge experimental design with USEARCH report statistics
  usearch_df <- cbind(dimsum_meta[['exp_design']][,c('aligned_pair', 'total_read_pairs')], usearch_df)
  usearch_df$cutadapt_pairs_too_short <- usearch_df$total_read_pairs - usearch_df$usearch_total_read_pairs
  #Plot 1: read pair count statistics
  usearch_df$pairname <- sapply(strsplit(usearch_df$aligned_pair, '.split'), '[', 1)
  usearch_df_collapse <- ddply(usearch_df, "pairname", summarise, 
    total_read_pairs = sum(total_read_pairs), 
    usearch_merged = sum(usearch_merged), 
    usearch_too_many_diffs = sum(usearch_too_many_diffs), 
    usearch_fwd_too_short = sum(usearch_fwd_too_short), 
    usearch_rev_too_short = sum(usearch_rev_too_short), 
    usearch_no_alignment_found = sum(usearch_no_alignment_found), 
    usearch_alignment_too_short = sum(usearch_alignment_too_short), 
    usearch_exp_errs_too_high = sum(usearch_exp_errs_too_high),
    usearch_min_Q_too_low = sum(usearch_min_Q_too_low),
    cutadapt_pairs_too_short = sum(cutadapt_pairs_too_short)
    )
  usearch_df_collapse_perc <- usearch_df_collapse
  usearch_df_collapse_perc[,3:11] <- as.data.frame(t(scale(t(usearch_df_collapse_perc[,3:11]), center = F, scale = usearch_df_collapse_perc$total_read_pairs)))*100
  usearch_df_collapse_perc <- usearch_df_collapse_perc[,c(1, 3:11)]
  #Plot
  plot_df <- melt(usearch_df_collapse_perc, id="pairname")
  d <- ggplot(plot_df, aes(pairname, value)) +
    geom_col(aes(fill = variable)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Sample names", y = "Pairs (percentage)", title = paste0("Read pair alignment statistics"))
  ggsave(file.path(report_outpath, paste0('dimsum_stage_usearch_report_paircounts.png')), d, width=12, height=8)
  #Plot2: read pair merge length statistics
  plot_df <- melt(usearch_df[,grep('pairname|usearch_merge_', colnames(usearch_df))], id="pairname")
  d <- ggplot(plot_df, aes(pairname, value)) +
    geom_boxplot(aes(color = variable)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    coord_cartesian(ylim = c(0, max(plot_df$value))) +
    labs(x = "Sample names", y = "Merged length", title = paste0("Alignment length distributions (over splits)"))
  ggsave(file.path(report_outpath, paste0('dimsum_stage_usearch_report_mergedlength.png')), d, width=12, height=8)
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Update fastq metadata
  dimsum_meta_new[['exp_design']] <- cbind(dimsum_meta_new[['exp_design']], usearch_df[,3:15])
  return(dimsum_meta_new)
}

#dimsum_stage_usearch
#
# Run USEARCH on all fastq files.
#
# dimsum_meta: an experiment metadata object (required)
# usearch_outpath: USEARCH output path (required)
# execute: whether or not to execute the system command (default: TRUE)
#
# Returns: an updated experiment metadata object.
#
dimsum_stage_usearch <- function(
  dimsum_meta,
  usearch_outpath,
  execute = TRUE,
  report = TRUE,
  report_outpath = NULL
  ){
  #Create cutadapt directory (if doesn't already exist)
  usearch_outpath <- gsub("/$", "", usearch_outpath)
  dir.create(usearch_outpath)
  #Sample names
  sample_names = paste0(
    dimsum_meta[["exp_design"]]$sample_name, '_e', 
    dimsum_meta[["exp_design"]]$experiment, '_s', 
    dimsum_meta[["exp_design"]]$selection_id, '_b', 
    dimsum_meta[["exp_design"]]$biological_replicate, '_t', 
    dimsum_meta[["exp_design"]]$technical_replicate, '_split', 
    dimsum_meta[["exp_design"]]$split, sep = "")
  #Additional usearch options related to alignment length
  temp_options = paste0(' -fastq_minovlen ', dimsum_meta[["usearchMinovlen"]])
  if( dimsum_meta[["usearchMinovlen"]] < 16 ){temp_options = paste0(temp_options, " -xdrop_nw ", dimsum_meta[["usearchMinovlen"]])}
  if( dimsum_meta[["usearchMinovlen"]] < 16 ){temp_options = paste0(temp_options, " -minhsp ", dimsum_meta[["usearchMinovlen"]])}
  if( dimsum_meta[["usearchMinovlen"]] < 16 ){temp_options = paste0(temp_options, " -band ", dimsum_meta[["usearchMinovlen"]])}
  if( dimsum_meta[["usearchMinovlen"]] < 5 ){temp_options = paste0(temp_options, " -hspw ", dimsum_meta[["usearchMinovlen"]])}
  #Run USEARCH on all fastq file pairs
  message("Merging paired-end FASTQ files with USEARCH:")
  all_fastq <- file.path(dimsum_meta[["exp_design"]]$pair_directory, c(dimsum_meta[['exp_design']]$pair1, dimsum_meta[['exp_design']]$pair2))
  print(all_fastq)
  message("Processing...")
  for(i in 1:length(sample_names)){
    #TODO: usearch binary path specifiable on commandline?
    #TODO: only run if usearch arguments specified
    print(dimsum_meta[["exp_design"]][i,c('pair1', 'pair2')])
    #Check if this system command should be executed
    if(execute){
      temp_out = system(paste0(
        "usearch -fastq_mergepairs ",
        file.path(dimsum_meta[["exp_design"]]$pair_directory, dimsum_meta[["exp_design"]][i,]$pair1),
        " -reverse ",
        file.path(dimsum_meta[["exp_design"]]$pair_directory, dimsum_meta[["exp_design"]][i,]$pair2),
        " -fastqout ",
        file.path(usearch_outpath, paste0(sample_names[i], '.usearch')),
        " -report ",
        file.path(usearch_outpath, paste0(sample_names[i], '.report')),
        " -fastq_minqual ",
        as.character(dimsum_meta[["usearchMinQual"]]),
        " -fastq_merge_maxee ",
        as.character(dimsum_meta[["usearchMaxee"]]),
        temp_options,
        " -threads ",
        num_cores,
        " > ",
        file.path(usearch_outpath, paste0(sample_names[i], '.usearch.stdout')),
        " 2> ",
        file.path(usearch_outpath, paste0(sample_names[i], '.usearch.stderr'))))
    }
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Merged fastq filenames
  dimsum_meta_new[["exp_design"]]$aligned_pair <- paste0(sample_names, ".usearch")
  dimsum_meta_new[['exp_design']]$aligned_pair_directory <- usearch_outpath
  #Generate usearch report
  if(report){
    dimsum_meta_new_report <- dimsum_stage_usearch_report(dimsum_meta_new, report_outpath)
    return(dimsum_meta_new_report)
  }else{
    return(dimsum_meta_new)
  }
}

#dimsum_stage_unique
#
# Run fastx_collapser on all fastq files.
#
# dimsum_meta: an experiment metadata object (required)
# unique_outpath: fastx_collapser output path (required)
# execute: whether or not to execute the system command (default: TRUE)
#
# Returns: an updated experiment metadata object.
#
dimsum_stage_unique <- function(
  dimsum_meta,
  unique_outpath,
  execute = TRUE
  ){
  #Create cutadapt directory (if doesn't already exist)
  unique_outpath <- gsub("/$", "", unique_outpath)
  dir.create(unique_outpath)
  #Run fastx_collapser on all aligned read pair fastq files
  message("Getting unique aligned read counts with fastx_collapser:")
  all_fasta <- file.path(dimsum_meta[["exp_design"]]$aligned_pair_directory, dimsum_meta[['exp_design']]$aligned_pair)
  print(all_fasta)
  message("Processing...")
  for(read_pair in dimsum_meta[["exp_design"]]$aligned_pair){
    #TODO: usearch binary path specifiable on commandline?
    #TODO: only run if usearch arguments specified
    print(read_pair)
    #Check if this system command should be executed
    if(execute){
      temp_out = system(paste0(
        "fastx_collapser -i ",
        file.path(dimsum_meta[["exp_design"]]$aligned_pair_directory, read_pair),
        " -o ",
        file.path(unique_outpath, paste0(read_pair, '.unique')),
        " > ",
        file.path(unique_outpath, paste0(read_pair, '.unique.stdout')),
        " 2> ",
        file.path(unique_outpath, paste0(read_pair, '.unique.stderr'))))
    }
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Unique fasta filenames
  dimsum_meta_new[["exp_design"]]$aligned_pair_unique <- paste0(dimsum_meta_new[["exp_design"]]$aligned_pair, ".unique")
  dimsum_meta_new[['exp_design']]$aligned_pair_unique_directory <- unique_outpath
  return(dimsum_meta_new)
}

#dimsum_stage_filter
#
# Run fastx_collapser2AAtable.py on all input files.
#
# dimsum_meta: an experiment metadata object (required)
# filter_outpath: fastx_collapser2AAtable.py output path (required)
# execute: whether or not to execute the system command (default: TRUE)
#
# Returns: an updated experiment metadata object.
#
dimsum_stage_filter <- function(
  dimsum_meta,
  filter_outpath,
  execute = TRUE
  ){
  #Create cutadapt directory (if doesn't already exist)
  filter_outpath <- gsub("/$", "", filter_outpath)
  dir.create(filter_outpath)
  #Construct filtered variant count table with corresponding amino acid sequences
  message("Constructing filtered variant count table with corresponding amino acid sequences")
  all_fasta <- file.path(dimsum_meta[["exp_design"]]$aligned_pair_unique_directory, dimsum_meta[['exp_design']]$aligned_pair_unique)
  print(all_fasta)
  message("Processing...")
  for(read_pair in dimsum_meta[['exp_design']]$aligned_pair_unique){
    print(read_pair)
    #Check if this system command should be executed
    if(execute){
      temp_out = system(paste0(
        "fastx_collapser2AAtable.py -i ", 
        file.path(dimsum_meta[["exp_design"]]$aligned_pair_unique_directory, read_pair), 
        " -o ", 
        file.path(filter_outpath, paste0(read_pair, ".tsv")), 
        " -n ", 
        nchar(dimsum_meta[['wildtypeSequence']]),
        " > ",
        file.path(filter_outpath, paste0(read_pair, '.tsv.stdout')),
        " 2> ",
        file.path(filter_outpath, paste0(read_pair, '.tsv.stderr'))))
    }
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Filtered fasta filenames
  dimsum_meta_new[["exp_design"]]$aligned_pair_unique_tsv <- paste0(dimsum_meta_new[["exp_design"]]$aligned_pair_unique, ".tsv")
  dimsum_meta_new[['exp_design']]$aligned_pair_unique_tsv_directory <- filter_outpath
  #Get filter results for all samples
  filter_files <- file.path(dimsum_meta_new[['exp_design']]$aligned_pair_unique_tsv_directory, gsub('.tsv$', '.tsv.stdout', dimsum_meta_new[['exp_design']][,'aligned_pair_unique_tsv']))
  filter_list <- list()
  for(i in 1:length(filter_files)){
    temp_out <- system(paste0("cat ", filter_files[i]), intern=TRUE)
    filter_list[[i]] <- as.integer(rev(unlist(strsplit(temp_out[grep('Sequences failed', temp_out)], '\\t')))[1])
  }
  dimsum_meta_new[['exp_design']]$filter_incorrect_length = unlist(filter_list)
  return(dimsum_meta_new)
}

#dimsum_stage_merge_report
#
# Generate final summary plots for all samples.
#
# dimsum_meta: an experiment metadata object (required)
# report_outpath: final report output path (required)
#
# Returns: an updated experiment metadata object
#
dimsum_stage_merge_report <- function(
  dimsum_meta,
  report_outpath
  ){
  #Create report directory (if doesn't already exist)
  report_outpath <- gsub("/$", "", report_outpath)
  dir.create(report_outpath)
  #Final statistics
  merge_df <- dimsum_meta[['exp_design']]
  merge_df$pairname <- sapply(strsplit(merge_df$aligned_pair, '.split'), '[', 1)
  merge_df_collapse <- ddply(merge_df, "pairname", summarise, 
    total_read_pairs = sum(total_read_pairs), 
    usearch_merged = sum(usearch_merged), 
    cutadapt_pairs_too_short = sum(cutadapt_pairs_too_short),
    filter_incorrect_length = sum(filter_incorrect_length),
    too_many_aa_mutations = sum(too_many_aa_mutations))
  merge_df_collapse$usearch_not_merged <- merge_df_collapse$total_read_pairs-merge_df_collapse$cutadapt_pairs_too_short-merge_df_collapse$usearch_merged
  merge_df_collapse$retained <- merge_df_collapse$total_read_pairs-merge_df_collapse$cutadapt_pairs_too_short-merge_df_collapse$usearch_not_merged-merge_df_collapse$filter_incorrect_length-merge_df_collapse$too_many_aa_mutations
  merge_df_collapse_perc <- merge_df_collapse
  merge_df_collapse_perc[,4:8] <- merge_df_collapse_perc[,4:8]/merge_df_collapse_perc$total_read_pairs*100
  merge_df_collapse_perc <- merge_df_collapse_perc[,c(1,4:8)]
  merge_df_collapse <- merge_df_collapse[,c(1,4:8)]
  #Plot 1: final results for all samples (total reads)
  plot_df <- melt(merge_df_collapse, id="pairname")
  d <- ggplot(plot_df, aes(pairname, value)) +
    geom_col(aes(fill = variable)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Sample names", y = "Total read pairs (variants)", title = paste0("Read pairs to filtered variant statistics"))
  ggsave(file.path(report_outpath, paste0('dimsum_stage_merge_report_variantcounts.png')), d, width=12, height=8)
  #Plot 2: final results for all samples (percentages)
  plot_df <- melt(merge_df_collapse_perc, id="pairname")
  d <- ggplot(plot_df, aes(pairname, value)) +
    geom_col(aes(fill = variable)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Sample names", y = "Percentage of initial read pairs (variants)", title = paste0("Read pairs to filtered variant statistics (percentage)"))
  ggsave(file.path(report_outpath, paste0('dimsum_stage_merge_report_variantpercentages.png')), d, width=12, height=8)
  #Plot 3: AA mutation counts
  aa_mut_df <- data.frame(
    'AA_mut_0'=sapply(dimsum_meta[['aa_mutation_counts']], '[', '0'),
    'AA_mut_1'=sapply(dimsum_meta[['aa_mutation_counts']], '[', '1'),
    'AA_mut_2'=sapply(dimsum_meta[['aa_mutation_counts']], '[', '2'),
    'AA_mut_sum'=sapply(dimsum_meta[['aa_mutation_counts']], sum))
  aa_mut_df[is.na(aa_mut_df)] <- 0
  aa_mut_df$pairname <- sapply(strsplit(merge_df$aligned_pair, '.split'), '[', 1)
  aa_mut_df_collapse <- ddply(aa_mut_df, "pairname", summarise, 
    AA_mut_0 = sum(AA_mut_0), 
    AA_mut_1 = sum(AA_mut_1), 
    AA_mut_2 = sum(AA_mut_2),
    AA_mut_sum = sum(AA_mut_sum))
  aa_mut_df_collapse$AA_mut_3plus = aa_mut_df_collapse$AA_mut_sum-aa_mut_df_collapse$AA_mut_0-aa_mut_df_collapse$AA_mut_1-aa_mut_df_collapse$AA_mut_2
  aa_mut_df_collapse_perc = aa_mut_df_collapse
  aa_mut_df_collapse_perc[,c(2:5, 6)] <- aa_mut_df_collapse_perc[,c(2:5, 6)]/aa_mut_df_collapse_perc$AA_mut_sum*100
  aa_mut_df_collapse_perc <- aa_mut_df_collapse_perc[,c(1:4, 6)]
  aa_mut_df_collapse <- aa_mut_df_collapse[,c(1:4, 6)]
  #Plot
  plot_df <- melt(aa_mut_df_collapse, id="pairname")
  d <- ggplot(plot_df, aes(pairname, value)) +
    geom_col(aes(fill = variable)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Sample names", y = "Total variants", title = paste0("Variant amino acid mutation statistics"))
  ggsave(file.path(report_outpath, paste0('dimsum_stage_merge_report_aamutationcounts.png')), d, width=12, height=8)
  #Plot 4: AA mutation percentages
  plot_df <- melt(aa_mut_df_collapse_perc, id="pairname")
  d <- ggplot(plot_df, aes(pairname, value)) +
    geom_col(aes(fill = variable)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Sample names", y = "Percentage of variants", title = paste0("Variant amino acid mutation statistics (percentage)"))
  ggsave(file.path(report_outpath, paste0('dimsum_stage_merge_report_aamutationpercentages.png')), d, width=12, height=8)
  #Plot 5: Nucleotide mutation counts
  nuc_mut_df <- data.frame(
    'nuc_mut_0'=sapply(dimsum_meta[['nuc_mutation_counts']], '[', '0'),
    'nuc_mut_1'=sapply(dimsum_meta[['nuc_mutation_counts']], '[', '1'),
    'nuc_mut_2'=sapply(dimsum_meta[['nuc_mutation_counts']], '[', '2'),
    'nuc_mut_sum'=sapply(dimsum_meta[['nuc_mutation_counts']], sum))
  nuc_mut_df[is.na(nuc_mut_df)] <- 0
  nuc_mut_df$pairname <- sapply(strsplit(merge_df$aligned_pair, '.split'), '[', 1)
  nuc_mut_df_collapse <- ddply(nuc_mut_df, "pairname", summarise, 
    nuc_mut_0 = sum(nuc_mut_0), 
    nuc_mut_1 = sum(nuc_mut_1), 
    nuc_mut_2 = sum(nuc_mut_2),
    nuc_mut_sum = sum(nuc_mut_sum))
  nuc_mut_df_collapse$nuc_mut_3plus = nuc_mut_df_collapse$nuc_mut_sum-nuc_mut_df_collapse$nuc_mut_0-nuc_mut_df_collapse$nuc_mut_1-nuc_mut_df_collapse$nuc_mut_2
  nuc_mut_df_collapse_perc = nuc_mut_df_collapse
  nuc_mut_df_collapse_perc[,c(2:5, 6)] <- nuc_mut_df_collapse_perc[,c(2:5, 6)]/nuc_mut_df_collapse_perc$nuc_mut_sum*100
  nuc_mut_df_collapse_perc <- nuc_mut_df_collapse_perc[,c(1:4, 6)]
  nuc_mut_df_collapse <- nuc_mut_df_collapse[,c(1:4, 6)]
  #Plot
  plot_df <- melt(nuc_mut_df_collapse, id="pairname")
  d <- ggplot(plot_df, aes(pairname, value)) +
    geom_col(aes(fill = variable)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Sample names", y = "Total variants", title = paste0("Variant nucleotide mutation statistics"))
  ggsave(file.path(report_outpath, paste0('dimsum_stage_merge_report_nucmutationcounts.png')), d, width=12, height=8)
  #Plot 6: Nucleotide mutation percentages
  plot_df <- melt(nuc_mut_df_collapse_perc, id="pairname")
  d <- ggplot(plot_df, aes(pairname, value)) +
    geom_col(aes(fill = variable)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Sample names", y = "Percentage of variants", title = paste0("Variant nucleotide mutation statistics (percentage)"))
  ggsave(file.path(report_outpath, paste0('dimsum_stage_merge_report_nucmutationpercentages.png')), d, width=12, height=8)
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  return(dimsum_meta_new)
}

#dimsum_stage_merge
#
# Merge all variant files.
#
# dimsum_meta: an experiment metadata object (required)
# merge_outpath: merged variant data output path (required)
# execute: whether or not to execute the system command (default: TRUE)
#
# Returns: an updated experiment metadata object.
#
dimsum_stage_merge <- function(
  dimsum_meta,
  merge_outpath,
  execute = TRUE,
  report = TRUE,
  report_outpath = NULL
  ){
  #Create cutadapt directory (if doesn't already exist)
  merge_outpath <- gsub("/$", "", merge_outpath)
  dir.create(merge_outpath)
  #WT nucleotide sequence
  wt_NTseq <- tolower(dimsum_meta[['wildtypeSequence']])
  #WT AA sequence
  wt_AAseq <- paste0(seqinr::translate(strsplit(wt_NTseq,split="")[[1]]),collapse="")
  #Initialise dictionary of AA mappings
  aa_dict <- list()
  #Initialise merged variant data object
  variant_data <- NULL
  #AA mutation filter dictionary
  aa_mutation_filter_dict <- list()
  #Nucleotide mutation count dictionary
  nuc_mutation_dict <- list()
  #AA mutation count dictionary
  aa_mutation_dict <- list()
  
  #Load filtered variant count files
  message("Loading filtered variant count files:")
  all_count <- file.path(dimsum_meta[['exp_design']]$aligned_pair_unique_tsv_directory, dimsum_meta[['exp_design']]$aligned_pair_unique_tsv)
  print(all_count)
  message("Processing...")
  for(count_file in all_count){
    print(count_file)
    #Check if this code should be executed
    if(execute){
      file_id <- gsub('.usearch.unique.tsv', '', basename(count_file))
      #Load file
      temp_file <- fread(count_file, header = T, sep="\t", stringsAsFactors = F)
      #Calculate number of aa mutations
      temp_file[,Nmut_aa := adist(aa_seq,wt_AAseq)]
      #Calculate number of nucleotide mutations
      temp_file[,Nmut_nt := adist(nt_seq,wt_NTseq)]
      #Number of sequences with greater than desired number of mutations
      aa_mutation_filter_dict[[count_file]] <- sum(temp_file[Nmut_aa > dimsum_meta[["maxAAMutations"]],]$count)
      #Save nucleotide mutation distribution
      nuc_mutation_dict[[count_file]] <- tapply(temp_file$count, temp_file$Nmut_nt, sum)
      #Save amino acid mutation distribution
      aa_mutation_dict[[count_file]] <- tapply(temp_file$count, temp_file$Nmut_aa, sum)
      #Subset to desired number of mutations
      temp_file <- temp_file[Nmut_aa <= dimsum_meta[["maxAAMutations"]],]
      #Save AA mapping to dictionary
      temp_dict <- as.list(temp_file$aa_seq)
      names(temp_dict) <- temp_file$nt_seq
      aa_dict <- c(aa_dict, temp_dict)
      aa_dict <- aa_dict[unique(names(aa_dict))]
      #First file loaded
      if(count_file == all_count[1]){
        variant_data <- temp_file
        names(variant_data)[grep(names(variant_data),pattern="count")] = paste0(file_id, "_count")
      #Not first file loaded (merge with previous data)
      }else{
        variant_data = merge(variant_data,temp_file[,.(nt_seq,.SD),,.SDcols=names(temp_file)[grep(names(temp_file),pattern="count")]],by="nt_seq",all = T)
        names(variant_data)[grep(names(variant_data),pattern=".SD")] = paste0(file_id, "_count")
      }
    }
  }
  #Check if this code should be executed
  if(execute){
    #Add AA sequence to merged variant data
    variant_data$aa_seq <- unlist(aa_dict[variant_data$nt_seq])
    #Calculate number of AA mutations of merged variant data
    variant_data[,Nmut_aa := adist(aa_seq,wt_AAseq)]
    #Replace NA counts with zeros
    variant_data[is.na(variant_data)] <- 0
    #Indicate WT sequence
    variant_data[nt_seq == wt_NTseq,WT := TRUE]
    #Indicate STOPs
    variant_data[,STOP := ifelse(length(grep(aa_seq,pattern="\\*"))==1,TRUE,FALSE),aa_seq]
    #Calculate number of nucleotide mutations of merged variant data
    variant_data[,Nmut_nt := adist(nt_seq,wt_NTseq)]
    #Merge split counts
    split_base <- unique(sapply(strsplit(colnames(variant_data)[grep('_count', colnames(variant_data))], "_split"), '[', 1))
    variant_data_merge <- sum_datatable_columns(dt=variant_data, column_patterns=split_base, suffix="_count")
    #Merge technical counts
    split_base <- unique(sapply(strsplit(colnames(variant_data_merge)[grep('_count', colnames(variant_data_merge))], "_t"), '[', 1))
    variant_data_merge <- sum_datatable_columns(dt=variant_data_merge, column_patterns=split_base, suffix="_count")
    #Save merged variant data
    message("Saving merged variant data...")
    save(variant_data_merge, file=file.path(merge_outpath, paste0(dimsum_meta[["project_name"]], '_variant_data_merge.RData')))
    message("Done")
  }
  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Merged variant data path
  dimsum_meta_new[["variant_data_merge_path"]] <- file.path(merge_outpath, paste0(dimsum_meta[["project_name"]], '_variant_data_merge.RData')) 
  #AA mutation results
  dimsum_meta_new[['exp_design']]$too_many_aa_mutations <- unlist(aa_mutation_filter_dict)  
  dimsum_meta_new[["aa_mutation_counts"]] <- aa_mutation_dict
  #Nucleotide mutation results
  dimsum_meta_new[["nuc_mutation_counts"]] <- nuc_mutation_dict 
  #Generate merge report
  if(report){
    dimsum_meta_new_report <- dimsum_stage_merge_report(dimsum_meta_new, report_outpath)
    return(dimsum_meta_new_report)
  }else{
    return(dimsum_meta_new)
  }
}

###########################
### GLOBALS
###########################

#Metadata object
exp_metadata <- list()

# #TEMP: set arguments manually (FOSintra_may2016)
# arg_list <- list()
# arg_list$fastqFileDir <- "/users/blehner/sequencing_data/Guillaume_Diss_EMBL/FOSintra_may2016/"
# arg_list$fastqFileExtension <- ".txt"
# arg_list$gzipped <- TRUE
# arg_list$stranded <- FALSE
# arg_list$experimentDesignPath <- "/users/blehner/afaure/DMS/pipelinetest_20180604/experimentDesign_FOSintra.txt"
# arg_list$cutadapt5First <- "AACCGGAGGAGGGAGCTG"
# arg_list$cutadapt5Second <- "GCTGCCAGGATGAACTC"
# arg_list$cutadapt3First <- "GAGTTCATCCTGGCAGC"
# arg_list$cutadapt3Second <- "CAGCTCCCTCCTCCGGTT"
# arg_list$cutadaptMinLength <- 50
# arg_list$cutadaptErrorRate <- 0.2
# arg_list$cutadaptDiscardUntrimmed <- F
# arg_list$usearchMinQual <- 30
# arg_list$usearchMaxee <- 0.5
# arg_list$usearchMinovlen <- 16
# arg_list$outputPath <- "/users/blehner/afaure/DMS/pipelinetest_20180604"
# arg_list$projectName <- "GD_FOSintra_2016-06-17"
# arg_list$wildtypeSequence <- "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA"
# arg_list$maxAAMutations <- 2
# arg_list$startStage <- 10
# arg_list$stopStage <- 0
# arg_list$numCores <- 10

# #TEMP: set arguments manually (XL_cI_low)
# arg_list <- list()
# arg_list$fastqFileDir <- "/users/blehner/sequencing_data/Xianghua_Li/cI_low_expression/"
# arg_list$fastqFileExtension <- ".fastq"
# arg_list$gzipped <- FALSE
# arg_list$stranded <- FALSE
# arg_list$barcodeDesignPath <- "/users/blehner/afaure/DMS/pipelinetest_20180603/barcodeDesign_cI_low_forcutadapt.txt"
# arg_list$barcodeErrorRate <- 0.25
# arg_list$experimentDesignPath <- "/users/blehner/afaure/DMS/pipelinetest_20180603/experimentDesign_cI_low.txt"
# arg_list$cutadapt5First <- "GCTTGAGGACGCACGTCGC"
# arg_list$cutadapt5Second <- "TCTGGCGATTGAAGGGCT"
# arg_list$cutadapt3First <- "AGCCCTTCAATCGCCAGA"
# arg_list$cutadapt3Second <- "GCGACGTGCGTCCTCAAGC"
# arg_list$cutadaptMinLength <- 50
# arg_list$cutadaptErrorRate <- 0.2
# arg_list$cutadaptDiscardUntrimmed <- T
# arg_list$usearchMinQual <- 20
# arg_list$usearchMaxee <- 0.5
# arg_list$usearchMinovlen <- 16
# arg_list$outputPath <- "/users/blehner/afaure/DMS/pipelinetest_20180603"
# arg_list$projectName <- "XL_cI_low_2015-09-15"
# arg_list$wildtypeSequence <- "CTTAAAGCAATTTATGAAAAAAAGAAAAATGAACTTGGCTTATCCCAGGAATCTGTCGCAGACAAGATGGGGATGGGGCAGTCAGGCGTTGGTGCTTTATTTAATGGCATCAATGCATTAAATGCTTATAACGCCGCATTGCTTGCAAAAATTCTCAAAGTTAGCGTTGAAGAATTT"
# arg_list$maxAAMutations <- 2
# arg_list$startStage <- 10
# arg_list$stopStage <- 0
# arg_list$numCores <- 10

# #TEMP: set arguments manually (XL_cI_high)
# arg_list <- list()
# arg_list$fastqFileDir <- "/users/blehner/sequencing_data/Xianghua_Li/cI_high_expression_2017-01-25-CA1FHANXX/"
# arg_list$fastqFileExtension <- ".fastq"
# arg_list$gzipped <- FALSE
# arg_list$stranded <- FALSE
# arg_list$barcodeDesignPath <- "/users/blehner/afaure/DMS/pipelinetest_20180523/barcodeDesign_cI_high.txt"
# arg_list$barcodeErrorRate <- 0.25
# arg_list$experimentDesignPath <- "/users/blehner/afaure/DMS/pipelinetest_20180523/experimentDesign_cI_high.txt"
# arg_list$cutadapt5First <- "GCTTGAGGACGCACGTCGC"
# arg_list$cutadapt5Second <- "TCTGGCGATTGAAGGGCT"
# arg_list$cutadapt3First <- "AGCCCTTCAATCGCCAGA"
# arg_list$cutadapt3Second <- "GCGACGTGCGTCCTCAAGC"
# arg_list$cutadaptMinLength <- 50
# arg_list$cutadaptErrorRate <- 0.2
# arg_list$cutadaptDiscardUntrimmed <- T
# arg_list$usearchMinQual <- 25
# arg_list$usearchMaxee <- 0.5
# arg_list$usearchMinovlen <- 16
# arg_list$outputPath <- "/users/blehner/afaure/DMS/pipelinetest_20180523"
# arg_list$projectName <- "XL_cI_high_2017-01-26"
# arg_list$wildtypeSequence <- "CTTAAAGCAATTTATGAAAAAAAGAAAAATGAACTTGGCTTATCCCAGGAATCTGTCGCAGACAAGATGGGGATGGGGCAGTCAGGCGTTGGTGCTTTATTTAATGGCATCAATGCATTAAATGCTTATAACGCCGCATTGCTTGCAAAAATTCTCAAAGTTAGCGTTGAAGAATTT"
# arg_list$maxAAMutations <- 2
# arg_list$startStage <- 10
# arg_list$stopStage <- 0
# arg_list$numCores <- 10

# #TEMP: set arguments manually (Li2016_101)
# arg_list <- list()
# arg_list$fastqFileDir <- "/users/blehner/jschmiedel/DMS2struct/datasets/tRNA_Li2016/FastQ/"
# arg_list$fastqFileExtension <- ".fastq"
# arg_list$gzipped <- FALSE
# arg_list$stranded <- TRUE
# arg_list$barcodeDesignPath <- NULL
# arg_list$experimentDesignPath <- "/users/blehner/afaure/DMS/pipelinetest_20180528/experimentDesign_tRNA_Li2016_101.txt"
# arg_list$cutadapt5First <- "AGTTCAACCAAGTTG...TTGATTATTTTTTTTT"
# arg_list$cutadapt5Second <- "AAAAAAAAATAATCAA...CAACTTGGTTGAACT"
# arg_list$cutadapt3First <- NULL
# arg_list$cutadapt3Second <- NULL
# arg_list$cutadaptMinLength <- 50
# arg_list$cutadaptErrorRate <- 0.2
# arg_list$cutadaptDiscardUntrimmed <- F
# arg_list$usearchMinQual <- 30
# arg_list$usearchMaxee <- 0.5
# arg_list$usearchMinovlen <- 16
# arg_list$outputPath <- "/users/blehner/afaure/DMS/pipelinetest_20180528"
# arg_list$projectName <- "tRNA_Li2016_101"
# arg_list$startStage <- 10
# arg_list$stopStage <- 0
# wildtypeSequence <- "GTTCCGTTGGCGTAATGGTAACGCGTCTCCCTCCTAAGGAGAAGACTGCGGGTTCGAGTCCCGTACGGAACG"
# maxAAMutations <- 2
# numCores <- 10

# #TEMP: set arguments manually (Li2016_126)
# arg_list <- list()
# arg_list$fastqFileDir <- "/users/blehner/jschmiedel/DMS2struct/datasets/tRNA_Li2016/FastQ/"
# arg_list$fastqFileExtension <- ".fastq"
# arg_list$gzipped <- FALSE
# arg_list$stranded <- TRUE
# arg_list$barcodeDesignPath <- NULL
# arg_list$experimentDesignPath <- "/users/blehner/afaure/DMS/pipelinetest_20180528/experimentDesign_tRNA_Li2016_126.txt"
# arg_list$cutadapt5First <- "AGTTCAACCAAGTTG...TTGATTATTTTTTTTT"
# arg_list$cutadapt5Second <- "AAAAAAAAATAATCAA...CAACTTGGTTGAACT"
# arg_list$cutadapt3First <- NULL
# arg_list$cutadapt3Second <- NULL
# arg_list$cutadaptMinLength <- 50
# arg_list$cutadaptErrorRate <- 0.2
# arg_list$cutadaptDiscardUntrimmed <- F
# arg_list$usearchMinQual <- 30
# arg_list$usearchMaxee <- 0.5
# arg_list$usearchMinovlen <- 16
# arg_list$outputPath <- "/users/blehner/afaure/DMS/pipelinetest_20180528"
# arg_list$projectName <- "tRNA_Li2016_126"
# arg_list$startStage <- 10
# arg_list$stopStage <- 0
# wildtypeSequence <- "GTTCCGTTGGCGTAATGGTAACGCGTCTCCCTCCTAAGGAGAAGACTGCGGGTTCGAGTCCCGTACGGAACG"
# maxAAMutations <- 2
# numCores <- 10

# #TEMP: set arguments manually (Li2018)
# arg_list <- list()
# arg_list$fastqFileDir <- "/users/blehner/jschmiedel/DMS2struct/datasets/tRNA_Li2018/FastQ/"
# arg_list$fastqFileExtension <- ".fastq"
# arg_list$gzipped <- FALSE
# arg_list$stranded <- TRUE
# arg_list$barcodeDesignPath <- NULL
# arg_list$experimentDesignPath <- "/users/blehner/afaure/DMS/pipelinetest_20180528/experimentDesign_tRNA_Li2018.txt"
# arg_list$cutadapt5First <- "AGTTCAACCAAGTTG...TTGATTATTTTTTTTT"
# arg_list$cutadapt5Second <- "AAAAAAAAATAATCAA...CAACTTGGTTGAACT"
# arg_list$cutadapt3First <- NULL
# arg_list$cutadapt3Second <- NULL
# arg_list$cutadaptMinLength <- 50
# arg_list$cutadaptErrorRate <- 0.2
# arg_list$cutadaptDiscardUntrimmed <- F
# arg_list$usearchMinQual <- 30
# arg_list$usearchMaxee <- 0.5
# arg_list$usearchMinovlen <- 16
# arg_list$outputPath <- "/users/blehner/afaure/DMS/pipelinetest_20180528"
# arg_list$projectName <- "tRNA_Li2018"
# arg_list$startStage <- 10
# arg_list$stopStage <- 0
# wildtypeSequence <- "GTTCCGTTGGCGTAATGGTAACGCGTCTCCCTCCTAAGGAGAAGACTGCGGGTTCGAGTCCCGTACGGAACG"
# maxAAMutations <- 2
# numCores <- 10

### Save metadata
#Remove trailing "/" if present
exp_metadata[["fastq_path_original"]] <- gsub("/$", "", arg_list$fastqFileDir)
exp_metadata[["fastq_file_extension"]] <- arg_list$fastqFileExtension
exp_metadata[["gzipped"]] <- arg_list$gzipped
exp_metadata[["stranded"]] <- arg_list$stranded
exp_metadata[["barcode_design_path"]] <- arg_list$barcodeDesignPath
exp_metadata[["barcodeErrorRate"]] <- arg_list$barcodeErrorRate
exp_metadata[["experiment_design_path"]] <- arg_list$experimentDesignPath
exp_metadata[["cutadaptCut5First"]] <- arg_list$cutadaptCut5First
exp_metadata[["cutadaptCut5Second"]] <- arg_list$cutadaptCut5Second
exp_metadata[["cutadaptCut3First"]] <- arg_list$cutadaptCut3First
exp_metadata[["cutadaptCut3Second"]] <- arg_list$cutadaptCut3Second
exp_metadata[["cutadapt5First"]] <- arg_list$cutadapt5First
exp_metadata[["cutadapt5Second"]] <- arg_list$cutadapt5Second
exp_metadata[["cutadapt3First"]] <- arg_list$cutadapt3First
exp_metadata[["cutadapt3Second"]] <- arg_list$cutadapt3Second
exp_metadata[["cutadaptMinLength"]] <- arg_list$cutadaptMinLength
exp_metadata[["cutadaptErrorRate"]] <- arg_list$cutadaptErrorRate
exp_metadata[["cutadaptDiscardUntrimmed"]] <- arg_list$cutadaptDiscardUntrimmed
exp_metadata[["usearchMinQual"]] <- arg_list$usearchMinQual
exp_metadata[["usearchMaxee"]] <- arg_list$usearchMaxee
exp_metadata[["usearchMinovlen"]] <- arg_list$usearchMinovlen
#Remove trailing "/" if present
exp_metadata[["output_path"]] <- gsub("/$", "", arg_list$outputPath)
exp_metadata[["project_name"]] <- arg_list$projectName
exp_metadata[["wildtypeSequence"]] <- arg_list$wildtypeSequence
exp_metadata[["maxAAMutations"]] <- arg_list$maxAAMutations

#First pipeline stage to run
first_stage <- arg_list$startStage
last_stage <- arg_list$stopStage
#Number of cores
num_cores <- arg_list$numCores

###########################
### MAIN
###########################

### Output file path, working and temp directories
###########################

#Create working directory (if doesn't already exist)
exp_metadata[["project_path"]] <- file.path(exp_metadata[["output_path"]], exp_metadata[["project_name"]])
dir.create(exp_metadata[["project_path"]])
#Set working directory
setwd(exp_metadata[["project_path"]])
#Create temp directory (if doesn't already exist)
exp_metadata[["tmp_path"]] <- file.path(exp_metadata[["project_path"]], "tmp")
dir.create(exp_metadata[["tmp_path"]])

### Get experiment design
###########################

#TODO: check if file exists
#TODO: check if all fastq files exist
exp_metadata[["exp_design"]] <- read.table(exp_metadata[["experiment_design_path"]], header = T, stringsAsFactors = F, sep="\t")
#Add original FASTQ directory
exp_metadata[["exp_design"]]$pair_directory <- exp_metadata[["fastq_path_original"]]

### Get barcode design (if provided)
###########################

if(!is.null(exp_metadata[["barcode_design_path"]])){
  exp_metadata[["barcode_design"]] <- read.table(exp_metadata[["barcode_design_path"]], header = T, stringsAsFactors = F, sep="\t")
}

### Pipeline stages
###########################

### Step 0: Start pipeline tracking
pipeline <- list()
pipeline[['0_original']] <- exp_metadata

### Step 1: Run demultiplex on all fastq files
pipeline[['1_demultiplex']] <- dimsum_stage_demultiplex(pipeline[['0_original']], file.path(pipeline[['0_original']][["tmp_path"]], "demultiplex"), 
  execute = (first_stage <= 1 & (last_stage == 0 | last_stage >= 1)))

### Step 2: Run FASTQC on all fastq files
pipeline[['2_fastqc']] <- dimsum_stage_fastqc(pipeline[['1_demultiplex']], file.path(pipeline[['1_demultiplex']][["tmp_path"]], "fastqc"), 
  execute = (first_stage <= 2 & (last_stage == 0 | last_stage >= 2)), report_outpath = file.path(pipeline[['1_demultiplex']][["tmp_path"]], "reports"))

### Step 3: Unzip FASTQ files if necessary
pipeline[['3_fastq']] <- dimsum_stage_unzip(pipeline[['2_fastqc']], file.path(pipeline[['2_fastqc']][["tmp_path"]], "fastq"), 
  execute = (first_stage <= 3 & (last_stage == 0 | last_stage >= 3)))

### Step 4: Split FASTQ files
pipeline[['4_split']] <- dimsum_stage_split(pipeline[['3_fastq']], file.path(pipeline[['3_fastq']][["tmp_path"]], "split"), 
  execute = (first_stage <= 4 & (last_stage == 0 | last_stage >= 4)))

### Step 5: Remove adapters from FASTQ files with cutadapt if necessary
pipeline[['5_cutadapt']] <- dimsum_stage_cutadapt(pipeline[['4_split']], file.path(pipeline[['4_split']][["tmp_path"]], "cutadapt"), 
  execute = (first_stage <= 5 & (last_stage == 0 | last_stage >= 5)), report_outpath = file.path(pipeline[['4_split']][["tmp_path"]], "reports"))

### Step 6: Merge paired-end reads with USEARCH
pipeline[['6_usearch']] <- dimsum_stage_usearch(pipeline[['5_cutadapt']], file.path(pipeline[['5_cutadapt']][["tmp_path"]], "usearch"), 
  execute = (first_stage <= 6 & (last_stage == 0 | last_stage >= 6)), report_outpath = file.path(pipeline[['5_cutadapt']][["tmp_path"]], "reports"))

### Step 7: Get unique aligned read counts with FASTX-Toolkit
pipeline[['7_unique']] <- dimsum_stage_unique(pipeline[['6_usearch']], file.path(pipeline[['6_usearch']][["tmp_path"]], "usearch"), 
  execute = (first_stage <= 7 & (last_stage == 0 | last_stage >= 7)))

### Step 8: Construct filtered variant count table with corresponding amino acid sequences
pipeline[['8_filter']] <- dimsum_stage_filter(pipeline[['7_unique']], file.path(pipeline[['7_unique']][["tmp_path"]], "usearch"), 
  execute = (first_stage <= 8 & (last_stage == 0 | last_stage >= 8)))

### Step 9: Merge variant count tables
pipeline[['9_merge']] <- dimsum_stage_merge(pipeline[['8_filter']], file.path(pipeline[['8_filter']][["tmp_path"]], "merge"), 
  execute = (first_stage <= 9 & (last_stage == 0 | last_stage >= 9)), report_outpath = file.path(pipeline[['8_filter']][["tmp_path"]], "reports"))

### Save workspace
###########################

message("Saving workspace image...")
save.image(file=file.path(pipeline[['9_merge']][["tmp_path"]], paste0(pipeline[['9_merge']][["project_name"]], '_workspace.RData')))
message("Done")


