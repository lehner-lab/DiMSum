
#' dimsum__fastqc_report
#'
#' Generate FASTQC summary plots for all fastq files.
#'
#' @param dimsum_meta an experiment metadata object (required)
#' @param report_outpath FASTQC report output path (required)
#'
#' @return an updated experiment metadata object
#' @export
dimsum__fastqc_report <- function(
  dimsum_meta,
  report_outpath
  ){

  #Create report directory (if doesn't already exist)
  report_outpath <- gsub("/$", "", report_outpath)
  suppressWarnings(dir.create(report_outpath))

  #Input files
  all_fastqc <- file.path(dimsum_meta[["exp_design"]][,'fastqc_directory'], unique(c(dimsum_meta[['exp_design']][,"pair1_fastqc"], dimsum_meta[['exp_design']][,"pair2_fastqc"])))
  #Check if all input files exist
  dimsum__check_files_exist(
    required_files = all_fastqc,
    execute = TRUE,
    exit = FALSE)

  #Initialise read length
  dimsum_meta[['exp_design']][,"pair1_length"] <- NA
  dimsum_meta[['exp_design']][,"pair2_length"] <- NA

  #Get results for all fastq files
  for(col_name in c('pair1_fastqc', 'pair2_fastqc')){
    fastqc_files <- file.path(dimsum_meta[['exp_design']][,'fastqc_directory'], dimsum_meta[['exp_design']][,col_name])
    fastqc_list <- list()
    encoding <- list()
    for(i in 1:length(fastqc_files)){
      temp_out <- readLines(fastqc_files[i], 500)
      filename <- gsub('Filename\\t', '', temp_out[4])
      encoding[[i]] <- gsub('Encoding\\t', '', temp_out[6])
      read_length <- gsub('Sequence length\\t', '', temp_out[9])
      if(grepl("-", read_length)){
        dimsum_meta[['exp_design']][i,gsub("_fastqc", "_length", col_name)] <- as.numeric(unlist(strsplit(read_length, "-"))[2])
      }else{
        dimsum_meta[['exp_design']][i,gsub("_fastqc", "_length", col_name)] <- as.numeric(read_length)
      }
      temp_nlines <- grep('>>END_MODULE', temp_out)[2]
      temp_out <- temp_out[c(13:(temp_nlines-1))]
      temp_out_data <- strsplit(temp_out[2:length(temp_out)], '\\t')
      temp_df <- as.data.frame(do.call('rbind', lapply(lapply(temp_out_data, '[', -1), as.numeric)))
      rownames(temp_df) <- sapply(temp_out_data, '[', 1)
      colnames(temp_df) <- unlist(strsplit(temp_out[1], '\\t'))[-1]
      fastqc_list[[filename]] <- temp_df
    }
    fastqc_df1 <- dimsum__cbind_fill(lapply(fastqc_list, '[', 'Mean'))
    #Temporarily code NAs as Inf (dimsum__cbind_fill assumes no NAs)
    fastqc_df1[is.na(fastqc_df1)] <- Inf
    fastqc_df2 <- dimsum__cbind_fill(lapply(fastqc_list, '[', '10th Percentile'))
    #Temporarily code NAs as Inf (dimsum__cbind_fill assumes no NAs)
    fastqc_df2[is.na(fastqc_df2)] <- Inf

    #convert to lists
    fastqc_df1_list <- lapply(lapply(as.list(fastqc_df1), as.data.frame), 
      function(x){
        rownames(x)<-rownames(fastqc_df1)
        return(x)})
    #convert to list
    fastqc_df2_list <- lapply(lapply(as.list(fastqc_df2), as.data.frame), 
      function(x){
        rownames(x)<-rownames(fastqc_df2)
        return(x)})

    #Join lists
    fastqc_df_list <- c(fastqc_df1_list, fastqc_df2_list)
    #Consistent ranges (dimsum__cbind_fill assumes no NAs)
    fastqc_df <- dimsum__cbind_fill(fastqc_df_list)
    #Decode Inf as NA
    fastqc_df[is.infinite(as.matrix(fastqc_df))] <- NA
    #Split
    fastqc_df1_plot <- data.frame(fastqc_df[,1:dim(fastqc_df1)[2]])
    rownames(fastqc_df1_plot) <- rownames(fastqc_df)
    colnames(fastqc_df1_plot) <- colnames(fastqc_df)[1:dim(fastqc_df1)[2]]
    fastqc_df1_plot[,'base_position'] <- 1:length(rownames(fastqc_df1_plot))
    fastqc_df2_plot <- data.frame(fastqc_df[,(dim(fastqc_df2)[2]+1):dim(fastqc_df)[2]])
    rownames(fastqc_df2_plot) <- rownames(fastqc_df)
    colnames(fastqc_df2_plot) <- colnames(fastqc_df)[(dim(fastqc_df2)[2]+1):dim(fastqc_df)[2]]    
    fastqc_df2_plot[,'base_position'] <- 1:length(rownames(fastqc_df2_plot))

    #Plot
    plot_df1 <- reshape2::melt(fastqc_df1_plot, id="base_position")
    plot_df1[,'statistic'] <- 'Mean'
    plot_df2 <- reshape2::melt(fastqc_df2_plot, id="base_position")
    plot_df2[,'statistic'] <- '10th Percentile'
    plot_df <- rbind(plot_df1, plot_df2)
    plot_df[,'Read_name'] <- factor(plot_df[,'variable'])
    #Remove NAs
    plot_df <- plot_df[!is.na(plot_df[,'value']),]
    # #Variable region boundaries
    # temp_adapt5 <- c("cutadapt5First", "cutadapt5Second")[as.numeric(gsub("pair|_fastqc", "", col_name))]
    # temp_cut5 <- c("cutadaptCut5First", "cutadaptCut5Second")[as.numeric(gsub("pair|_fastqc", "", col_name))]
    # #Positions of 5' and 3' boundaries of constant regions
    # vr_5 <- apply(cbind(nchar(dimsum_meta[['exp_design']][,temp_adapt5]), dimsum_meta[['exp_design']][,temp_cut5]), 1, sum, na.rm = T)
    # vr_3 <- vr_5 + nchar(dimsum_meta[["wildtypeSequence"]])
    # #Only show sequenced boundaries
    # vr_5 <- vr_5[vr_5<=dimsum_meta[['exp_design']][,gsub("_fastqc", "_length", col_name)]]
    # vr_3 <- vr_3[vr_3<=dimsum_meta[['exp_design']][,gsub("_fastqc", "_length", col_name)]]
    # #Min and max boundary positions
    # vr_5 <- range(vr_5)
    # if(length(vr_3)!=0){vr_3 <- range(vr_3)}
    # vr_boundaries <- unique(c(vr_5, vr_3))
    # #Plot axis position
    # pos_start <- as.numeric(sapply(strsplit(rownames(fastqc_df1_plot), "-"), '[', 1))
    # vr_boundaries_pos <- NULL
    # for(b in vr_boundaries){
    #   vr_boundaries_pos <- c(vr_boundaries_pos, which(pos_start/b>=1)[1])
    # }
    #Encoding
    encoding_format <- paste(unique(unlist(encoding)), collapse = ", ")
    #Plot (if not pair2 or paired design)
    if(col_name != "pair2_fastqc" | dimsum_meta[["paired"]]){
      d <- ggplot2::ggplot(plot_df, ggplot2::aes(base_position, value, color = Read_name)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept=c(20, 28), linetype = 2) +
        ggplot2::theme_bw() +
        ggplot2::coord_cartesian(ylim = c(0, max(plot_df[,'value']))) + 
        # ggplot2::geom_vline(xintercept = vr_boundaries_pos, linetype = 2) +
        #ggplot2::annotate("text", label = "variable region" , x = median(unique(plot_df[,"base_position"])), y = 0) + 
        ggplot2::scale_x_continuous(
        breaks = (1:length(rownames(fastqc_df1_plot)))[seq(1, length(rownames(fastqc_df1_plot)), 5)],
        label = rownames(fastqc_df1_plot)[seq(1, length(rownames(fastqc_df1_plot)), 5)]) +
        #ggplot2::labs(x = "Position in read (bp)", y = "Quality score", title = paste0("Read ", gsub("pair|_fastqc", "", col_name), " quality scores across all bases (", encoding_format, ")"))
        ggplot2::labs(x = "Position in read (bp)", y = "Quality score", title = paste0("Base quality score encoding: ", encoding_format))
      d <- d + ggplot2::facet_wrap(~statistic, nrow=2, ncol=1)
      dimsum__save_png(file.path(report_outpath, paste0('dimsum__fastqc_report_', col_name, '.png')), d, width=12, height=8)
    }
  }

  #Render report
  dimsum__render_report(dimsum_meta = dimsum_meta)

  #New experiment metadata
  dimsum_meta_new <- dimsum_meta
  #Update fastq metadata
  return(dimsum_meta_new)
}

