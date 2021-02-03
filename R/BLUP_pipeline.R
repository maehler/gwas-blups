suppressPackageStartupMessages({
  library(dplyr)
  library(here)
})
source(here("R/generate_BLUP.R"))

#######################################################################
# Set arguments 
#######################################################################

input_data <- "/mnt/picea/projects/aspseq/nstreet/swasp/Sara/GWAS_buds/Metabolite_data/raw/s12.Buds.LC.targeted.conc.normalised.Jan2018_Updated_Aug2020.txt"
out_dir <- "/mnt/picea/projects/aspseq/nstreet/swasp/Sara/GWAS_buds/Metabolite_data/BLUPs/BLUPs_SPGs_pipeline_test"
data_name <- "SPGs" # File name prefix
line_name <- "Genotype" # Has to match line column name in input_data. E.g., "Genotype", "Clone". Required for boxplots. 

# BLUP arguments - specify column number. See generate_BLUP.R function for more info.    
random_effect = 4
sample_column = 1
fixed_effect = 7
start_column = 16
col_convert = TRUE

#######################################################################
# Generate BLUPs
#######################################################################

# Create stats folder
stats_out_dir <- file.path(out_dir, "stats")
dir.create(stats_out_dir, recursive = TRUE, showWarnings = FALSE)

# Read file 
raw_data <- readr::read_tsv(input_data)

# Print arguments 
message("Arguments provided")
message(paste("\tInput data:", input_data))
message(paste("\tOutput directory:", out_dir))
message(paste("\tRandom effect:", colnames(raw_data)[random_effect]))
message(paste("\tFixed effect:", colnames(raw_data)[fixed_effect]))
message(paste("\tTrait start column:", colnames(raw_data)[start_column]))

# Generate BLUP
results <-
  generate_BLUP(dat = raw_data, sample_column = sample_column, random_effect = random_effect, start_column = start_column, fixed_effect = fixed_effect, col_convert = col_convert)

if (is.list(results)) {
  readr::write_tsv(results$BLUP, file = file.path(out_dir, paste0(data_name, "_BLUP.txt")))
  saveRDS(results$Outliers_residuals, file = file.path(stats_out_dir, paste0(data_name, "_Outliers_residuals.rds")))
  writeLines(paste0(results$Not_converge_columns, collapse = ", "), con = file.path(stats_out_dir, paste0(data_name, "_Traits_not_converge.txt")))
  readr::write_tsv(results$Outlier_data, file = file.path(stats_out_dir, paste0(data_name, "_Outlier_data.txt")))
  readr::write_tsv(results$Outlier_removed_data, file = file.path(stats_out_dir, paste0(data_name, "_Outlier_removed_data.txt")))
  readr::write_tsv(results$Shapiro_raw_normtest, file = file.path(stats_out_dir, paste0(data_name, "_raw_shapiro.txt")))
  readr::write_tsv(results$BLUP_shapiro, file = file.path(stats_out_dir, paste0(data_name, "_BLUP_shapiro.txt")))
}

#######################################################################
# Generate plots
#######################################################################

# Histogram BLUPs
if (is.list(results)) {
  BLUP <- results$BLUP
  plot_lst <- lapply(names(BLUP[-1]), function(t) {
    hist_plot <- ggplot2::ggplot(BLUP, ggplot2::aes(.data[[t]])) + 
      ggplot2::geom_histogram(bins = 20) +
      ggplot2::ggtitle(t)
  })
  message("Creating BLUP histogram plots")
  ml <- gridExtra::marrangeGrob(plot_lst, nrow = 6, ncol = 6)
  ggplot2::ggsave(file.path(stats_out_dir, "BLUP_histogram.pdf"), ml, width = 20, height = 15)
}

# Non-transformed qqplot BLUPs
if (is.list(results)) {
  BLUP <- results$BLUP
  trans <- results$Shapiro_raw_normtest
  trans_TRUE <- trans %>% subset(Shapiro_OK == TRUE) %>% select(c(Trait)) %>% pull()
  BLUP_trans_TRUE <- BLUP %>% dplyr::select(any_of(trans_TRUE))

  plot_lst <- lapply(names(BLUP_trans_TRUE), function(t) {
    q_plot <- ggplot2::ggplot(BLUP_trans_TRUE, ggplot2::aes(sample = .data[[t]])) + 
      ggplot2::stat_qq() +
      ggplot2::stat_qq_line() + 
      ggplot2::labs(title = t)
  })
  message("Creating BLUP (non-transformed) qqplots")
  ml <- gridExtra::marrangeGrob(plot_lst, nrow = 6, ncol = 6)
  ggplot2::ggsave(file.path(stats_out_dir, "BLUP_qqplot_nontrans.pdf"), ml, width = 20, height = 15)
}

# Transformed qqplot BLUPs
if (is.list(results)) {
  BLUP <- results$BLUP
  trans <- results$Shapiro_raw_normtest
  trans_FALSE <- trans %>% subset(Shapiro_OK == FALSE) %>% select(c(Trait)) %>% pull()
  BLUP_trans_FALSE <- BLUP %>% dplyr::select(any_of(trans_FALSE))
  
  plot_lst <- lapply(names(BLUP_trans_FALSE), function(t) {
      q_plot <- ggplot2::ggplot(BLUP_trans_FALSE, ggplot2::aes(sample = .data[[t]])) + 
        ggplot2::stat_qq() +
        ggplot2::stat_qq_line() +
        ggplot2::labs(title = t)
  })
    message("Creating BLUP (transformed) qqplots")
    ml <- gridExtra::marrangeGrob(plot_lst, nrow = 6, ncol = 6)
    ggplot2::ggsave(file.path(stats_out_dir, "BLUP_qqplot_trans.pdf"), ml, width = 20, height = 15)
}

# Boxplots - color outliers
if (is.list(results)) {
  outlier_list <- results$Outliers_residuals
  trait_names <- outlier_list[lapply(outlier_list,length)>0] %>% names()
  raw_data$my_color <- "NULL"
  raw_data[[line_name]] <- as.factor(raw_data[[line_name]])
  
  cbPalette = c("darkcyan", "coral1")
  plot_lst <- lapply(trait_names, function(t) {
    my_rows <- outlier_list[[t]]
    for(i in 1:nrow(raw_data)){
      if(i %in% my_rows){
        raw_data[i, "my_color"] <- "rm_outliers"
      } else {
        raw_data[i, "my_color"] <- "kept_samples"
      }
    }
    p <- ggplot2::ggplot(raw_data, ggplot2::aes(x = .data[[line_name]], y = .data[[t]])) +
      ggplot2::geom_boxplot() + 
      ggplot2::geom_jitter(ggplot2::aes(color = `my_color`)) + 
      ggplot2::scale_fill_manual(values = cbPalette) +
      ggplot2::scale_color_manual(values = cbPalette) +
      ggplot2::labs(title = t) 
    return(p)
  })
  message("Creating boxplots")
  ml <- gridExtra::marrangeGrob(plot_lst, nrow = 1, ncol = 1)
  ggplot2::ggsave(file.path(stats_out_dir, "Boxplots_residual_outlier_rm.pdf"), ml, width = 25, height = 12)
}

#######################################################################
# Save session info
#######################################################################

writeLines(capture.output(sessionInfo()), con = file.path(stats_out_dir, paste("sessionInfo",Sys.Date() ,".txt", sep="")))
message("Done")