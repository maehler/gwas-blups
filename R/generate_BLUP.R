#' generate BLUP function
#'
#' @description The goal of generate_BLUP is to run Best Linear Unbiased Predictions. Based on modified code from Happi GWAS generate_BLUP function. 
#' @importFrom magrittr %>%
#' @importFrom data.table as.data.table
#' @param dat An input dataset.
#' @param random_effect The random effect column(s). The first column number in vector will be the selected BLUP intercept. Random effect will be included in the model as: (1 | random effect).
#' @param sample_column The sample ID column. Needed for orderNorm calculation. 
#' @param start_column The start column index for traits.
#' @param fixed_effect The fixed effect column(s).
#' @param col_convert If TRUE, converts model effect columns to factor type.
#' @return blue or NULL if something missing.
#' @keywords BLUP, Pre-GWAS
#' @export
#'
generate_BLUP <- function(dat, random_effect, sample_column, start_column, fixed_effect = NULL, col_convert = TRUE) {
  
  # Check if input arguments are given
  if (missing(dat)) {
    stop("Data needs to be provided")
  }
  
  if (missing(random_effect)) {
    stop("Random effect column(s) needs to be provided")
  } else if (!is.numeric(random_effect) && any(random_effect < 1)) {
    stop("Random effect needs to be > 0")
  }
  
  if (missing(sample_column)) {
    stop("Sample column needs to be provided")
  } else if (!is.numeric(sample_column) && sample_column < 1) {
    stop("Sample column needs to be > 0")
  }
  
  if (missing(start_column)) {
    stop("Start column needs to be provided")
  } else if (!is.numeric(start_column) && start_column < 1) {
    stop("Start column needs to be > 0")
  }
  
  if (!is.null(fixed_effect)) {
    if(!is.numeric(fixed_effect) && any(fixed_effect) < 1) {
      stop("Fixed effect column number(s) must be integers > 0")
    }
  }
  
  if (col_convert) {
    cat("\tConverting lmer model effect columns to factor type\n")
    dat <- dplyr::mutate_at(dat, random_effect, as.factor)
    if (!is.null(fixed_effect)) {
      dat <- dplyr::mutate_at(dat, fixed_effect, as.factor)
    }
  }
  
  # Convert the response columns to numeric
  dat <- dplyr::mutate_at(dat, seq(start_column, ncol(dat), as.numeric))
  
  # Extract fixed effect name(s) (if provided)
  if (length(fixed_effect) > 0) {
  fixed_ef <- colnames(dat)[fixed_effect]
  }
  
  # Use this sanity check to see if column names are consistent throughout the analysis.
  start_colnames <- colnames(dat)
  
  #######################################################################
  ## Outlier Removal
  #######################################################################
  print("Remove outliers")
  
  # Create lmer formula
  if (length(random_effect) > 0 & length(fixed_effect) > 0) { # with fixed effect
    termlabels <- c(fixed_ef)
    for (i in 1:length(random_effect)) {
      my_col <- random_effect[i]
      temp <- paste("(1|", colnames(dat)[my_col], ")", sep = "")
      termlabels <- c(termlabels, temp)
    }
  } else if (length(random_effect) > 0 & is.null(fixed_effect)){ # no fixed effect
    termlabels <- c()
    for (i in 1:length(random_effect)) {
      my_col <- random_effect[i]
      temp <- paste("(1|", colnames(dat)[my_col], ")", sep = "")
      termlabels <- c(termlabels, temp)
    }
  }
  
  cat(paste("\tlmer model, effect included in formula:",termlabels,"\n"))
  
  # Find outliers
  trait_names <- colnames(dat)[start_column:ncol(dat)]
  outliers_residuals <- lapply(trait_names, function(t){
    res <- car::outlierTest(lme4::lmer(formula = reformulate(termlabels = termlabels, response = t), data = dat, REML = TRUE), n.max
                            = nrow(dat), cutoff = 0.05)
    rm_row <- res$bonf.p %>% rmutil::as.data.frame() %>% dplyr::filter(. < 0.05) %>% row.names()
    return(rm_row)
  })
  names(outliers_residuals) <- trait_names
  
  if(!identical(outliers_residuals, character(0))){
    temp_outliers_residuals <- outliers_residuals
    
    outlier_dat <- data.frame()
    
    # Remove outliers 
    for (i in 1:length(temp_outliers_residuals)) {
      rows <- match(temp_outliers_residuals[[i]], row.names(dat[names(temp_outliers_residuals)[i]]))
      columns <- which(grepl(names(temp_outliers_residuals)[i], colnames(dat)))
      temp_outliers_residuals[[i]] <- dat[rows, c(sample_column, columns)]

      dat[rows, columns] <- NA
      
      if(nrow(temp_outliers_residuals[[i]]) > 0){
        if(nrow(outlier_dat) == 0){
          outlier_dat <- temp_outliers_residuals[[i]]
        }else{
          outlier_dat <- merge(outlier_dat, temp_outliers_residuals[[i]],
                               by = intersect(colnames(outlier_dat), colnames(temp_outliers_residuals[[i]])),
                               all = TRUE)
        }
      }
    }
    # Re-arrange row names
    row.names(dat) <- seq(from = 1, to = nrow(dat), by = 1)
  }

  outlier_removed_dat = dat
  
  # Sanity check 
  if(!identical(start_colnames, colnames(dat))){
    stop(paste("Column names are changed. The column assumption no longer hold."))
  }
  
  #######################################################################
  # orderNorm
  #######################################################################
  print("Calculate orderNorm transformed")
  
  # Extract trait names 
  trait_names <- colnames(dat)[c(start_column:ncol(dat))]
    
  # Store shapiro test 
  shap_norm <- lapply(trait_names, function(t) {
    shap_res <- shapiro.test(dat[[t]]) 
    shap_p <- shap_res$p.value
    data.frame(Trait = t,
               Shapiro_pval = shap_p,
               Shapiro_OK = ifelse(shap_p > 0.05, "TRUE", "FALSE"),
               stringsAsFactors = FALSE)
  }) %>% bind_rows %>% as.data.frame()
  
  # Extract names of trait in need of transfromation
  trans_trait <- shap_norm %>% subset(Shapiro_OK == FALSE) %>% dplyr::select(c(Trait)) %>% dplyr::pull()
  
  # Transform those traits
  ## transform 
  for(t in trans_trait){
    orderNorm_obj <- bestNormalize::orderNorm(dat[[t]])
    dat[[t]] <- orderNorm_obj$x.t
  }
  
  ## catch warning messages 
   for(t in trans_trait){
     tryCatch({
       orderNorm_obj <- bestNormalize::orderNorm(dat[[t]])
     }, warning = function(w){
       cat(paste(t, ":", w))
     })
   }
  
  # Sanity check 
  if(!identical(start_colnames, colnames(dat))){
    stop(paste("Column names are changed. The column assumption no longer hold."))
  }
  
  #######################################################################
  ## generate BLUP
  #######################################################################
  print("Calculate BLUPs")
  
  # Create lmer formula
  if (length(random_effect) > 0 & length(fixed_effect) > 0) { # with fixed effect(s)
    termlabels <- c(fixed_ef)
    for (i in 1:length(random_effect)) {
      my_col <- random_effect[i]
      temp <- paste("(1|", colnames(dat)[my_col], ")", sep = "")
      termlabels <- c(termlabels, temp)
    }
  } else if (length(random_effect) > 0 & is.null(fixed_effect)){ # no fixed effect(s)
    termlabels <- c()
    for (i in 1:length(random_effect)) {
      my_col <- random_effect[i]
      temp <- paste("(1|", colnames(dat)[my_col], ")", sep = "")
      termlabels <- c(termlabels, temp)
    }
  }
  cat(paste("\tlmer model, effect included in formula:",termlabels,"\n"))
  
  # fit the model
  BLUP_out <- list()
  for(i in start_column:ncol(dat)){
    
    dat[is.infinite(dat[,i]),i] = NA
    dat[is.nan(dat[,i]),i] = NA
    
    lme <- lme4::lmer(formula = reformulate(termlabels = termlabels, response = colnames(dat)[i]), data = dat, REML = TRUE)

    # estimate BLUP
    modelblup <- lme4::ranef(lme)
    cat(paste("\tBLUP selected:",names(modelblup[1]),"\n"))
    
    # extract BLUP and add grand mean when only one repetition present
    BLUP_out[[colnames(dat)[i]]] <- modelblup[[1]] + summary(lme)$coefficients[1]
  }

  if(length(BLUP_out) > 0){
    for (i in 1:length(BLUP_out)) {
      temp <- as.data.frame(BLUP_out[[i]])
      temp$New <- row.names(temp)
      temp$id <- names(BLUP_out)[i]
      temp <- temp[,c(3, 2, 1)]
      colnames(temp) <- c("id", "Line", "Intercept")
      row.names(temp) <- seq(from = 1, to = nrow(temp), by = 1)
      BLUP_out[[i]] <- temp
      
      if(i==1){
        BLUP_out_df <- BLUP_out[[i]]
      } else{
        BLUP_out_df <- rbind(BLUP_out_df, BLUP_out[[i]])
      }
    }
    blup <- reshape2::dcast(BLUP_out_df, Line ~ id, value.var = "Intercept")
    colnames(blup)[1] <- colnames(dat)[random_effect[1]]
    row.names(blup) <- seq(from = 1, to = nrow(blup), by = 1)
  }
  
  not_converge_columns = c()
  for(i in 1:ncol(blup)){
    if(is.numeric(blup[,i])){
      if(var(blup[,i], na.rm = TRUE) < 1e-6){
        not_converge_columns = c(not_converge_columns, colnames(blup)[i])
      }
    }
  }
  
  if(!is.null(not_converge_columns)){
    blup = blup[, -(match(not_converge_columns, colnames(blup)))]
  }
  
  # Shapiro test BLUP
  blup_stats <- lapply(names(blup)[-1], function(i) {
    shap_BLUPs <- shapiro.test(blup[[i]]) 
    shap_BLUPs_P <- shap_BLUPs$p.value
    data.frame(Trait = i,
               Shapiro_BLUPs_pval = shap_BLUPs_P, 
               Shapiro_BLUPs_OK = ifelse(shap_BLUPs_P > 0.05, "TRUE", "FALSE"),
               stringsAsFactors = FALSE)
  }) %>% bind_rows %>% as.data.frame()
  
  # Store data   
  if(exists("blup")){
    return(
      list(
        "Outlier_removed_data" = outlier_removed_dat,
        "Outlier_data" = outlier_dat,
        "Outliers_residuals" = outliers_residuals,
        "Shapiro_raw_normtest" = shap_norm,
        "BLUP" = blup,
        "BLUP_shapiro" = blup_stats,
        "Not_converge_columns" = not_converge_columns
      )
    )
  } else{
    return(-1)
  }
}
