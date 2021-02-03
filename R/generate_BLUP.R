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
  
  `%>%` <- dplyr::`%>%`
  
  # Check if input arguments are given
  if (missing(dat)) {
    stop("Data needs to be provided")
  }
  
  if (missing(random_effect)) {
    stop("Random effect column(s) needs to be provided")
  } else if (!is.numeric(random_effect) & any(random_effect < 1) | length(random_effect) == 0) {
    stop("Random effect needs to be > 0")
  }
  
  if (missing(sample_column)) {
    stop("Sample column needs to be provided")
  } else if (!is.numeric(sample_column) & sample_column < 1) {
    stop("Sample column needs to be > 0")
  }
  
  if (missing(start_column)) {
    stop("Start column needs to be provided")
  } else if (!is.numeric(start_column) & start_column < 1) {
    stop("Start column needs to be > 0")
  }
  
  if (!is.null(fixed_effect)) {
    if(!is.numeric(fixed_effect) & any(fixed_effect < 1) | length(fixed_effect) == 0) {
      stop("Fixed effect column number(s) must be integers > 0")
    }
  }
  
  if (col_convert) {
    message("\tConverting lmer model effect columns to factor type")
    dat <- dplyr::mutate(dat, across(all_of(random_effect), as.factor))
    if (!is.null(fixed_effect)) {
      dat <- dplyr::mutate(dat, across(all_of(fixed_effect), as.factor))
    }
  }
  
  # Convert the response columns to numeric
  dat <- dplyr::mutate(dat, across(seq(start_column, ncol(dat)), as.numeric))
  
  # Extract fixed effect name(s) (if provided)
  if (!is.null(fixed_effect)) {
    fixed_ef <- colnames(dat)[fixed_effect]
  }
  
  # Use this sanity check to see if column names are consistent throughout the analysis.
  start_colnames <- colnames(dat)
  
  #######################################################################
  ## Outlier Removal
  #######################################################################
  message("Remove outliers")
  
  # Create lmer formula
  termlabels <- c(paste0("(1|", colnames(dat)[random_effect], ")"))
  if (length(random_effect) > 0 & !is.null(fixed_effect)) { # with fixed effect
    termlabels <- c(fixed_ef, termlabels)
  }
  
  message(paste("\tlmer model, effect included in formula:", paste(termlabels, collapse = ", ")))
  
  # Find outliers
  trait_names <- colnames(dat)[start_column:ncol(dat)]
  
  outliers_residuals <- purrr::map(setNames(trait_names, trait_names), function(t) {
    res <- car::outlierTest(
      lme4::lmer(formula = reformulate(termlabels = termlabels, response = t),
                 data = dat, REML = TRUE),
      n.max = nrow(dat),
      cutoff = 0.05
    )
    as.integer(names(res$bonf.p)[res$bonf.p < 0.05])
  })
  
  outlier_dat <- purrr::map2_dfr(names(outliers_residuals), outliers_residuals, ~ {
    if (length(.y) > 0) {
      dplyr::select(dat, all_of(sample_column), all_of(.x)) %>% dplyr::slice(.y)
    }
  })
  
  purrr::walk2(names(outliers_residuals), outliers_residuals, ~ {
    if (length(.y) > 0) {
      dat <<- dplyr::mutate(dat, across(all_of(.x), function(column) replace(column, .y, NA)))
    }
  })

  outlier_removed_dat <- dat
  
  # Sanity check 
  if(!identical(start_colnames, colnames(dat))){
    stop("Column names are changed. The column assumption no longer hold.")
  }
  
  #######################################################################
  # orderNorm
  #######################################################################
  message("Calculate orderNorm transformed")
  
  # Store shapiro test
  shap_norm <- purrr::map_df(trait_names, function(t) {
    shap_res <- shapiro.test(dplyr::pull(dat, t))
    tibble(Trait = t,
           Shapiro_pval = shap_res$p.value,
           Shapiro_OK = Shapiro_pval > 0.05)
  })
  
  # Extract names of trait in need of transfromation
  trans_trait <- shap_norm$Trait[!shap_norm$Shapiro_OK]
  
  # Transform those traits
  ## transform 
  dat <- dplyr::mutate(dat, dplyr::across(all_of(trans_trait), ~ {
    output <- purrr::quietly(bestNormalize::orderNorm)(.)
    if (length(output$warnings) > 0) {
      warning(paste0(dplyr::cur_column(), ": ", output$warnings))
    }
    output$result$x.t
  }))
  
  # Sanity check 
  if(!identical(start_colnames, colnames(dat))){
    stop("Column names are changed. The column assumption no longer hold.")
  }
  
  #######################################################################
  ## generate BLUP
  #######################################################################
  message("Calculate BLUPs")
  message(paste("\tlmer model, effect included in formula:", paste(termlabels, collapse = ", ")))
  
  # Check all columns for NaN and Inf
  dat <- dplyr::mutate(dat, across(seq(start_column, ncol(dat)),
                                   ~ ifelse(is.infinite(.) | is.nan(.), NA, .)))
  
  # fit the model
  blup <- purrr::map2_dfc(seq_along(trait_names), setNames(trait_names, trait_names), ~ {
    lme <- lme4::lmer(formula = reformulate(termlabels = termlabels, response = .y), data = dat, REML = TRUE)
    # estimate BLUP
    modelblup <- lme4::ranef(lme)
    message(paste0("\tBLUP selected (", .y, "): ", names(modelblup[1])))
    
    # extract BLUP and add grand mean when only one repetition present
    if (.x == 1) {
      # Only add a genotype column for the first trait in order to not have
      # more than one of them in the final tibble.
      dplyr::tibble(Genotype = rownames(modelblup[[1]]),
                    !!(.y) := modelblup[[1]][[1]] + summary(lme)$coefficients[1])
    } else {
      dplyr::tibble(!!(.y) := modelblup[[1]][[1]] + summary(lme)$coefficients[1])
    }
  })
  
  blup_var <- purrr::map_dbl(dplyr::select(blup, -Genotype), var, na.rm = TRUE)
  not_converge_columns <- names(blup_var)[blup_var < 1e-6]
  blup <- dplyr::select(blup, -all_of(not_converge_columns))
  
  # Shapiro test BLUP
  purrr::map_df(dplyr::select(blup, -Genotype), ~ {
    tibble(Shapiro_BLUPs_pval = shapiro.test(.)$p.value,
           Shapiro_BLUPs_OK = Shapiro_BLUPs_pval > 0.05)
  }, .id = "Trait")

  list(
    "Outlier_removed_data" = outlier_removed_dat,
    "Outlier_data" = outlier_dat,
    "Outliers_residuals" = outliers_residuals,
    "Shapiro_raw_normtest" = shap_norm,
    "BLUP" = blup,
    "BLUP_shapiro" = blup_stats,
    "Not_converge_columns" = not_converge_columns
  )
}
