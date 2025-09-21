# =======================================================
# Script: pls_model.R
# Purpose: Implement Partial Least Squares (PLS) regression 
#          for brain connectivity analysis.
#          This file includes two distinct functions:
#          (1) Fit a PLS model and extract connection-level 
#              coefficients.
#          (2) Perform resampling tests (permutation/bootstrap) 
#              to compute p-values for coefficients.
#
# Author: Majid Saberi, University of Michigan
# Contact: majidsa@umich.edu
# Date: 2025-09-19
# Version: 1.0
#
# Dependencies:
#   - pls
#
# Citation:
#   If you use this code, please cite:
#   Saberi et al. Neural synchrony in the pain connectome is 
#   associated with pain severity and its interactions with 
#   mental health outcomes: a transdiagnostic study using 
#   magnetoencephalography and multivariate modeling. 2025.
# =======================================================

library(pls)


# -------------------------------------------------------
# Function: coef_pls()
# -------------------------------------------------------
# Purpose:
#   Fit a PLS regression model and extract coefficients for 
#   brain connections (predictors).
#
# Use case:
#   - Get connection-level weights linking connectivity (X) 
#     to outcome (Y).
#   - Supports different cross-validation strategies.
#
# Arguments:
#   data       : data frame with predictors (X) and outcome (Y)
#   cv_method  : cross-validation method
#                ("none", "LOO" = leave-one-out, "CV" = k-fold)
#   segments   : number of CV folds if cv_method = "CV"
#   max_comp   : maximum number of components to evaluate
#
# Returns:
#   - Numeric vector of connection-level coefficients from 
#     the optimal number of components.
# -------------------------------------------------------
coef_pls <- function(data, cv_method = c("none", "LOO", "CV"),
                     segments = 10, max_comp = 10) {
  cv_method <- match.arg(cv_method)
  
  # Fit model depending on CV method
  model <- switch(cv_method,
                  "none" = plsr(Y ~ ., data = data, ncomp = max_comp,
                                scale = TRUE, center = TRUE, validation = "none"),
                  "LOO"  = plsr(Y ~ ., data = data, ncomp = max_comp,
                                scale = TRUE, center = TRUE, validation = "LOO"),
                  "CV"   = plsr(Y ~ ., data = data, ncomp = max_comp,
                                scale = TRUE, center = TRUE, validation = "CV",
                                segments = segments)
  )
  
  # Select optimal number of components by RMSEP
  opt_comp <- which.min(RMSEP(model)[[1]][1, , -1])
  
  # Extract coefficients at optimal components
  coefs <- model$coefficients[,,opt_comp]
  return(coefs)
}


# -------------------------------------------------------
# Function: coef_test_pls()
# -------------------------------------------------------
# Purpose:
#   Assess statistical significance of PLS connection-level 
#   coefficients using resampling (permutation or bootstrap).
#
# Use case:
#   - Estimate empirical p-values for coefficients.
#   - Test robustness of PLS model against null distributions.
#
# Arguments:
#   data          : predictors (X) + outcome (Y)
#   cv_method     : "none", "LOO", or "CV"
#   segments      : folds if cv_method = "CV"
#   max_comp      : maximum number of components
#   n_resamples   : number of permutations/bootstrap samples
#   resample_type : "permutation" or "bootstrap"
#
# Returns:
#   - Numeric vector of p-values for each connection coefficient
# -------------------------------------------------------
coef_test_pls <- function(data, cv_method = c("none", "LOO", "CV"),
                          segments = 10, max_comp = 10,
                          n_resamples = 100,
                          resample_type = c("permutation", "bootstrap")) {
  cv_method <- match.arg(cv_method)
  resample_type <- match.arg(resample_type)
  
  # 1. Fit actual model
  coef_actual <- coef_pls(data, cv_method, segments, max_comp)
  
  # 2. Resampling loop
  coef_resampled <- matrix(NA, nrow = n_resamples, ncol = length(coef_actual))
  for (i in 1:n_resamples) {
    if (resample_type == "permutation") {
      # Shuffle Y while keeping X intact
      data_resample <- data
      data_resample$Y <- sample(data$Y)
    } else {
      # Bootstrap subjects (rows)
      idx <- sample(1:nrow(data), replace = TRUE)
      data_resample <- data[idx, ]
    }
    coef_resampled[i, ] <- coef_pls(data_resample, cv_method, segments, max_comp)
  }
  
  # 3. Compute empirical p-values
  p_values <- sapply(1:length(coef_actual), function(j) {
    mean(abs(coef_resampled[, j]) >= abs(coef_actual[j]))
  })
  
  return(p_values)
}
