# =======================================================
# Script: regress.R
# Purpose: Residualize predictors (X) and outcome (Y) by 
#          regressing out covariates.
#          This file includes two distinct functions:
#          (1) Remove covariate effects from each predictor 
#              connection (X).
#          (2) Remove covariate effects from the outcome (Y).
#
# Author: Majid Saberi, University of Michigan
# Contact: majidsa@umich.edu
# Date: 2025-09-19
# Version: 1.0
#
# Dependencies:
#   - base R
#
# Citation:
#   If you use this code, please cite:
#   Saberi et al. Neural synchrony in the pain connectome is 
#   associated with pain severity and its interactions with 
#   mental health outcomes: a transdiagnostic study using 
#   magnetoencephalography and multivariate modeling. 2025.
# =======================================================


# -------------------------------------------------------
# Function: residualize_X()
# -------------------------------------------------------
# Purpose:
#   Remove the influence of covariates from each predictor 
#   connection (column of X).
#
# Use case:
#   - Adjust brain connectivity features (X) for nuisance 
#     factors before multivariate modeling.
#
# Arguments:
#   X          : data frame of predictors 
#                [subjects x connections]
#   covariates : data frame of covariates [subjects x variables]
#
# Returns:
#   - Data frame containing residualized X (same dimensions)
# -------------------------------------------------------
residualize_X <- function(X, covariates) {
  X_res <- matrix(NA, nrow = nrow(X), ncol = ncol(X))
  colnames(X_res) <- colnames(X)
  rownames(X_res) <- rownames(X)
  
  # Loop over each connection (column of X)
  for (j in 1:ncol(X)) {
    df <- data.frame(covariates, Conn = X[, j])
    mdl <- lm(Conn ~ ., data = df)      # regress connection on covariates
    X_res[, j] <- residuals(mdl)        # extract residuals
  }
  
  return(as.data.frame(X_res))
}


# -------------------------------------------------------
# Function: residualize_Y()
# -------------------------------------------------------
# Purpose:
#   Remove the influence of covariates from the outcome (Y).
#
# Use case:
#   - Adjust the response variable (e.g., pain score) for 
#     nuisance factors such as age, sex, or site.
#
# Arguments:
#   Y          : numeric response vector
#   covariates : data frame of covariates [subjects x variables]
#
# Returns:
#   - Data frame containing the residualized Y
# -------------------------------------------------------
residualize_Y <- function(Y, covariates) {
  mdl <- lm(Y ~ ., data = covariates)  # regress Y on covariates
  Y_res <- residuals(mdl)              # extract residuals
  return(as.data.frame(Y_res))
}
