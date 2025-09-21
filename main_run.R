# =======================================================
# Script: main_run.R
# Purpose: End-to-end pipeline runner for PLS modeling
#          including reformatting, residualization, 
#          modeling, and statistical testing.
#
# Author: Majid Saberi, University of Michigan
# Contact: majidsa@umich.edu
# Date: 2025-09-19
# Version: 1.0
#
# Dependencies:
#   - pls
#   - base R
#
# Input:
#   - Connectivity arrays (e.g., wPLI_Conn)
#   - Phenotype/covariate data frame
#
# Output:
#   - Coefficient and p-value matrices (ROI x ROI)
#
# Usage:
#   source("R/main_run.R")
#   # Example inside script demonstrates full workflow
#
# Citation:
#   If you use this code, please cite:
#   Saberi et al. Neural synchrony in the pain connectome is 
#   associated with pain severity and its interactions with 
#   mental health outcomes: a transdiagnostic study using 
#   magnetoencephalography and multivariate modeling. 2025.
# =======================================================


# -------------------------------
# Step 0. Setup
# -------------------------------
setwd("your working directory")

# Load custom functions
source("reformat.R")
source("regress.R")
source("pls_model.R")


# -------------------------------
# Step 1. Load your data
# -------------------------------
# Example assumes an RData file containing:
#   - wPLI_Conn: 4D connectivity array [roi, roi, condition, subject]
#   - Pheno:     phenotype/covariate data frame
load("your RData file") 

# Select condition-specific connectivity (example: first band)
Connec <- wPLI_Conn[,,1,]  #Example
# This produces an ROI x ROI x subject array (e.g., 36 x 36 x N)


# -------------------------------
# Step 2. Convert connectivity to data frame
# -------------------------------
# Convert array (ROI x ROI x subj) → data frame (subjects x connections).
# If your connectivity is already in subj x connection format, skip this step.
X <- connectivity_to_df(Connec)

# Assign outcome variable (dependent variable)
# Here: Pain, but you could substitute any phenotype or combined scores
Y <- Pheno$Pain

# -------------------------------
# Step 3. Residualize against covariates (optional)
# -------------------------------
# Define covariates dateframe (example: age and sex are columns 5 and 6), must be numerical
Z <- Pheno[,c("Age","Sex")]

# Residualize predictors (connections) and outcome
X_ <- residualize_X(X, Z)
Y_ <- residualize_Y(Y, Z)


# -------------------------------
# Step 4. Merge dependent and independent variables
# -------------------------------
# Combine into one dataset for PLS regression
data <- data.frame(X_, Y_)
colnames(data)[ncol(data)] <- "Y"


# -------------------------------
# Step 5. Run PLS regression + resampling tests
# -------------------------------

# ---- Define analysis parameters ----
# Cross-validation method:
#   "none" = no CV
#   "LOO"  = leave-one-out CV
#   "CV"   = k-fold CV (requires 'segments')
cv_method <- "CV"

# Number of folds if cv_method = "CV"
segments <- 10

# Maximum number of PLS components to evaluate
max_comp <- 5

# Resampling parameters
n_resamples <- 100                   # number of permutation/bootstrap samples
resample_type <- "permutation"       # "permutation" or "bootstrap"

# ---- Run PLS regression ----
# Fit PLS model and extract connection-level coefficients
coefs <- coef_pls(
  data, 
  cv_method = cv_method, 
  segments  = segments, 
  max_comp  = max_comp
)

# ---- Resampling-based testing ----
# Assess significance of coefficients using chosen resampling strategy
pvals <- coef_test_pls(
  data, 
  cv_method     = cv_method, 
  segments      = segments,
  max_comp      = max_comp, 
  n_resamples   = n_resamples,
  resample_type = resample_type
)

# ---- Convert results to matrices ----
# (adjust n_roi = number of brain regions in your data)
n_roi <- 36
coefs_mat <- vec_to_sqmatrix(coefs, n_roi)
pvals_mat <- vec_to_sqmatrix(pvals, n_roi)

# -------------------------------
# Step 6. Save results
# -------------------------------
write.csv(coefs_mat, "pls_coefficients_matrix.csv")
write.csv(pvals_mat, "pls_pvalues_matrix.csv")
