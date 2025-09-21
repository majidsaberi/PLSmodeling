# =======================================================
# Script: reformat.R
# Purpose: Reformat input data and results.
#          This file includes two distinct functions:
#          (1) Reformat connectivity arrays to subject x connection 
#              data frames (input for statistical modeling).
#          (2) Reconstruct ROI x ROI matrices from connection 
#              vectors (e.g., model coefficients or results).
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
# Function: connectivity_to_df()
# -------------------------------------------------------
# Purpose:
#   Used at the *start* of analysis.
#   Converts a 3D connectivity array [roi x roi x subj] 
#   into a subject x connection data frame.
#
# Use case:
#   - Input for regression/PLS modeling, where each column 
#     is a unique connection and each row is a subject.
#
# Arguments:
#   conn : 3D array [roi, roi, subj]
#
# Returns:
#   data frame [subjects x connections]
# -------------------------------------------------------
connectivity_to_df <- function(conn) {
  data_conn <- t(apply(conn, 3, function(x) x[upper.tri(x)])) # selected upper-triangular indices per subject
  return(data_conn)
}


# -------------------------------------------------------
# Function: vec_to_sqmatrix()
# -------------------------------------------------------
# Purpose:
#   Used at the *end* of analysis.
#   Converts a vector of connection values (e.g., PLS regression 
#   coefficients or p-values) back 
#   into a full symmetric adjacency matrix [roi x roi].
#
# Use case:
#   - Map model results (coefficients) back to ROI space
#     for visualization in connectome/brain plots.
#
# Arguments:
#   vec   : numeric vector of upper-triangular connection values
#   n_roi : number of ROIs
#
# Returns:
#   symmetric adjacency matrix [n_roi x n_roi]
# -------------------------------------------------------
vec_to_sqmatrix <- function(vec, n_roi) {
  mat <- matrix(0, n_roi, n_roi)   # initialize empty matrix
  ut_idx <- which(upper.tri(mat))  # find upper-triangular indices
  mat[ut_idx] <- vec               # fill upper-triangle
  mat <- mat + t(mat)              # mirror to lower-triangle
  diag(mat) <- NA                  # remove diagonal values
  return(mat)
}
