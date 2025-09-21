# PLS Modeling of Brain Connectivity and Clinical Outcomes

This repository provides R scripts for applying **Partial Least Squares (PLS) regression** to identify associations between **brain connectivity features** (independent variables) and **continuous clinical outcomes** (dependent variables, such as pain severity).  

### Why PLS Modeling?
PLS regression is particularly well-suited for neuroimaging and clinical data because:  

- **Handles high-dimensional, collinear predictors**  
  Brain connectivity matrices often contain hundreds of connections, many of which are strongly correlated. PLS can effectively reduce dimensionality while preserving meaningful variance.  

- **Integrates brain–behavior associations**  
  Unlike univariate tests, PLS simultaneously considers the entire connectivity pattern, capturing distributed neural synchrony linked to symptom burden.  

- **Compatible with resampling-based inference**  
  The method integrates well with permutation and bootstrap tests, providing empirical p-values that strengthen interpretability in small or heterogeneous clinical samples.  

---

### Repository Structure

R/
├── reformat.R     # Convert connectivity arrays ↔ data frames / matrices
├── regress.R      # Residualize predictors (connectivity) and outcome (clinical variable)
├── pls_model.R    # Fit PLS model, extract coefficients, resampling p-values
├── main_run.R     # End-to-end pipeline runner
Data/              # (Optional) demo dataset
README.md          # Project documentation
LICENSE            # Usage license


### Citation

If you use this code, please cite:

Saberi et al. (2025). Neural synchrony in the pain connectome is associated with pain severity and its interactions with mental health outcomes: a transdiagnostic study using magnetoencephalography and multivariate modeling.
