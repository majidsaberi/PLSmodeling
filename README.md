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

* `reformat.R`: Converts connectivity data between arrays, data frames, and matrices.
* `regress.R`: Residualizes predictors (connectivity) and the outcome (a clinical variable).
* `pls_model.R`: Fits the PLS model, extracts coefficients, and calculates resampling p-values.
* `main_run.R`: The primary script to run the entire analysis pipeline from start to finish.
* `README.md`: This documentation file.

### Citation

If you use this code, please cite:

Saberi et al. (2025). Neural synchrony in the pain connectome is associated with pain severity and its interactions with mental health outcomes: a transdiagnostic study using magnetoencephalography and multivariate modeling.
