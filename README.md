# multi-VAR
The multiple-subject vector autoregressive (multi-VAR) model is a VAR-based high-dimensional time series framework designed to identify commonly shared paths (common components) and subject-specific paths (unique components) by decomposing the estimated entries of the VAR transition matrices (individual paths). The current version enhances the identification performance of common and unique paths by leveraging debiased Lasso estimators through nodewise regression and a communication-efficient data integration framework based on robust estimation of location parameters. In addition, it provides statistical tools to assess the nullity and homogeneity of individual paths, as well as the significance of the common paths. The corresponding manuscript is Kim, Fisher, and Pipiras (2025). This version also includes the R code used in the simulation studies and data application, along with the associated fMRI dataset.

## References

<div id="ref-multivarse">

Younghoon Kim, Zachary F. Fisher, and Vladas Pipiras. 2025.
“Joint modeling and inference of multiple-subject high-dimensional sparse vector autoregressive models.” <https://arxiv.org/abs/2510.14044>.

<div id="ref-multivar2">

Zachary F Fisher, Younghoon Kim, Vladas Pipiras, Christopher Crawford, Daniel J Petrie, Michael D Hunter, Charles F Geier. 2024.
“Structured Estimation of Heterogeneous Time Series.” *Multivariate Behavioral Research* 59 (6): 1270-1289.
<https://www.tandfonline.com/doi/abs/10.1080/00273171.2023.2283837>.

<div id="ref-multivar1">

Zachary F Fisher, Younghoon Kim, Barbara L Fredrickson, Vladas Pipiras. 2022.
“Penalized estimation and forecasting of multiple subject intensive longitudinal data.” *Psychometrika* 87 (2): 403-431.
<https://www.cambridge.org/core/journals/psychometrika/article/abs/penalized-estimation-and-forecasting-of-multiple-subject-intensive-longitudinal-data/A5933DB71583B091EB3D65CBA2E28498>.

