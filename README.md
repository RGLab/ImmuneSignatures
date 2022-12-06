# ImmuneSignatures

[![R-CMD-check](https://github.com/RGLab/ImmuneSignatures/workflows/R-CMD-check/badge.svg)](https://github.com/RGLab/ImmuneSignatures/actions)
[![docker](https://github.com/RGLab/ImmuneSignatures/actions/workflows/docker-build.yaml/badge.svg)](https://hub.docker.com/r/rglab/immunesignatures)


The purpose of this package is to reproduce the original results from [the HIPC ImmuneSignatures study](http://immunology.sciencemag.org/content/2/14/eaal4656) on baseline transcriptional predictors of response to the Influenza vaccine in cohorts of young and old adults.  

Before installing this package, you will need to install five dependencies only available on bioconductor: qusage, metaIntegrator, preProcessorCore, DESeq, BioBase. You can do this with the `biocInstaller::biocLite("pkg_to_install")` function.

Once the bioconductorTo see a report with the results, you can first install the package via `devtools::install_github("rglab/ImmuneSignatures", build_vignettes = TRUE)` then view the vignette with `browseVignettes(package = "ImmuneSignatures")`.  Note that it can take about 30 minutes to build the vignette on your local machine.

Since the code used to generate the vignette is a heavily refactored version of the original code, the original code is provided within the package and can be viewed on [GitHub](https://github.com/rglab/ImmuneSignatures/tree/main/inst/OrigCode). 
