About ImmuneSignatures
=========
  
*Last Updated:* May 2017  
*Maintainer:* Evan Henrich, Sys. Analyst / Programmer, Gottard Lab @ Fred Hutchinson Cancer Research Center 
*Contact:* ehenrich@fredhutch.org

The purpose of this package is to reproduce the original results from the HIPC ImmuneSignatures study on baseline transcriptional predictors of response to the Influenza vaccine in cohorts of young and old adults.  

Before installing this package, you will need to install five dependencies only available on bioconductor: qusage, metaIntegrator, preProcessorCore, DESeq, BioBase. You can do this with the `biocInstaller::biocLite("pkg_to_install")` function.

Once the bioconductorTo see a report with the results, you can first install the package via `devtools::install_github("rglab/ImmuneSignatures", build_vignettes = TRUE)` then view the vignette with `browseVignettes(package = "ImmuneSignatures")`.  Note that it can take about 30 minutes to build the vignette on your local machine.

Since the code used to generate the vignette is a heavily refactored version of the original code, the original code is provided within the package and can be viewed on github at <https://github.com/rglab/ImmuneSignatures/tree/main/OrigCode>. 
