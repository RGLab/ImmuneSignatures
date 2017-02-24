About ImmSigPkg
=========
  
*Last Updated:* February 2017  
*Maintainer:* Evan Henrich  
*Contact:* ehenrich@fredhutch.org

A package of functions that allow the user to run the entire HIPC ImmuneSignatures Project 
Data Analysis Pipeline with original parameters from the published manuscript or with 
user-defined options.

Getting Started
===============

**To download and install the package for interactive work**

```R
library(devtools)
install_github("ehfhcrc/ImmSigPkg")
library(ImmSigPkg)
```

**To see vignette:** 
```R
# in console
git clone https://github.com/ehfhcrc/ImmSigPkg.git

# in R
setwd(<path_to_ImmSigPkg_source_code_dir>)
library(devtools)
build_vignettes()
# open vignette at <working_dir>/ImmSigPkg/Inst/Doc/Manuscript_Pipeline.html
```

Note: the vignette shows the output of the pipeline using the original parameters and takes approximately 
45 minutes to build.  At this time, it is recommended to build from source without compression by cloning 
the repo locally because there are unexplained differences between output for the old cohort when compressed 
versus built / knit directly. This may be due to an R environmental option used in the `qusage()` call within `meta_analysis()`.

Running the Pipeline
====================

To run the entire pipeline you can execute the function `hipc_full_pipeline()` and
follow the prompts.  It is recommended to use the original parameters for the first run
to ensure things flow smoothly.  This may take about 1 hour on your local machine. 
Once you have run everything through once and generated rawdata files (.txt) as well 
as .rds expressionSet files, you may run the meta analysis with different parameters 
using the command `hipc_meta_analysis()`.  One point, this will write over your current 
result files.  

**Notes as of February 2017:**  

1. SDY80 rawdata is not yet able to be processed directly, pipeline uses original 
rawdata files provided by collaborators. Therefore, selecting "process SDY80 rawdata" 
when prompted will error out.  This functionality will be incorporated as time permits.  
2. SDY400 genetic expression data is not publicly available and the pipeline uses
rawdata for baseline expression provided by the collaborators.  This will be made public
in Summer 2017 and will be incorporated as time permits.
