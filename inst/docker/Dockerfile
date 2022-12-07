FROM bioconductor/bioconductor_docker:latest
RUN R -q -e 'install.packages("remotes")'
RUN R -q -e 'remotes::install_github("RGLab/ImmuneSignatures")'
RUN R -q -e 'remotes::install_url("https://www.bioconductor.org/packages//2.10/bioc/src/contrib/DESeq_1.8.3.tar.gz")'

