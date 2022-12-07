FROM bioconductor/bioconductor_docker:latest
RUN R -q -e 'install.packages("remotes")'
RUN R -q -e 'remotes::install_github("RGLab/ImmuneSignatures")'
