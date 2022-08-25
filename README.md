# LAIP-Simulator
R-based data generator to simulate blast-events with user defined marker expressions.
The simulated data can be integrated in an existing fcs-file to simulate a MRD-positive bone matter sample. 

Installation:
This package uses several Bioconductor packages that are not available in CRAN-repository. They can be installed using the BiocManager.

install.packages("BiocManager")
BiocManager::install(c("flowCore", "Biobase", "BiocGeneric"))