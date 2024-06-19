# ATOMIC: Analytic Tools for Omics Mining and Integrative Computation
Welcome to the GitHub page for the ATOMIC package!

# Overview
ATOMIC (Analytic Tools for Omics Mining and Integrative Computation) is a comprehensive R package designed to facilitate the analysis and integration of omics data. Whether you're working with genomics, transcriptomics, proteomics, or metabolomics data, ATOMIC provides a suite of powerful tools to streamline your analytical workflows.

# Features
ATOMIC is a new R package designed to smoothly integrate heterogeneous -datasets for efficient mining and analysis. Currently in its early phases, it provides functionalities specifically designed for analyzing massive datasets and is connected to the NIBMG, Kalyani LAN for beta testing. While still in development, its extensive tools promise to improve multiomic analysis. Future goals include making it available to all users, allowing a larger community to take advantage of its capabilities for complex data interpretation and integration.

# Installation
To install ATOMIC from GitHub, you'll need to use the remotes package. If you haven't already installed remotes, you can do so by running:

```R

install.packages("remotes")
#Once you have remotes installed, you can install ATOMIC by running the following command:

remotes::install_github("sb-genome/ATOMIC")
```
# Usage
Here's a quick example to get you started with using ATOMIC:

```R
# Load the ATOMIC package
library(ATOMIC)

# Load the dependances for the genre of functions that would be used.
g_initialize(c("Gtex_Summary" , "Haploreg")));
# g_initialize( c("Gtex_Summary", "1000GP_LD_Summary", "Gtex_LD_Summary" , "UKBB_BP_Enrichment" , "Haploreg" , "Genomic_Regions" ));

# Perform Gtex Summary Analysis
GTEx_Summary(SNP = "chr3_49773614_G_A_b38", Gene="ENSG00000182179.12", Tissue="Whole_Blood");

# Perform Haploreg Data Mining
snpToInfo(rsids = "rs10")
snpToInfo(rsids = c("rs10","rs3","rs4","rs10047249"))

# Close the connection to the Servers
g_free(connec)
```
For detailed documentation and examples, please refer to the official ATOMIC documentation.

# Contributing
We welcome contributions from the community to make ATOMIC even better! Whether it's bug fixes, new features, or documentation improvements, your contributions are valuable. Please see the contribution guidelines for more information on how to contribute.

# Contact
If you have any questions, feedback, or suggestions, feel free to contact us. We'd love to hear from you!
sb1@nibmg.ac.in , db2@nibmg.ac.in

Happy analyzing with ATOMIC! 🧬🔬
