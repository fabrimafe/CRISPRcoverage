# CRISPRcoverage

Repository for coverage analyses of large scale chromosomal rearrangements induced by Cas9/CRISPR as presented in 

**"A CRISPR-induced DNA break can trigger crossover, chromosomal loss and chromothripsis-like rearrangements"**, 
*Aviva Samach<sup>1</sup>, Fabrizio Mafessoni<sup>1</sup>, Or Gross<sup>1</sup>, Cathy Melamed-Bessudo, Shdema Filler-Hayut, Tal Dahan-Meir, Ziva Amsellem, Wojciech P. Pawlowski, Avraham A. Levy*

<sup>1</sup> equal contributors

biorxiv, https://doi.org/10.1101/2023.05.22.541757

### USAGE

**utilities_LSF** contains bash scripts that load modules on a LSF cluster and can be used to wrap the R scripts and to compute coverage using GATK

**R** contains R scripts for statistical tests and analyses

.
 * [R](./R)
   * [allchr_coverageplot.R](./R/allchr_coverageplot.R)
   * [chrarm_deletions2table.R](./R/chrarm_deletions2table.R)
   * [chrarm_deletions2table_tt4.R](./R/chrarm_deletions2table_tt4.R)
   * [chrarm_deletions_plotcustome.R](./R/chrarm_deletions_plotcustome.R)
   * [chrarm_deletions_plot.R](./R/chrarm_deletions_plot.R)
   * [chrarm_deletions.R](./R/chrarm_deletions.R)
   * [chrarm_deletions_scan.R](./R/chrarm_deletions_scan.R)
 * [utilities_LSF](./utilities_LSF)
   * [chrarm_deletions2table.sh](./utilities_LSF/chrarm_deletions2table.sh)
   * [script_cov.sh](./utilities_LSF/script_cov.sh)
 * [README.md](./README.md)

