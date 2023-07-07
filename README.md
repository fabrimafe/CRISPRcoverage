# CRISPRcoverage

Repository for coverage analyses of large scale chromosomal rearrangements induced by Cas9/CRISPR as presented in 

**"A CRISPR-induced DNA break can trigger crossover, chromosomal loss and chromothripsis-like rearrangements"**, 
*Aviva Samach<sup>1</sup>, Fabrizio Mafessoni<sup>1</sup>, Or Gross<sup>1</sup>, Cathy Melamed-Bessudo, Shdema Filler-Hayut, Tal Dahan-Meir, Ziva Amsellem, Wojciech P. Pawlowski, Avraham A. Levy*

<sup>1</sup> equal contributors

biorxiv, https://doi.org/10.1101/2023.05.22.541757

## USAGE

**utilities_LSF** contains bash scripts that load modules on a LSF cluster and can be used to wrap the R scripts and to compute coverage using GATK

**R** contains R scripts for statistical tests and analyses

#### SCRIPTS

*script_cov.sh*
#script to compute genome-wide coverage using GATK
#you can run this by:
#./script_cov.sh myinput.bam output.tab chrX
#ARGUMENTS:
#myinput.bam: a bam file
#output.tab: the output file
#chrX: the chromosome for which to compute coverage

*allchr_coverageplot.R*
script to generate genome-wide coverage plots as in fig.6. of the manuscript

*chrarm_deletions.R*
script to test the significance of putative chromosomal arm deletions. The output is a list of statistics computed for different window-sizes, before and after the putative CRISPR cut sites, for different window-sizes before and after the cut site (1st column). A window-size of 0 indicates the full chromosome. For instance, for plant 1-1 described in the paper:

| window_size | coverage_pre/coverage_post | pvalue.complete.deletion | pvalue.1chr.deletion | pvalue.no.deletion | loglik.complete.deletion | loglik.1chr.deletion | loglik.no.deletion | best_model |
| ----------- | -------------------------- | ------------------------ | -------------------- | ------------------ | ------------------------ | -------------------- | ------------- | ---------- |
| 0 | 0.682996165094607 | 0 | 1 | 0 | -1009489.02991463 | -10571.6563606919 | -13328.352691852 | 1chr.deletion |
| 1000 | 0.458333333333333 | 4.06718181990874e-54 | 1 | 0.0148006050561844 | -125.949921546219 | -3.01327685661561 | -7.22636407348871 | 1chr.deletion|
| 20000 | 0.495229007633588 | 0 | 1 | 4.90826790182424e-29 | -2954.14520522009 | -29.730283024759 | -94.9143296115036 |1chr.deletion |
| 1e+06 | 0.598689201719632 | 0 | 8.34214803165462e-197 | 1 | -165796.633530543 | -2873.18983097139 | -2421.70188839253 | no.deletion |

The second column indicates the fold-dosage compared to the proximal part of the chromosome, which is ~0.5 for a 1-chromosome deletion. p-values are computed with a likelihood ratio test comparing the model with the highest log likelihood ( last column) versus the three different models, i.e. complete deletion, 1-chromosomal deletion and no deletion. The log-likelihood is reported in the loglik columns. Smallest window sizes better control for local coverage, but might have reduced power. Here, the 1Mb window size does not have sufficient power to reject a model with no rearrangement, while the full chromosome analysis and smaller window sizes indicate the loss of one chromosomal arm.

##### Repository structure

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

