# CRISPRcoverage

Repository for coverage analyses of large scale chromosomal rearrangements induced by Cas9/CRISPR as presented in 

"A CRISPR-induced DNA break can trigger crossover, chromosomal loss and chromothripsis-like rearrangements", Aviva Samach*, Fabrizio Mafessoni*, Or Gross*, Cathy Melamed-Bessudo, Shdema Filler-Hayut, Tal Dahan-Meir, Ziva Amsellem, Wojciech P. Pawlowski, Avraham A. Levy

*equal contributors

biorxiv, https://doi.org/10.1101/2023.05.22.541757

### USAGE

utilities_LSF contains bash scripts that load modules on a LSF cluster and can be used to wrap the R scripts and to compute coverage using GATK

R contains R scripts for statistical tests and analyses

.
├── R
│   ├── allchr_coverageplot.R
│   ├── chrarm_deletions2table.R
│   ├── chrarm_deletions2table_tt4.R
│   ├── chrarm_deletions_plotcustome.R
│   ├── chrarm_deletions_plot.R
│   ├── chrarm_deletions.R
│   └── chrarm_deletions_scan.R
├── README.md
└── utilities_LSF
    ├── chrarm_deletions2table.sh
    └── script_cov.sh

