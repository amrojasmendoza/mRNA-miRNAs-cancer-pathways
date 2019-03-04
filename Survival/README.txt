### What is this repository for? ###

Script to reproduce the correlation and survival analyses of the paper:


### How do I get set up? ###

Download the repository. 

Run the script runAnalysis.sh

The script requires R to be installed in the computer, and the following R libraries. 
If these libraries are not previously installed, the script will install them.

survival
data.table


Two folders will be created:

figures - containing 6 figures:
generalPairsCorrelations01.pdf : results of correlation analysis in 41 general gene-miRNA pairs
generalPairsCorrelations02.pdf : results of correlation analysis in 36-selected general gene-miRNA pairs
generalPairsSurvival.pdf       : results of survival analysis in 36-selected general gene-miRNA pairs
lungPairsCorrelations01.pdf    : results of correlation analysis in 74 lung-exclusive gene-miRNA pairs
lungPairsCorrelations02.pdf    : results of correlation analysis in 40 lung-exclusive gene-miRNA pairs
lungPairsSurvival.pdf          : results of survival analysis in 40 lung-exclusive gene-miRNA pairs

tables - containing 
general.correlation.contribution.csv : contributions to gene expression by CNA, methylation and miRNA
general.correlation.csv				 : correlation coefficients in 41 general gene-miRNA pairs
general.correlation.pval.csv	     : correlation P-values in 41 general gene-miRNA pairs
general.coxph.coefs.csv              : coxph coefficients in 36 selected general gene-miRNA pairs
general.coxph.pvals.csv              : coxph P-values in 36 selected general gene-miRNA pairs
lung.correlation.csv                 : correlation coefficients in 74 lung-exclusive gene-miRNA pairs
lung.correlation.pval.csv            : correlation P-values in 74 lung-exclusive gene-miRNA pairs
lung.coxph.coefs.csv                 : coxph coefficients in 40 selected lung-exclusives gene-miRNA pairs
lung.coxph.pvals.csv                 : coxph P-values in 40 selected lung-exclusives gene-miRNA pairs



### Who do I talk to? ###

Sergio Alonso
sergio.alonso.phd@gmail.com
====
Amrojasmendoza added a script to address issue 1. Credit of the code to sergio.alonso.phd@gmail.com (Sergio Alonso) 02-22-2019
