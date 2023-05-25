# causal-mediation-case-study
Brief introduction

case_CHR_MDS.R is the main file to run MDS on real data for each CHR. 3 paremeters must be provided: Chr(index of CHR), MdsTime(number of data splitting), Fdr(FDR level). 

Example:

Rscript --vanilla case_CHR_MDS.R --Chr 1 --MdsTime 250 --Fdr 0.05
