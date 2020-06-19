# RA-OA

The python version used is 3.7.6, and the R version is 3.6.1

ATAC-seq:

1. ATAC_PreProcess_Disease_ImmuneCells.ipynb: ATAC-seq normalization and significant difference peak analysis

2. FisherExactTest.py: Fisher's exact test is used to calculate p to test the significance of associations between gene sets

3. KmeansCluster.py: Unsupervised clustering of differential peaks 

4. LinearRegressionOLS.py: We used ordinary least squares (OLS) using statsmodels in Python to performed linear regression as an approach to measure the correlation of RAAS with the clinical data of RA patients including DAS28-CRP, DAS28-ESR, CRP, ESR, TCJ, SCJ, RF, and Anti-CCP. 

5. PairwiseSigPeakFilter.py: Differential analysis. For ATAC-seq, each cell type was compared with all other cell types, and cell type-specific peaks were filtered with abs(log2 fold change) > 4, p < 0.001, and FDR < 0.01. 

RNA-seq:

6. DEseq.r: RNA-seq normalization

7. RNA_PreProcess_Disease_ImmuneCells.ipynb: RNA-seq Significant difference peak analysis

8. GSEAanalysis.R: To calculate the enrichment score and p of the gene set downloaded from MSigDB in the gene matrix. The values of p less than 0.05 are considered to be significantly enriched. 

DATA:

The “DATA” folder contains the raw counts matrix and peaks bed file required for 'ATAC_PreProcess_Disease_ImmuneCells.ipynb', 'DEseq.r' and 'RNA_PreProcess_Disease_ImmuneCells.ipynb'.
