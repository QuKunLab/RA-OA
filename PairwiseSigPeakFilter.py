# coding=utf-8
# -*- coding: utf-8 -*-
"""
Created on 2020.01.08

@author: Cecilia
"""


import pandas as pd
import numpy as np
import scipy.stats
import os
#Define f(x)


def fdr_BH(p_vals):
    from scipy.stats import rankdata
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1
    return fdr



F='/home/qukun/huangbb/test/AllSample_PairwiseCompare/ATACseq_QNorm_All.log2.txt'


QNormDFlog=pd.read_table(F,sep='\t',index_col=0)


#del QNormDFlog['RAMO124']

#BF='/home/huangbeibei/test/filter/PeakBed.removeY.bed'
#BedDF.index=BedDF[3]

Samples=list(QNormDFlog)

RABSamples=[i for i in Samples if ('RAB' in i )]
OABSamples=[i for i in Samples if ('OAB' in i )]
HCBSamples=[i for i in Samples if ('HCB' in i )]
RAMSamples=[i for i in Samples if ('RAM' in i)]
OAMSamples=[i for i in Samples if ('OAM' in i)]
HCMSamples=[i for i in Samples if ('HCM' in i)]
RAT4Samples=[i for i in Samples if ('RATA' in i )]
OAT4Samples=[i for i in Samples if ('OATA' in i )]
HCT4Samples=[i for i in Samples if ('HCTA' in i )]
RAT8Samples=[i for i in Samples if ('RATB' in i )]
OAT8Samples=[i for i in Samples if ('OATB' in i )]
HCT8Samples=[i for i in Samples if ('HCTB' in i )]



def PairwiseCompare(ASamples,BSamples,logfc,pval,fdr):
    Dir='/home/qukun/huangbb/test/AllSample_PairwiseCompare/'
    outDir=os.path.join(Dir,'logfc{}_p{}_q{}_Batch21/'.format(logfc,pval,fdr))
    if not os.path.exists(outDir)
        os.mkdir(DirX)

    #QNormBedDF=BedDF.loc[list(QNormDFlog.index)]
    #OtherSamples=list(set(Samples).difference(set(CellSamples)))
    ACount=QNormDFlog[ASamples]
    BCount=QNormDFlog[BSamples]
    AvsB_FD=np.abs(ACount.apply(np.mean,axis=1)-BCount.apply(np.mean,axis=1))>logfc
    #AvsB_FD_up=(ACount.apply(np.mean,axis=1)-BCount.apply(np.mean,axis=1))>logfc
    #AvsB_FD_down=(BCount.apply(np.mean,axis=1)-ACount.apply(np.mean,axis=1))>logfc
    AvsB_PVAL=(pd.Series(scipy.stats.ttest_ind(ACount,BCount,axis=1)[1],index=ACount.index))<pval
    AvsB_PVAL1=pd.Series(scipy.stats.ttest_ind(ACount,BCount,axis=1)[1],index=ACount.index)
    AvsB_QVVAL=(pd.Series(fdr_BH(AvsB_PVAL1),index=ACount.index))<fdr
    #FilterDF_up=QNormDFlog[AvsB_FD_up & AvsB_PVAL & AvsB_QVVAL]
    #FilterDF_down=QNormDFlog[AvsB_FD_down & AvsB_PVAL & AvsB_QVVAL]
    FilterDF=QNormDFlog[AvsB_FD & AvsB_PVAL & AvsB_QVVAL]
    FilterDF.to_csv('{}_vs_{}_sigpeak.txt'.format(logfc,pval,fdr,ASamples[0][0:4],BSamples[0][0:4]),sep='\t')



LISTsamples=[RABSamples,OABSamples,HCBSamples,RAMSamples,OAMSamples,HCMSamples,RAT4Samples,OAT4Samples,HCT4Samples,RAT8Samples,OAT8Samples,HCT8Samples]
for i in LISTsamples:
    for j in LISTsamples:
        if i!=j:
            PairwiseCompare(i,j,4.0,0.001,0.01)

#awk '{print $0}' *.txt > combine.txt
#awk '!x[$0]++' combine.txt > TwoTwoCompare_Merge.txt

#RAvsOA
#RAvsOA
#CellTypeSpecific_Number(RAs,OAs,Rname,Oname,1.0,0.05,0.05)
#CellTypeSpecific_Number(RAs,OAs,Rname,Oname,1.0,0.001,0.1)
#CellTypeSpecific_Number(RAs,OAs,Rname,Oname,1.0,0.01,0.1)
#RAvsHC
#CellTypeSpecific_Number(RAs,HCs,Rname,Hname,1.0,0.05,0.05)
#CellTypeSpecific_Number(RAs,HCs,Rname,Hname,1.0,0.001,0.1)
#CellTypeSpecific_Number(RAs,HCs,Rname,Hname,1.0,0.01,0.1)
#OAvsHC
