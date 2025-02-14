#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import statsmodels.api as sm
import sys
import scipy.stats
import math
import scipy.stats
from scipy.stats.mstats import gmean
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams["font.sans-serif"] = "Arial"
from pylab import savefig
import matplotlib.pyplot as plt
import pandas as pd
from math import pi
from sklearn.linear_model import LinearRegression


###OLS
ClusterN='3'
dataF='clinical_v1.txt'
Data=pd.read_csv(dataF,sep='\t',index_col=0)
X='C{}_score'.format(ClusterN)
#Y=['ESR','CRP','RF','ACCP','TSS','Reads','TCJ','SCJ','Pain','DAS28_ESR4','DAS28_CRP4','DAS28_ESR3','DAS28_CRP3']
#Y=['CRP']
Y=['ESR','CRP','RF','ACCP','TSS','Reads','TCJ','SCJ','GH','DAS28-ESR','DAS28-CRP']

Data.sort_values(by=[X], inplace=True)
for i in Y:
    import statsmodels.api as sm
    DataN=Data.dropna(axis = 0, subset = [i] )

    XX=DataN[X]
    YY=DataN[i]
    XX=sm.add_constant(XX)
    est1 = sm.OLS(YY,XX).fit()

    #拟合Y值
    y_fitted = est1.fittedvalues
    #fig, ax = plt.subplots(figsize=(8,6))
    #ax.plot(DataN[X], DataN[i], 'o', label='data')
    #ax.plot(DataN[X], y_fitted, 'r--.',label='OLS')
    #ax.legend(loc='best')

    #sns.set_style('ticks')

    ax=sns.lmplot(x=X, y=i,ci=95,data=Data,height=2,aspect=1.2,markers='o',line_kws={'color': 'blue'},scatter_kws={"s": 4, "alpha": 1,'color': 'blue'})


    #ax=sns.lmplot(x=X, y=i,ci=95,data=Data,height=2,fit_reg=True,aspect=1.2,markers='o',line_kws={'color': 'black'},scatter_kws={"s": 3, "alpha": 0.55,'color': 'black'})
    plt.title('R2:{},P:{}'.format(est1.rsquared,est1.pvalues[1]), fontsize=9)
    #print (est1.pvalues)
    plt.savefig('Cluster{}_{}.OLS.pdf'.format(ClusterN,i))
    plt.show()

    #print est1.summary()
    #print (est1.rsquared)



###CRP related Peaks:
###step1:
F=Dir+'KmeansCluster_c3.txt'
DF=pd.read_csv(F,sep='\t',index_col=0)

ClinalF=Dir+'KmeansCluster_c3.SumScore.AllPeak.clinicalM.txt'
ClinalDF=pd.read_csv(ClinalF,sep='\t',index_col=0)
ClinalDF

RADF=DF[list(ClinalDF.index)]
CRP=list(ClinalDF['CRP'])

import statsmodels.api as sm
def Pearson_R(x):
    x=sm.add_constant(x)
    return (sm.OLS(CRP,x).fit()).rsquared
def Pearson_P(x):
    x=sm.add_constant(x)
    return (sm.OLS(CRP,x).fit()).pvalues[1]

#XX=sm.add_constant(XX)
#sm.add_constant(x)
#est1 = sm.OLS(YY,XX).fit()

DF_R=RADF.copy()
DF_R['Rsquare']=DF_R.apply(Pearson_R,axis=1)
DF_R.to_csv(Dir+'RA_Rsquare.txt',sep='\t')

DF_P=RADF.copy()
DF_P['Pvalue']=DF_P.apply(Pearson_P,axis=1)
DF_P.to_csv(Dir+'RA_Pvalue.txt',sep='\t')


###step2: plot
PlotF=Dir+'RApateint_OLS.txt'
Data=ReadTable(PlotF)
Data_SortRank = Data.sort_values(by=['Rsquare'],ascending=True)
Data_SortRank.head()

sns.set_style('ticks')
fig=plt.figure(figsize=(1,3))
plt.scatter(Data_SortRank['Rsquare'],Data_SortRank['logP'], s=1,marker='o',linestyle='-',c=Data_SortRank['Rsquare'],cmap='viridis_r', alpha=1)
plt.axhline(y=1.302954452,ls='--',color='grey',alpha=0.5,lw=0.5)
plt.axvline(x=0.306234788,ls='--',color='grey',alpha=0.5,lw=0.5)
plt.tick_params(axis='y',labelsize=10,direction='out')
plt.tick_params(axis='x',labelsize=10,direction='out')
#plt.title('UniversityRank vs Score', fontsize=16)
plt.xlabel('R square', fontsize=12)
plt.ylabel('-log10(p value)', fontsize=12)
#plt.text(-49,72, 'univeristy1' , ha='center', va= 'bottom',fontsize=15)
#plt.plot(-49,72,'ro',color='darkred')
#plt.text(-11,89, 'univeristy2' , ha='center', va= 'bottom',fontsize=15)
#plt.plot(-11,89,'ro',color='darkred')
ax=plt.gca()
ax.spines['bottom'].set_linewidth(0.25)
ax.spines['left'].set_linewidth(0.25)
ax.spines['right'].set_linewidth(0.0)
ax.spines['top'].set_linewidth(0.0)
#ax.set_xlim([0.3,0.9])
#ax.set_ylim([1.5,5.5])
fig.savefig(Dir+'PeaksPerson_inC3.pdf')
plt.show()


#
