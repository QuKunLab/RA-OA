#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import statsmodels.api as sm
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import os
import sys
import scipy.stats
from scipy.stats.mstats import gmean
import scipy.stats as stats
import math
import matplotlib as mpl
from sklearn.cluster import KMeans
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams["font.sans-serif"] = "Arial"


#1.Z-score Normalzie DiseaseSP_DF:
Cell='Monocytes'
outDir=os.path.join('{}/DiffPeaks/mean3_fc2_p0.001_fdr0.05/KmeansCluster'.format(Cell))
if not os.path.exists(outDir):
    os.mkdir(outDir)

DiseaseSP_F='{}/DiffPeaks/mean3_fc2_p0.001_fdr0.05/TwoTwoCompare_Merge.sortCol.txt'.format(Cell)
DiseaseSP_DF=pd.read_table(DiseaseSP_F,sep='\t',index_col=0)
DiseaseSP_DFz= DiseaseSP_DF.apply(scipy.stats.zscore,axis=1,result_type='broadcast')

#decide K:1.手肘法(误差平方法SSE);2.轮廓系数法
SSE = []  # 存放每次结果的误差平方和
for k in range(1,10):
    estimator = KMeans(n_clusters=k)
    estimator.fit(DiseaseSP_DFz)
    SSE.append(estimator.inertia_)
X = range(1,10)
plt.style.use('seaborn-white')
fig=plt.figure(figsize=(3.5,2))
ax=fig.add_axes([0.2,0.2,0.7,0.7])
ax.set_ylabel('Sum of the squared errors',fontsize=10)
ax.set_xlabel('k number',fontsize=10)
ax.tick_params(axis='y',length=7,labelsize=8,direction='out')
ax.tick_params(axis='x',length=7,labelsize=8,direction='out')
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
ax.spines['right'].set_linewidth(0.5)
ax.spines['top'].set_linewidth(0.5)
plt.plot(X,SSE,color='purple', marker='o', linestyle='dashed',linewidth=1, markersize=5)
fig.savefig(outDir+'/Kvalue_SSE.pdf')
#print '误差平方和：'
plt.show()

2.#根据最佳K值进行KMeans聚类 (Kmeans聚类用的ZscoreNorm后的DF!!!)
def KMean_Cluster(DF,outDirPrefix,k):
    #print 'Do KMean Cluster, k={}'.format(k)
    kmeans=KMeans(n_clusters=k)
    kmeans.fit(DF)
    Kcluster=pd.DataFrame(kmeans.labels_,index=list(DF.index),columns=['Cluster'])
    Kcluster.to_csv(outDir+'/TwoTwoCompareMerge_zscore_k{}.txt'.format(k),sep='\t')
    #return Kcluster
KMean_Cluster(DiseaseSP_DFz,outDir,2)
KMean_Cluster(DiseaseSP_DFz,outDir,3)

print ('K-means Done !')


# In[5]:


k='3'
Cell='Monocytes'
DiseaseSP_F='{}/DiffPeaks/mean3_fc2_p0.001_fdr0.05/TwoTwoCompare_Merge.sortCol.txt'.format(Cell)
DiseaseSP_DF=pd.read_table(DiseaseSP_F,sep='\t',index_col=0)
RAs=[i for i in list(DiseaseSP_DF) if 'RA' in i]
OAs=[i for i in list(DiseaseSP_DF) if 'OA' in i]
HCs=[i for i in list(DiseaseSP_DF) if 'HC' in i]

BedF= '{}/RAOAHC.removeY.bed'.format(Cell)  #read PeakBed
BedDF=pd.read_table(BedF,sep='\t',header=None)
BedDF.index=BedDF[3]

def PlotKmeanCluster_K3(k):
    kmeansDir=os.path.join('{}/DiffPeaks/mean3_fc2_p0.001_fdr0.05/KmeansCluster/kvalue_k{}/'.format(Cell,k))
    if not os.path.exists(kmeansDir):
        os.mkdir(kmeansDir)


    KClusterF='{}/DiffPeaks/mean3_fc2_p0.001_fdr0.05/KmeansCluster/TwoTwoCompareMerge_zscore_k{}.txt'.format(Cell,k)
    KCluster=pd.read_table(KClusterF,sep='\t',index_col=0)

    k1=KCluster[KCluster['Cluster']==0]
    k2=KCluster[KCluster['Cluster']==1]
    k3=KCluster[KCluster['Cluster']==2]
    k1DF=DiseaseSP_DF.loc[k1.index]
    k2DF=DiseaseSP_DF.loc[k2.index]
    k3DF=DiseaseSP_DF.loc[k3.index]
    k1Bed=BedDF.loc[k1DF.index]
    k2Bed=BedDF.loc[k2DF.index]
    k3Bed=BedDF.loc[k3DF.index]

    a1=k1DF.iloc[:,-2:-1].mean(axis=0)[0]
    a2=k2DF.iloc[:,-2:-1].mean(axis=0)[0]
    a3=k3DF.iloc[:,-2:-1].mean(axis=0)[0]
    if (a1 < a2) & (a2 < a3):
        KclusterDF_c1=k1DF.copy()
        KclusterDF_c2=k2DF.copy()
        KclusterDF_c3=k3DF.copy()
    elif (a1 < a3) & (a3 < a2):
        KclusterDF_c1=k1DF.copy()
        KclusterDF_c2=k3DF.copy()
        KclusterDF_c3=k2DF.copy()
    elif (a2 < a1) & (a1 < a3):
        KclusterDF_c1=k2DF.copy()
        KclusterDF_c2=k1DF.copy()
        KclusterDF_c3=k3DF.copy()
    elif (a2 < a3) & (a3 < a1):
        KclusterDF_c1=k2DF.copy()
        KclusterDF_c2=k3DF.copy()
        KclusterDF_c3=k1DF.copy()
    elif (a3 < a1) & (a1 < a2):
        KclusterDF_c1=k3DF.copy()
        KclusterDF_c2=k1DF.copy()
        KclusterDF_c3=k2DF.copy()
    elif (a3 < a2) & (a2 < a1):
        KclusterDF_c1=k3DF.copy()
        KclusterDF_c2=k2DF.copy()
        KclusterDF_c3=k1DF.copy()

    KclusterBed_c1=BedDF.loc[KclusterDF_c1.index]
    KclusterBed_c2=BedDF.loc[KclusterDF_c2.index]
    KclusterBed_c3=BedDF.loc[KclusterDF_c3.index]
    KclusterBed_c1.to_csv(kmeansDir+'KmeansCluster_c1.bed',sep='\t',header=False,index=False)
    KclusterBed_c2.to_csv(kmeansDir+'KmeansCluster_c2.bed',sep='\t',header=False,index=False)
    KclusterBed_c3.to_csv(kmeansDir+'KmeansCluster_c3.bed',sep='\t',header=False,index=False)
    KclusterDF_c1.to_csv(kmeansDir+'KmeansCluster_c1.txt',sep='\t')
    KclusterDF_c2.to_csv(kmeansDir+'KmeansCluster_c2.txt',sep='\t')
    KclusterDF_c3.to_csv(kmeansDir+'KmeansCluster_c3.txt',sep='\t')
    KclusterDF_c1c2c3=pd.concat([KclusterDF_c1,KclusterDF_c2,KclusterDF_c3],axis=0)
    KclusterDF_c1c2c3.to_csv(kmeansDir+'KmeansCluster_all.txt',sep='\t')
    KclusterBed_c1c2c3=BedDF.loc[KclusterDF_c1c2c3.index]
    KclusterBed_c1c2c3.to_csv(kmeansDir+'KmeansCluster_all.bed',sep='\t',header=False,index=False)

    def DFmean(inputDF,C):
        Df=DiseaseSP_DF.loc[inputDF.index]
        hc=Df[HCs]
        oa=Df[OAs]
        ra=Df[RAs]
        hcmean=hc.mean(axis=1)
        hcmeanDF = hcmean.to_frame()
        hcmeanDF.rename(columns={0:'HC'}, inplace = True)
        oamean=oa.mean(axis=1)
        oameanDF = oamean.to_frame()
        oameanDF.rename(columns={0:'OA'}, inplace = True)
        ramean=ra.mean(axis=1)
        rameanDF = ramean.to_frame()
        rameanDF.rename(columns={0:'RA'}, inplace = True)
        MergeM = pd.concat([hcmeanDF,oameanDF,rameanDF],axis=1)
        MergeM.to_csv(kmeansDir+'KmeansCluster_{}.average.txt'.format(C),sep='\t')
        #Boxplot
        plt.style.use('seaborn-white')
        fig=plt.figure(figsize=(1.5,2))
        ax=fig.add_axes([0.2,0.2,0.75,0.75])
        #sns.violinplot(data=AA,ax=ax1,palette=(['steelblue','gold','orangered']))
        sns.boxplot(data=MergeM,ax=ax,palette=(['steelblue','gold','orangered']),whis=0.5,fliersize=0.5,width=0.7,showfliers=False,medianprops={'linewidth':0.5},whiskerprops={'linewidth':0.5},boxprops={'linewidth':0.5},capprops={'linewidth':0.5})
        ax.tick_params(labelsize=8,width=0.5,direction='out')
        #ax.set_ylim([0,10])
        ax.spines['bottom'].set_linewidth(0.25)
        ax.spines['left'].set_linewidth(0.25)
        ax.spines['right'].set_linewidth(0.25)
        ax.spines['top'].set_linewidth(0.25)
        fig.savefig(kmeansDir+'KmeansCluster_{}_average.boxplot.pdf'.format(C))
        plt.show()
    DFmean(KclusterDF_c1,'c1')
    DFmean(KclusterDF_c2,'c2')
    DFmean(KclusterDF_c3,'c3')

    #zcore，plot heatmap:
    KclusterDFall_Z=KclusterDF_c1c2c3.apply(scipy.stats.zscore,axis=1,result_type='broadcast')
    KclusterDFc1_Z=KclusterDF_c1.apply(scipy.stats.zscore,axis=1,result_type='broadcast')
    KclusterDFc2_Z=KclusterDF_c2.apply(scipy.stats.zscore,axis=1,result_type='broadcast')
    KclusterDFc3_Z=KclusterDF_c3.apply(scipy.stats.zscore,axis=1,result_type='broadcast')

    fig1=sns.clustermap(KclusterDFall_Z,figsize=(4,5),center=0,vmin=-2,vmax=2,col_cluster=False,row_cluster=False,cmap='RdYlBu_r')
    fig1.savefig(kmeansDir+'KmeansCluster_all.heatmap.png',dpi=200)
    plt.show()
    plt.close('all')
    fig2=sns.clustermap(KclusterDFc1_Z,figsize=(4,0.0009*len(KclusterDFc1_Z)),center=0,vmin=-2,vmax=2,col_cluster=False,row_cluster=False,cmap='RdYlBu_r')
    fig2.savefig(kmeansDir+'KmeansCluster_c1.heatmap.png',dpi=500)
    plt.show()
    plt.close('all')
    fig3=sns.clustermap(KclusterDFc2_Z,figsize=(4,0.0009*len(KclusterDFc2_Z)),center=0,vmin=-2,vmax=2,col_cluster=False,row_cluster=False,cmap='RdYlBu_r')
    fig3.savefig(kmeansDir+'KmeansCluster_c2.heatmap.png',dpi=500)
    plt.show()
    plt.close('all')
    fig4=sns.clustermap(KclusterDFc3_Z,figsize=(4,0.0009*len(KclusterDFc3_Z)),center=0,vmin=-2,vmax=2,col_cluster=False,row_cluster=False,cmap='RdYlBu_r')
    fig4.savefig(kmeansDir+'KmeansCluster_c3.heatmap.png',dpi=500)
    plt.show()
    plt.close('all')

    HCz=KclusterDFall_Z[HCs]
    OAz=KclusterDFall_Z[OAs]
    RAz=KclusterDFall_Z[RAs]
    HCmean=HCz.mean(axis=1)
    HCmeanDF = HCmean.to_frame()
    HCmeanDF.rename(columns={0:'HC'}, inplace = True)
    OAmean=OAz.mean(axis=1)
    OAmeanDF = OAmean.to_frame()
    OAmeanDF.rename(columns={0:'OA'}, inplace = True)
    RAmean=RAz.mean(axis=1)
    RAmeanDF = RAmean.to_frame()
    RAmeanDF.rename(columns={0:'RA'}, inplace = True)
    KclusterDFall_Z_average = pd.concat([HCmeanDF,OAmeanDF,RAmeanDF],axis=1)

    fig4=sns.clustermap(KclusterDFall_Z_average,figsize=(1,6),center=0,vmin=-2,vmax=2,col_cluster=False,row_cluster=False,cmap='RdYlBu_r')
    fig4.savefig(kmeansDir+'KmeansCluster_all.heatmap.average.pdf')
    plt.show()
    plt.close('all')


# In[6]:


k='3'
PlotKmeanCluster_K3(k)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:
