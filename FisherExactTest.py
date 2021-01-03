#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import scipy.stats


#a=len(ret)
#b=len(tf)-len(ret)
#c=len(retNo)
#d=len(tfno)-len(retNo)
#ChiTest = np.array([[a, b], [c, d]])
#ChiP=chi2_contingency(ChiTest)[1]
#oddsratio, FisherP = stats.fisher_exact([[a, b], [c, d]])
#hyper=stats.hypergeom.pmf(a, a+b+c+d, a+b, a+c)


# In[15]:


C3与OC分化term
#c是Peak GREAT对应的基因
a=40
b=95-a
c=3201-a
d=122890-b

ChiTest = np.array([[a, b], [c, d]])
ChiP=chi2_contingency(ChiTest)[1]
oddsratio, FisherP = stats.fisher_exact([[a, b], [c, d]])
hyper=stats.hypergeom.pmf(a, a+b+c+d, a+b, a+c)

print ('ChiTest', ChiP)
print ('FisherP', FisherP)
print ('hyper', hyper)


# In[6]:


#C3与调节OC分化term
#c是Peak GREAT对应的基因
a=25
b=63-a
c=3201-a
d=122890-b

ChiTest = np.array([[a, b], [c, d]])
ChiP=chi2_contingency(ChiTest)[1]
oddsratio, FisherP = stats.fisher_exact([[a, b], [c, d]])
hyper=stats.hypergeom.pmf(a, a+b+c+d, a+b, a+c)

print ('ChiTest', ChiP)
print ('FisherP', FisherP)
print ('hyper', hyper)


# In[11]:


#C3与促炎term
a=44
b=148-a
c=3201-a
d=122890-b

ChiTest = np.array([[a, b], [c, d]])
ChiP=chi2_contingency(ChiTest)[1]
oddsratio, FisherP = stats.fisher_exact([[a, b], [c, d]])
hyper=stats.hypergeom.pmf(a, a+b+c+d, a+b, a+c)

print ('ChiTest', ChiP)
print ('FisherP', FisherP)
print ('hyper', hyper)




# In[ ]:





# In[ ]:





# In[ ]:
