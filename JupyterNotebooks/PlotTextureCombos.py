#!/usr/bin/env python
# coding: utf-8

# # PoleFigurePhaseFractions.ipynb
# Written by Adam Creuziger (adam.creuziger@nist.gov)
# 
# Oct 2017
# 
#     This data was developed by employees of the National Institute of Standards and Technology (NIST), an agency of the Federal Government. Pursuant to title 17 United States Code Section 105, works of NIST employees are not subject to copyright protection in the United States and are considered to be in the public domain.
# 
#     The data is provided by NIST as a public service and is expressly provided "AS IS." NIST MAKES NO WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR STATUTORY, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT AND DATA ACCURACY. NIST does not warrant or make any representations regarding the use of the data or the results thereof, including but not limited to the correctness, accuracy, reliability or usefulness of the data. NIST SHALL NOT BE LIABLE AND YOU HEREBY RELEASE NIST FROM LIABILITY FOR ANY INDIRECT, CONSEQUENTIAL, SPECIAL, OR INCIDENTAL DAMAGES (INCLUDING DAMAGES FOR LOSS OF BUSINESS PROFITS, BUSINESS INTERRUPTION, LOSS OF BUSINESS INFORMATION, AND THE LIKE), WHETHER ARISING IN TORT, CONTRACT, OR OTHERWISE, ARISING FROM OR RELATING TO THE DATA (OR THE USE OF OR INABILITY TO USE THIS DATA), EVEN IF NIST HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
# 
#     To the extent that NIST may hold copyright in countries other than the United States, you are hereby granted the non-exclusive irrevocable and unconditional right to print, publish, prepare derivative works and distribute the NIST data, in any medium, or authorize others to do so on your behalf, on a royalty-free basis throughout the world.
# 
#     You may improve, modify, and create derivative works of the data or any portion of the data, and you may copy and distribute such modifications or works. Modified works should carry a notice stating that you changed the data and should note the date and nature of any such change. Please explicitly acknowledge the National Institute of Standards and Technology as the source of the data: Data citation recommendations are provided below.
# 
#     Permission to use this data is contingent upon your acceptance of the terms of this agreement and upon your providing appropriate acknowledgments of NIST's creation of the data.
# 
# 
# See: https://www.nist.gov/director/licensing

# In[28]:


import numpy as np
import pandas as pd
import scipy as scipy
from scipy import interpolate
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import mplstereonet
import math
import os
import glob
import statistics


# In[29]:


# Get the current working directory path
cwd=os.getcwd()
#print cwd
xpcdatapath=os.path.abspath(os.path.join(os.path.dirname( cwd)))
#print xpcdatapath
#Folder=cwd+'/AveragedIntensites'
Folder='C:\Research\Texture-Sampling-PhaseMeasurement-BiasErrors-master\JupyterNotebooks\AveragedIntensites'

print (Folder)

#DFA=pd.read_excel((os.path.join(Folder, "alpha1F-HW20.xlsx"))


# ## Plot the effect of sample orientation and number of peaks (Fig 7)

# In[30]:



#read in the excel files with the averaged intensities
DFF=pd.read_excel((os.path.join(Folder,"TRIP700F.xlsx")),header=1,skipfooter=0)
DFA=pd.read_excel((os.path.join(Folder,"TRIP700A.xlsx")),header=1,skipfooter=0)

#Switch position of Tilt -row 7 and Rotate -row 8 data rows, matches Figure 4 better
#http://stackoverflow.com/questions/32929927/pandas-swap-rows-between-dataframes
tempF=DFF.loc[8]
DFF.loc[8,:]=DFF.loc[7,:].values
DFF.loc[7,:]=tempF.values

tempA=DFA.loc[8]
DFA.loc[8,:]=DFA.loc[7,:].values
DFA.loc[7,:]=tempA.values

#### Changed all .loc in node to .iloc   #unchanged
VF=.25
tolerr=0.05

#  Changed DF1-3 to DF1f-3f to do a temporary fix
DF1f=(VF*DFA["2Pairs"]/(VF*DFA["2Pairs"]+((1.0-VF)*DFF["2Pairs"])))-VF
DF2f=(VF*DFA["4Pairs"]/(VF*DFA["4Pairs"]+((1.0-VF)*DFF["4Pairs"])))-VF
DF3f=(VF*DFA["MaxUnique"]/(VF*DFA["MaxUnique"]+((1.0-VF)*DFF["MaxUnique"])))-VF
##nparrange(13)
ind = np.arange(5)  # the x locations for the groups
width = 0.15       # the width of the bars


#added as check
DF1f=DF1f.dropna()
DF2f=DF2f.dropna()
DF3f=DF3f.dropna()

DF1f=DF1f.drop([1,2,3,4,5,6,7,8,9])
DF2f=DF2f.drop([1,2,3,4,5,6,7,8,9])
DF3f=DF3f.drop([1,2,3,4,5,6,7,8,9])
DF1=DF1f.values
DF2=DF2f.values
DF3=DF3f.values
#print (DF1)
# print (ind)

# end check

## Actually start plotting

fig = plt.figure(figsize=(8,8), dpi=600)

ax1 = fig.add_subplot(211)

#fig, ax = plt.subplots()
ax1.axhline(y=VF*(1+tolerr),color='.4',ls="--",zorder=1)
ax1.axhline(y=VF*(1-tolerr),color='.4',ls="--",zorder=1)

rects1 = ax1.bar(ind, DF1, width, bottom=.25, facecolor='r', edgecolor='k',zorder=4)
rects2 = ax1.bar(ind+width, DF2, width,bottom=.25, facecolor='y', edgecolor='k',zorder=4)
rects3 = ax1.bar(ind+2*width, DF3, width,bottom=.25, facecolor='b', edgecolor='k',zorder=4)

#tried as stem plot, but this doesn't work as well..
#rects1 = ax.stem(ind, DF1, linefmt='r-',markerfmt='rs')
#rects2 = ax.stem(ind+width,DF2, linefmt='y-',markerfmt='rs')
#rects3 = ax.stem(ind+2*width,DF3, linefmt='b-',markerfmt='rs')


ax1.plot([0, 12], [VF, VF], color='k', linestyle='-', linewidth=.25)

ax1.axvline(x=3.75,color='k')
ax1.axvline(x=6.75,color='k')
ax1.axvline(x=9.75,color='k')



# add some text for labels, title and axes ticks
ax1.set_ylabel('Austenite Phase Fraction',fontsize=16)
#ax1.set_title('Effect of Sampling on Measured Phase Fraction - TRIP 700')
ax1.set_xticks(ind + width)
ax1.xaxis.set_tick_params(labelsize=16)
ax1.yaxis.set_tick_params(labelsize=16)
# Fix the order
#ax.set_xticklabels(DFF['HKL'], rotation=90)
ax1.set_xticklabels(["","","",""])  ##MCchange, 4 pairs of "" from 13
#ax1.set_xticklabels(["ND (Single)","<- RD ->","TD","Morris","ND","RD","TD","Rotation","60$^\circ$ Tilt",
#                    "Tilt & Rotation","5$^\circ$ Full",
#                    "22.5$^\circ$ Full","5$^\circ$ Partial",], rotation=90, va='center', y=-0.21)
#ax1.set_xticklabels(["ND","15$^\circ$ Partial","30$^\circ$ Partial","20.5$^\circ$ Partial",], rotation=90)


#fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)
fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.1)


bbox_props = dict(boxstyle="square,pad=0.3", fc="white", ec="k", lw=1)
ax1.text(1.25, 0.35, "TRIP 700", bbox=bbox_props,  ha='center', fontsize=18)
#11.25,0.35
# ax1.text(1.05, 0.12, "Single", bbox=bbox_props,  ha='center', fontsize=14)
# ax1.text(5.25, 0.12, "Ring", bbox=bbox_props,  ha='center', fontsize=14) 
# ax1.text(11.25, 0.12, "Hex", bbox=bbox_props,  ha='center', fontsize=14) 
# ax1.text(8.25, 0.12, "Tilt &\nRotate", bbox=bbox_props,  ha='center', fontsize=14) 


##MCchange, adjusted to 3 non-overlapping, non-duplicate pairs
##ax1.legend((rects1[0], rects2[0], rects3[0]), ('2Pairs', '4Pairs', 'MaxUnique')
ax1.legend((rects1[0], rects2[0], rects3[0]), ('2Pairs', '3Pairs', 'MaxUnique'), 
          bbox_to_anchor=(0., 1.12, 1., .102), loc=3,
           ncol=3, mode="expand", borderaxespad=0.,fontsize=18)

ax1.set_ylim([VF-.15,VF+.15])
ax1.set_xlim([-.25, 12.75])

ax1.set_axisbelow(True)
plt.grid()

#plt.tight_layout()


#### Now for the second file

DFF=pd.read_excel((os.path.join(Folder,"TRIP780F.xlsx")),header=1,skipfooter=0)
DFA=pd.read_excel((os.path.join(Folder,"TRIP780A.xlsx")),header=1,skipfooter=0)

## Switch 


##  Changed DF1-3 to DF1f-3f to do a temporary fix
DF1f1=(VF*DFA["2Pairs"]/(VF*DFA["2Pairs"]+((1.0-VF)*DFF["2Pairs"])))-VF
DF2f1=(VF*DFA["4Pairs"]/(VF*DFA["4Pairs"]+((1.0-VF)*DFF["4Pairs"])))-VF
DF3f1=(VF*DFA["MaxUnique"]/(VF*DFA["MaxUnique"]+((1.0-VF)*DFF["MaxUnique"])))-VF

ind = np.arange(5)  # the x locations for the groups
width = 0.15       # the width of the bars


#added as check
DF1f1=DF1f1.dropna()
DF2f1=DF2f1.dropna()
DF3f1=DF3f1.dropna()

DF1f1=DF1f1.drop([1,2,3,4,5,6,7,8,9])
DF2f1=DF2f1.drop([1,2,3,4,5,6,7,8,9])
DF3f1=DF3f1.drop([1,2,3,4,5,6,7,8,9])
DF11=DF1f1.values
DF21=DF2f1.values
DF31=DF3f1.values
# print (DF1)
# print (ind)

## end check

# DF1=(VF*DFA["2Pairs"]/(VF*DFA["2Pairs"]+((1.0-VF)*DFF["2Pairs"])))-VF
# DF2=(VF*DFA["4Pairs"]/(VF*DFA["4Pairs"]+((1.0-VF)*DFF["4Pairs"])))-VF
# DF3=(VF*DFA["MaxUnique"]/(VF*DFA["MaxUnique"]+((1.0-VF)*DFF["MaxUnique"])))-VF

# ind = np.arange(13)  # the x locations for the groups
# width = 0.15       # the width of the bars

ax2 = fig.add_subplot(212)

ax2.axhline(y=VF*(1+tolerr),color='.4',ls="--",zorder=1)
ax2.axhline(y=VF*(1-tolerr),color='.4',ls="--",zorder=1)

rects1 = ax2.bar(ind, DF11, width, bottom=.25, facecolor='r', edgecolor='k',zorder=4)
rects2 = ax2.bar(ind+width, DF21, width,bottom=.25,  facecolor='y', edgecolor='k',zorder=4)
rects3 = ax2.bar(ind+2*width, DF31, width,bottom=.25,  facecolor='b', edgecolor='k',zorder=4)

ax2.plot([0, 12], [VF, VF], color='k', linestyle='-', linewidth=.25)

ax2.axvline(x=3.75,color='k')
ax2.axvline(x=6.75,color='k')
ax2.axvline(x=9.75,color='k')

# add some text for labels, title and axes ticks
ax2.set_ylabel('Austenite Phase Fraction', fontsize=16)
#ax2.set_title('Effect of Sampling on Measured Phase Fraction - TRIP 780')
ax2.set_xticks(ind + width)
ax2.xaxis.set_tick_params(labelsize=16)
ax2.yaxis.set_tick_params(labelsize=16)

#ax.set_xticklabels(DFF['HKL'], rotation=90)
#ax2.set_xticklabels(["","","","","","","","","","","","",""])
# ax2.set_xticklabels(["ND","RD","TD","Morris","ND","RD","TD","Rotation","60$^\circ$ Tilt",
#                      "Tilt & Rotation","15$^\circ$ Partial",
#                    "30$^\circ$ Partial","20.5$^\circ$ Partial",], rotation=90)
ax2.set_xticklabels(["ND","5$^\circ$ Partial","30$^\circ$ Partial","20.5$^\circ$ Partial",], rotation=90)

bbox_props = dict(boxstyle="square,pad=0.3", fc="white", ec="k", lw=1)
ax2.text(11.25, 0.35, "TRIP 780", bbox=bbox_props,  ha='center', fontsize=18)
ax2.text(1.05, 0.12, "Single", bbox=bbox_props,  ha='center', fontsize=16)
ax2.text(5.25, 0.12, "Ring", bbox=bbox_props,  ha='center', fontsize=16) 
ax2.text(11.25, 0.12, "Hex", bbox=bbox_props,  ha='center', fontsize=16) 
ax2.text(8.25, 0.12, "Tilt &\nRotate", bbox=bbox_props,  ha='center', fontsize=16) 


#ax2.legend((rects1[0], rects2[0], rects3[0]), ('2Pairs', '4Pairs', 'MaxUnique'), 
#          bbox_to_anchor=(0., 1.12, 1., .102), loc=3,
#           ncol=3, mode="expand", borderaxespad=0.)

ax2.set_ylim([VF-.15,VF+.15])
ax2.set_xlim([-.25, 12.75])

ax2.set_axisbelow(True)
plt.grid()
plt.savefig("ExampleTex-Draft1.png",dpi=600,format="png",orientation='portrait',bbox_inches='tight')

##New additions: Duplex and TRIP780b begin here ########

#TRIP780b start

DFF=pd.read_excel((os.path.join(Folder,"TRIP7802F.xlsx")),header=1,skipfooter=0)
DFA=pd.read_excel((os.path.join(Folder,"TRIP7802A.xlsx")),header=1,skipfooter=0)

#Switch position of Tilt -row 7 and Rotate -row 8 data rows, matches Figure 4 better
#http://stackoverflow.com/questions/32929927/pandas-swap-rows-between-dataframes
tempF=DFF.loc[8]
DFF.loc[8,:]=DFF.loc[7,:].values
DFF.loc[7,:]=tempF.values

tempA=DFA.loc[8]
DFA.loc[8,:]=DFA.loc[7,:].values
DFA.loc[7,:]=tempA.values


VF=.25
tolerr=0.05


#  Changed DF1-3 to DF1f-3f to do a temporary fix, added subscripts to attempt single plot
DF1f2=(VF*DFA["2Pairs"]/(VF*DFA["2Pairs"]+((1.0-VF)*DFF["2Pairs"])))-VF
DF2f2=(VF*DFA["4Pairs"]/(VF*DFA["4Pairs"]+((1.0-VF)*DFF["4Pairs"])))-VF
DF3f2=(VF*DFA["MaxUnique"]/(VF*DFA["MaxUnique"]+((1.0-VF)*DFF["MaxUnique"])))-VF

ind = np.arange(5)  # the x locations for the groups
width = 0.15       # the width of the bars


#added as check
DF1f2=DF1f2.dropna()
DF2f2=DF2f2.dropna()
DF3f2=DF3f2.dropna()

DF1f2=DF1f2.drop([1,2,3,4,5,6,7,8,9])
DF2f2=DF2f2.drop([1,2,3,4,5,6,7,8,9])
DF3f2=DF3f2.drop([1,2,3,4,5,6,7,8,9])
DF12=DF1f2.values
DF22=DF2f2.values
DF32=DF3f2.values
#print (DF1)
#print (ind)

# end check

# Actually start plotting

fig = plt.figure(figsize=(8,8), dpi=600)
##Commented below to produce a single plot
# ax3 = fig.add_subplot(211)

# #fig, ax = plt.subplots()
# ax3.axhline(y=VF*(1+tolerr),color='.4',ls="--",zorder=1)
# ax3.axhline(y=VF*(1-tolerr),color='.4',ls="--",zorder=1)

# rects1 = ax3.bar(ind, DF12, width, bottom=.25, facecolor='r', edgecolor='k',zorder=4)
# rects2 = ax3.bar(ind+width, DF22, width,bottom=.25, facecolor='y', edgecolor='k',zorder=4)
# rects3 = ax3.bar(ind+2*width, DF32, width,bottom=.25, facecolor='b', edgecolor='k',zorder=4)

# #tried as stem plot, but this doesn't work as well..
# #rects1 = ax.stem(ind, DF1, linefmt='r-',markerfmt='rs')
# #rects2 = ax.stem(ind+width,DF2, linefmt='y-',markerfmt='rs')
# #rects3 = ax.stem(ind+2*width,DF3, linefmt='b-',markerfmt='rs')


# ax3.plot([0, 12], [VF, VF], color='k', linestyle='-', linewidth=.25)

# ax3.axvline(x=3.75,color='k')
# ax3.axvline(x=6.75,color='k')
# ax3.axvline(x=9.75,color='k')

# # add some text for labels, title and axes ticks
# ax3.set_ylabel('Austenite Phase Fraction',fontsize=16)
# #ax1.set_title('Effect of Sampling on Measured Phase Fraction - TRIP 700')
# ax3.set_xticks(ind + width)
# ax3.xaxis.set_tick_params(labelsize=16)
# ax3.yaxis.set_tick_params(labelsize=16)
# # Fix the order
# #ax.set_xticklabels(DFF['HKL'], rotation=90)
# ax3.set_xticklabels(["","","","","","","","","","","","",""])
# #ax1.set_xticklabels(["ND (Single)","<- RD ->","TD","Morris","ND","RD","TD","Rotation","60$^\circ$ Tilt",
# #                    "Tilt & Rotation","5$^\circ$ Full",
# #                    "22.5$^\circ$ Full","5$^\circ$ Partial",], rotation=90, va='center', y=-0.21)

# #fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)
# fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.1)


# bbox_props = dict(boxstyle="square,pad=0.3", fc="white", ec="k", lw=1)
# ax3.text(11.25, 0.35, "TRIP 780b", bbox=bbox_props,  ha='center', fontsize=18)
# ax3.text(1.05, 0.12, "Single", bbox=bbox_props,  ha='center', fontsize=14)
# ax3.text(5.25, 0.12, "Ring", bbox=bbox_props,  ha='center', fontsize=14) 
# ax3.text(11.25, 0.12, "Hex", bbox=bbox_props,  ha='center', fontsize=14) 
# ax3.text(8.25, 0.12, "Tilt &\nRotate", bbox=bbox_props,  ha='center', fontsize=14) 

##MCchange, adjusted to not duplicate/3pairs
##ax3.legend((rects1[0], rects2[0], rects3[0]), ('2Pairs', '4Pairs', 'MaxUnique')
# ax3.legend((rects1[0], rects2[0], rects3[0]), ('2Pairs', '3Pairs', 'MaxUnique'), 
#           bbox_to_anchor=(0., 1.12, 1., .102), loc=3,
#            ncol=3, mode="expand", borderaxespad=0.,fontsize=18)

# ax3.set_ylim([VF-.15,VF+.15])
# ax3.set_xlim([-.25, 15])

# ax3.set_axisbelow(True)
# plt.grid()

#TRIP 780b end
#duplex starts

DFF=pd.read_excel((os.path.join(Folder,"DuplexF.xlsx")),header=1,skipfooter=0)
DFA=pd.read_excel((os.path.join(Folder,"DuplexA.xlsx")),header=1,skipfooter=0)

## Switch 


#  Changed DF1-3 to DF1f-3f to do a temporary fix
DF1f3=(VF*DFA["2Pairs"]/(VF*DFA["2Pairs"]+((1.0-VF)*DFF["2Pairs"])))-VF
DF2f3=(VF*DFA["4Pairs"]/(VF*DFA["4Pairs"]+((1.0-VF)*DFF["4Pairs"])))-VF
DF3f3=(VF*DFA["MaxUnique"]/(VF*DFA["MaxUnique"]+((1.0-VF)*DFF["MaxUnique"])))-VF
#MCchange ind 16 to 4
ind = np.arange(5)  # the x locations for the groups
ind2=np.vstack((ind,ind,ind,ind,ind))
width = 0.8       # the width of the bars


#added as check
DF1f3=DF1f3.dropna()
DF2f3=DF2f3.dropna()
DF3f3=DF3f3.dropna()

DF1f3=DF1f3.drop([1,2,3,4,5,6,7,8,9])
DF2f3=DF2f3.drop([1,2,3,4,5,6,7,8,9])
DF3f3=DF3f3.drop([1,2,3,4,5,6,7,8,9])
DF13=DF1f3.values
DF23=DF2f3.values
DF33=DF3f3.values

##MCchanges to produce ordered scatter (least to most coarse)
DF2=np.append(DF2,DF2[0])
DF21=np.append(DF21,DF21[0])
DF22=np.append(DF22,DF22[0])
DF23=np.append(DF23,DF23[0])

DF2=np.delete(DF2,[0])
DF21=np.delete(DF21,[0])
DF22=np.delete(DF22,[0])
DF23=np.delete(DF23,[0])

#print (DF1)
#print (ind)
# DF1 = np.concatenate([DF1,DF11,DF12,DF13], axis=0)
# DF2 = np.vstack((DF2,DF21,DF22,DF23))
# DF3 = np.concatenate([DF3,DF31,DF32,DF33], axis=0)
# end check
print(DF2)
ax4 = fig.add_subplot(212)

ax4.axhline(y=VF*(1+tolerr),color='.4',ls="--",zorder=1)
ax4.axhline(y=VF,color='k',ls="-",zorder=1)
ax4.axhline(y=VF*(1-tolerr),color='.4',ls="--",zorder=1)
##MCchange, commented out 2+4 pairs, Changed to scatter?
#rects1 = ax4.bar(ind, DF1, width, bottom=.25, facecolor='r', edgecolor='k',zorder=4)
#rects2 = ax4.bar(ind+width, DF2, width,bottom=.25,  facecolor='y', edgecolor='k',zorder=4)
#rects3 = ax4.bar(ind+2*width, DF3, width,bottom=.25,  facecolor='b', edgecolor='k',zorder=4)
scatter1=ax4.scatter(ind+width-.12, DF2+.25, color='g', label="TRIP 700", s=40, alpha=0.80, marker='o',zorder=4)
scatter2=ax4.scatter(ind+width-.04, DF21+.25, color='b', label="TRIP 780", s=40, alpha=0.80, marker='o',zorder=4)
scatter3=ax4.scatter(ind+width+.04, DF22+.25, color='r', label="TRIP 780b", s=40, alpha=0.80, marker='o',zorder=4)
scatter4=ax4.scatter(ind+width+.12, DF23+.25, color='#6a0dad', label="Duplex", s=40, alpha=0.80, marker='o',zorder=4)
# scatter2=ax.scatter(np.full(elem.shape, i) ,  list(df4pair.loc[i])[1:] ,
#                             color='y',label="3Pairs", s=18, alpha=0.33, marker='o',zorder=4)


ax4.plot([0, 16], [VF, VF], color='k', linestyle='-', linewidth=.25)
##MChange to one plot
# ax4.axvline(x=3.7,color='k')
# ax4.axvline(x=7.55,color='k')
# ax4.axvline(x=11.75,color='k')

# add some text for labels, title and axes ticks
ax4.set_ylabel('Austenite Phase Fraction', fontsize=16)
#ax2.set_title('Effect of Sampling on Measured Phase Fraction - TRIP 780')
ax4.set_xticks(ind + width)
ax4.xaxis.set_tick_params(labelsize=16)
ax4.yaxis.set_tick_params(labelsize=16)

#ax.set_xticklabels(DFF['HKL'], rotation=90)
#ax2.set_xticklabels(["","","","","","","","","","","","",""])

##MCchange adjusted to 4 bars
#ax4.set_xticklabels(["ND","RD","TD","Morris","ND","RD","TD","Rotation","60$^\circ$ Tilt",
#                     "Tilt & Rotation","15$^\circ$ Partial",
#                    "30$^\circ$ Partial","20.5$^\circ$ Partial",], rotation=90)
# ax4.set_xticklabels(["ND","15$^\circ$ Partial","30$^\circ$ Partial","20.5$^\circ$ Partial",
#                      "ND","15$^\circ$ Partial","30$^\circ$ Partial","20.5$^\circ$ Partial",
#                      "ND","15$^\circ$ Partial","30$^\circ$ Partial","20.5$^\circ$ Partial",
#                      "ND","15$^\circ$ Partial","30$^\circ$ Partial","20.5$^\circ$ Partial"], rotation=90)
ax4.set_xticklabels(["5$^\circ$ Full","15$^\circ$ Full","22.5$^\circ$ Full","30$^\circ$ Full","ND"], rotation=90)



bbox_props = dict(boxstyle="square,pad=0.3", fc="white", ec="k", lw=1)
# ax4.text(11.25, 0.35, "Duplex", bbox=bbox_props,  ha='center', fontsize=18)
# ax4.text(1.05, 0.12, "Single", bbox=bbox_props,  ha='center', fontsize=16)
# ax4.text(5.25, 0.12, "Ring", bbox=bbox_props,  ha='center', fontsize=16) 
# ax4.text(11.25, 0.12, "Hex", bbox=bbox_props,  ha='center', fontsize=16) 
# ax4.text(8.25, 0.12, "Tilt &\nRotate", bbox=bbox_props,  ha='center', fontsize=16) 
# ax4.text(13.5, 0.3, "Duplex", bbox=bbox_props,  ha='center', fontsize=18)
# ax4.text(1.6, 0.3, "TRIP 700", bbox=bbox_props,  ha='center', fontsize=16)
# ax4.text(5.6, 0.3, "TRIP 780", bbox=bbox_props,  ha='center', fontsize=16) 
# ax4.text(9.8, 0.3, "TRIP 780b", bbox=bbox_props,  ha='center', fontsize=16) 

#MCchange, dedited out old legend
# ax4.legend((rects1[0], rects2[0], rects3[0]), ('2Pairs', '3Pairs', 'MaxUnique'), 
#           bbox_to_anchor=(0., 1.12, 1., .102), loc=3,
#            ncol=3, mode="expand", borderaxespad=0.,fontsize=18)

ax4.legend((scatter1,scatter2,scatter3,scatter4), ('TRIP 700','TRIP 780','TRIP 780B','Duplex'), 
    bbox_to_anchor=(0., 1.12, .5, .802), loc=3,
    ncol=3, mode="expand", borderaxespad=0.,fontsize=10)



ax4.set_ylim([VF-.125,VF+.10])
ax4.set_xlim([-.25, 5.5])

ax4.set_axisbelow(True)
plt.grid()


#Duplex end
#New additions end

# Remove comments to save files
#plt.savefig("ExampleTex-Draft.pdf",dpi=600,format="pdf",orientation='portrait',bbox_inches='tight')
#plt.savefig("ExampleTex-Draft.eps",dpi=600,format="eps",orientation='portrait',bbox_inches='tight')
plt.savefig("ExampleTex-Draft2.png",dpi=600,format="png",orientation='portrait',bbox_inches='tight')
plt.show()


# # Figure 8: Plot Combinations of different textures

# In[31]:


# Get the current working directory path
cwd=os.getcwd()
#print cwd
xpcdatapath=os.path.abspath(os.path.join(os.path.dirname( cwd)))
#print xpcdatapath
#Folder=cwd+'/AveragedIntensites'
Folder='C:\Research\Texture-Sampling-PhaseMeasurement-BiasErrors-master\JupyterNotebooks\AveragedIntensites'

print (Folder)


# ## Search for HW20 files

# In[32]:


import fnmatch



# Set Volume Fraction
VF=.25

#Read in Texture files matching HW of 20 degrees
AusteniteTextures=[]
FerriteTextures=[]  


for file in os.listdir(Folder):
    #print file
    if fnmatch.fnmatch(file, '*A-HW20*'):
        #print file
        AusteniteTextures.append(file)
    elif fnmatch.fnmatch(file, '*F-HW20*'):
        FerriteTextures.append(file)
    else: ()
        
#print FerriteTextures
#print AusteniteTextures

# create a dataframe shape from existing data
#DFA=pd.read_excel((os.path.join(Folder, "alpha1F-HW20.xlsx")),header=1,skipfooter=0) ##original
DFA=pd.read_excel((os.path.join(Folder, "alpha1F-HW20.xlsx")),header=1,skipfooter=1)  ##fix/test
#copy the HKL reflection indexes
##fix/text .drop([13])
df2pair=DFA["HKL"]
df4pair=DFA["HKL"]
dfMaxUnique=DFA["HKL"]

for AustOrient in AusteniteTextures:
    for FerrOrient in FerriteTextures:
        
##fix/test skipfooter=1, orig 0
        DFF=pd.read_excel(os.path.join(Folder,FerrOrient),header=1,skipfooter=1)
        DFA=pd.read_excel(os.path.join(Folder,AustOrient),header=1,skipfooter=1)

        #Switch position of Tilt -row 7 and Rotate -row 8 data rows, matches Figure 4 better
        #http://stackoverflow.com/questions/32929927/pandas-swap-rows-between-dataframes
        tempF=DFF.loc[8]
        DFF.loc[8,:]=DFF.loc[7,:].values
        DFF.loc[7,:]=tempF.values

        tempA=DFA.loc[8]
        DFA.loc[8,:]=DFA.loc[7,:].values
        DFA.loc[7,:]=tempA.values
        
        #print DFF
        #print DFA
        
        # add minus VF for plotting
        DF1=(VF*DFA["2Pairs"]/(VF*DFA["2Pairs"]+((1.0-VF)*DFF["2Pairs"])))#-VF
        DF2=(VF*DFA["4Pairs"]/(VF*DFA["4Pairs"]+((1.0-VF)*DFF["4Pairs"])))#-VF
        DF3=(VF*DFA["MaxUnique"]/(VF*DFA["MaxUnique"]+((1.0-VF)*DFF["MaxUnique"])))#-VF
        
        Fname=FerrOrient.split(".")
        Aname=AustOrient.split(".")
        
        MixName= Fname[0]+"-"+Aname[0]
        #print MixName
        #print DF1, DF2, DF3
        
        data1 = pd.DataFrame({MixName: DF1})
        data2 = pd.DataFrame({MixName: DF2})
        data3 = pd.DataFrame({MixName: DF3})
        #print data
        
        
        df2pair = pd.concat([df2pair, data1], axis=1)
        df4pair = pd.concat([df4pair, data2], axis=1)
        dfMaxUnique = pd.concat([dfMaxUnique, data3], axis=1)
        
#print df2pair

#TexNames=list(df2pair.columns.values)
#TexNames[1:]
#print list(df[TexNames[0]])
#df1[TexNames[0]]
#print list(df2pair.columns.values)


# ## Plot 

# In[33]:


fig, ax = plt.subplots(figsize=(8,4), dpi=600)

#list of textures
TexNames=list(df2pair.columns.values)

#list of sampling methods
SamplingMethods=list(df2pair[TexNames[0]])

elem=np.arange(len(TexNames[1:]))


#for each sampling method
for i, item in enumerate(df2pair[TexNames[0]]):
    #print i
    #print item
    #Row slice from dataframe, without name
    #print list(df.ix[i])[1:]

    ##MCchange, adjusted to 3 pairs
    #print len(elem),len(list(df1.ix[i])[1:])
    #elem.fill(i)
    if i==0:
        #scatter1=ax.scatter(np.full(elem.shape, i-0.15) ,  list(df2pair.ix[i])[1:] ,
        #                    edgecolors='r',label="2Pair", s=10, alpha=1, marker="o",facecolors='none')
        scatter1=ax.scatter(np.full(elem.shape, i-0.15) ,  list(df2pair.loc[i])[1:] ,
                            color='r',label="2Pairs", s=18, alpha=0.33,marker='o',zorder=4)
        scatter2=ax.scatter(np.full(elem.shape, i) ,  list(df4pair.loc[i])[1:] ,
                            color='y',label="3Pairs", s=18, alpha=0.33, marker='o',zorder=4)
        scatter3=ax.scatter(np.full(elem.shape, i+0.15) ,  list(dfMaxUnique.loc[i])[1:] ,
                            color='b',label="MaxUnique", s=18, alpha=0.33, marker='o',zorder=4)
    else:
        scatter1=ax.scatter(np.full(elem.shape, i-0.15) ,  list(df2pair.loc[i])[1:] ,
                            color='r',label="", s=18, alpha=0.33,marker='o',zorder=4)
        scatter2=ax.scatter(np.full(elem.shape, i) ,  list(df4pair.loc[i])[1:] ,
                            color='y',label="", s=18, alpha=0.33, marker='o',zorder=4)
        scatter3=ax.scatter(np.full(elem.shape, i+0.15) ,  list(dfMaxUnique.loc[i])[1:] ,
                            color='b',label="", s=18, alpha=0.33, marker='o',zorder=4)
        #elem.fill(i)   
    
    
ax.legend(bbox_to_anchor=(0., 1.12, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.,fontsize=18)

ax.xaxis.set_tick_params(labelsize=16)
ax.yaxis.set_tick_params(labelsize=16)

#plt.xticks(np.arange(len(SamplingMethods)),SamplingMethods, rotation=90)

bbox_props = dict(boxstyle="square,pad=0.3", fc="white", ec="k", lw=1)
ax.text(2.75, 0.05, "Single", bbox=bbox_props,  ha='center',fontsize=14)
ax.text(5.0, 0.05, "Ring", bbox=bbox_props,  ha='center',fontsize=14) 
ax.text(11.0, 0.05, "Hex", bbox=bbox_props,  ha='center',fontsize=14) 
ax.text(8.5, 0.025, "Tilt &\nRotate", bbox=bbox_props,  ha='center',fontsize=14) 

plt.xticks(np.arange(len(SamplingMethods)),["ND","RD","TD","Morris","ND","RD","TD",
                                             "Rotation","60$^\circ$ Tilt", "Tilt & Rotation",
                                           "15$^\circ$ Full","30$^\circ$ Full","20.5$^\circ$ Full"], rotation=90,
                                           fontsize=16)

#ax.set_xticklabels(["ND","RD","TD","Morris","ND","RD","TD","5$^\circ$ Full",
                 #   "22.5$^\circ$ Full","5$^\circ$ Partial","60$^\circ$ Tilt", "Rotation", "Tilt & Rotation"], rotation=90)


#ax.plot([-1, 13], [VF, VF], color='k', linestyle='-', linewidth=.5)
ax.axhline(y=VF,color='k',ls="-",zorder=5)
ax.axhline(y=VF*(1+tolerr),color='.4',ls="--",zorder=1)
ax.axhline(y=VF*(1-tolerr),color='.4',ls="--",zorder=1)


ax.axvline(x=3.5,color='k')
ax.axvline(x=6.5,color='k')
ax.axvline(x=9.5,color='k')

ax.set_axisbelow(True)
plt.grid()
ax.set_ylabel('Austenite Phase Fraction', fontsize=16)
#ax.set_title('Effect of Texture Combinations')
ax.set_xlim([-.5, 12.5])

# Remove comments to save files
#plt.savefig("TextureCombinations-Draft.pdf", dpi=600,format="pdf",orientation='portrait',bbox_inches='tight')
#plt.savefig("TextureCombinations-Draft.eps",dpi=600,format="eps",orientation='portrait',bbox_inches='tight')
plt.savefig("TextureCombinations-Draft.png",dpi=600,format="png",orientation='portrait',bbox_inches='tight')
plt.show()


# ## Table 2: Save the number of combinations within the ±5% range
# - Attempted to plot the median and as a histograms, but a table seems to be a better method
# - Plot commands are commented out
# - The median value shows some skew to the bias errors, but this table was not deemed to add sufficently to the paper.

# In[34]:



#VF=0.25
AcceptRange=0.05

#fig, ax = plt.subplots(figsize=(8,4), dpi=600)

#list of textures
TexNames=list(df2pair.columns.values)

#list of sampling methods
SamplingMethods=list(df2pair[TexNames[0]])

elem=np.arange(len(TexNames[1:]))


#sales = [('Nothing', 0, 0, 0)]
##MCchange, adjusted to 3 pairs
labels = ['Scheme', '2peaks', '3peaks', 'MaxUnique']

#print dfMedian
#dfinRange

#for each sampling method
for i, item in enumerate(df2pair[TexNames[0]]):
    #print i
    #print item
    #Row slice from dataframe, without name
    #print list(df.ix[i])[1:]
    
    a=statistics.median(list(df2pair.loc[i])[1:] )
    b=statistics.median(list(df4pair.loc[i])[1:] )
    c=statistics.median(list(dfMaxUnique.loc[i])[1:])
 
    df2inRange=sum(VF*(1-AcceptRange) <= x <= VF*(1+AcceptRange) for x in (df2pair.loc[i])[1:])
    df4inRange=sum(VF*(1-AcceptRange) <= x <= VF*(1+AcceptRange) for x in (df4pair.loc[i])[1:])
    dfMaxUniqueinRange=sum(VF*(1-AcceptRange) <= x <= VF*(1+AcceptRange) for x in (dfMaxUnique.loc[i])[1:])
    
    if i==0:    
        dfMedian=pd.DataFrame.from_records([(item, a, b, c)], columns=labels)
        dfinRange=pd.DataFrame.from_records([(item, df2inRange, df4inRange, dfMaxUniqueinRange)], columns=labels)
    else:  
        #print list([item, a, b, c])
        dfMedianI=pd.DataFrame.from_records([(item, a, b, c)], columns=labels)
        dfinRangeI=pd.DataFrame.from_records([(item, df2inRange, df4inRange, dfMaxUniqueinRange)], columns=labels)
        #dfMedian = pd.DataFrame.from_records(sales, columns=labels)    
        dfMedian=dfMedian.append(dfMedianI, ignore_index=True)
        dfinRange=dfinRange.append(dfinRangeI, ignore_index=True)
    
    #print len(list(df2pair.ix[i])[1:])

    #print item, df2inRange, df4inRange, dfMaxUniqueinRange
    
    #print len(elem),len(list(df1.ix[i])[1:])
    #elem.fill(i)

        #scatter1=ax.scatter(np.full(elem.shape, i-0.15) ,  list(df2pair.ix[i])[1:] ,
        #                    edgecolors='r',label="2Pair", s=10, alpha=1, marker="o",facecolors='none')
        #scatter1=ax.hist(  list(df2pair.ix[i])[1:] , bins,
        #                    color='r',label="2Pair",alpha=0.33)
        #scatter2=ax.hist(  list(df4pair.ix[i])[1:] , bins,
        #                    color='y',label="4Pair",alpha=0.33)
        #scatter3=ax.hist(  list(dfMaxUnique.ix[i])[1:] , bins,
        #                    color='b',label="MaxUnique",alpha=0.33)


print (dfMedian)
print (dfinRange)

# Comment out to save tables
#dfMedian.to_csv("MedianTable.csv")
#dfinRange.to_csv("InRangeTable.csv")



#ax.axvline(x=.25,color='k',ls="dashed")
#ax.axvline(x=.25+0.0125,color='c',ls="dashed")
#ax.axvline(x=.25-0.0125,color='c',ls="dashed")
#plt.grid()
#ax.set_ylabel('Number of Cases')
#ax.set_title('Effect of Texture Combinations')
#ax.set_xlim([0, 1])
#ax.set_xlim([.23, .27])

#plt.savefig("TextureCombinations-Draft.pdf", dpi=600,format="pdf",orientation='portrait',bbox_inches='tight')
#plt.savefig("TextureCombinations-Draft.eps",dpi=600,format="eps",orientation='portrait',bbox_inches='tight')
#plt.show()


# ## Plot Effect of Texture Intensity (Fig 9)

# In[35]:


# Get the current working directory path
cwd=os.getcwd()
#print cwd
xpcdatapath=os.path.abspath(os.path.join(os.path.dirname( cwd)))
#print xpcdatapath
#Folder=cwd+'/AveragedIntensites'
Folder='C:\Research\Texture-Sampling-PhaseMeasurement-BiasErrors-master\JupyterNotebooks\AveragedIntensites'

print (Folder)


# ## Retreive texture index file

# In[36]:


# Get the current working directory path
cwd=os.getcwd()
#print cwd
xpcdatapath=os.path.abspath(os.path.join(os.path.dirname( cwd)))
#print xpcdatapath
#TIfile=xpcdatapath+'/Matlab/MtexData/ComputedTextureIndexValues.txt'
TIfile='C:\Research\Texture-Sampling-PhaseMeasurement-BiasErrors-master\Matlab\MtexData\ComputedTextureIndexValues.txt'

print (TIfile)


# ### Diagnostics to check file status

# In[37]:


#read in TIfile

TextureIndexes=pd.read_table(TIfile)
#print TextureIndexes

#strip whitespace
TextureIndexes.columns=TextureIndexes.columns.str.strip()

#TextureIndexes.loc[TextureIndexes['Name'] == 'Uniform']
#TI=TextureIndexes[TextureIndexes['Name'].str.contains("SA-HW10")].as_matrix() #The .as_matrix method will be removed soon, .values replaces; alternative is interactiveshellapp.init_path()
TI=TextureIndexes[TextureIndexes['Name'].str.contains("SA-HW10")].values



#TI=HW10[HW10['Name'].str.contains("SA")]["TI"].as_matrix()
print (TI)


# In[38]:


## changed skipfooter=0 to skipfooter=1

DFF=pd.read_excel((os.path.join(Folder, "ShearF-HW40.xlsx")),header=1,skipfooter=1)
DFA=pd.read_excel((os.path.join(Folder, "CubeA-HW40.xlsx")),header=1,skipfooter=1)
DFF

#[0] is the ND Single Row
#print 


if (((DFF["[2, 1, 1]"][0]/DFF["[2, 0, 0]"][0])>=0.869) and ((DFF["[2, 1, 1]"][0]/DFF["[2, 0, 0]"][0])<=1.195)):
    print (True, DFF["[2, 1, 1]"][0]/DFF["[2, 0, 0]"][0])
else:
    print (False, DFF["[2, 1, 1]"][0]/DFF["[2, 0, 0]"][0])

if (((DFA["[2, 2, 0]"][0]/ DFA["[2, 0, 0]"][0])>=0.800) and ((DFA["[2, 2, 0]"][0]/ DFA["[2, 0, 0]"][0])<=1.235)):
    print (True, DFA["[2, 2, 0]"][0]/ DFA["[2, 0, 0]"][0])
else:
    print (False, DFA["[2, 2, 0]"][0]/ DFA["[2, 0, 0]"][0])


# ## Search for all files, includes many HW values (this does take a bit...)
How ASTM E975 applicability is flagged:

 R (theoretical intensity) values from ASTM E975 were used
"()" indicates bounds stated in E975
10% variation -> 1.1/0.9, 0.9/1.1
 
 Ferrite peak ratio = 190.8/20.73 = 9.204
Results in bounds of 7.531 (8) < 9.204 < 11.249 (11) for ferrite 
Normailized on rounding 0.869 < 1 < 1.195

 Austenite peak ratio = 47.88/34.78 = 1.377
Results in bounds of 1.126 (1.1) < 1.377 < 1.682 (1.7) for ferrite
Normailized on rounding 0.800 < 1 < 1.235

# In[39]:


## Ensure the excel files are not open, errors will result if they are
## This step can take some time

import fnmatch

# Set Volume Fraction
VF=.25

#Read in Texture files, sorting by Austenite and Ferrite
AusteniteTextures=[]
FerriteTextures=[]  

for file in os.listdir(Folder):
    #print file
    if fnmatch.fnmatch(file, '*A-*'):
        #print file
        AusteniteTextures.append(file)
    elif fnmatch.fnmatch(file, '*F-*'):
        FerriteTextures.append(file)
    else: ()
        
#print FerriteTextures
#print AusteniteTextures

# create a dataframe shape from existing data
##test/fix, skipfooter=1, orig 0
DFA=pd.read_excel((os.path.join(Folder, "alpha1F-HW20.xlsx")),header=1,skipfooter=1)

#copy the HKL reflection indexes
df2pair=DFA["HKL"]
df4pair=DFA["HKL"]
dfMaxUnique=DFA["HKL"]

df2pair.loc[len(df2pair)]="TI"
df2pair.loc[len(df2pair)]="ASTM E975 AND"
df2pair.loc[len(df2pair)]="ASTM E975 Austenite Only"

for AustOrient in AusteniteTextures:
    for FerrOrient in FerriteTextures:
        
##may need to undo this one, changed skipfooter=0 to skipfooter=1
        DFF=pd.read_excel(os.path.join(Folder,FerrOrient),header=1,skipfooter=1)
        DFA=pd.read_excel(os.path.join(Folder,AustOrient),header=1,skipfooter=1)

#         print (DFF)
#         print (DFA)
        
        # add minus VF for plotting
        DF1=(VF*DFA["2Pairs"]/(VF*DFA["2Pairs"]+((1.0-VF)*DFF["2Pairs"])))#-VF
        DF2=(VF*DFA["4Pairs"]/(VF*DFA["4Pairs"]+((1.0-VF)*DFF["4Pairs"])))#-VF
        DF3=(VF*DFA["MaxUnique"]/(VF*DFA["MaxUnique"]+((1.0-VF)*DFF["MaxUnique"])))#-VF
        
        Fname=FerrOrient.split(".")
        Aname=AustOrient.split(".")
        
        #print Fname[0]
        
#         TIA=TextureIndexes[TextureIndexes['Name'].str.contains(Aname[0])]["TI"].as_matrix()
#         TIF=TextureIndexes[TextureIndexes['Name'].str.contains(Fname[0])]["TI"].as_matrix()
        TIA=TextureIndexes[TextureIndexes['Name'].str.contains(Aname[0])]["TI"].values
        TIF=TextureIndexes[TextureIndexes['Name'].str.contains(Fname[0])]["TI"].values
    
    
        MixName= Fname[0]+"-"+Aname[0]
        
        #print TIA, TIF
        
        AverageTI=0.5*(TIA+TIF)
        #print MixName
        #print DF1, DF2, DF3
        #print AverageTI
        
        # Add TI to datafile
        DF1.loc[len(DF1)]=AverageTI[0]
        DF2.loc[len(DF2)]=AverageTI[0]
        DF3.loc[len(DF3)]=AverageTI[0]

        # Add ASTM Flag to Datafile, only checking ND, index [0] !!!
        
        if (((DFF["[2, 1, 1]"][0]/DFF["[2, 0, 0]"][0])>=0.869) and ((DFF["[2, 1, 1]"][0]/DFF["[2, 0, 0]"][0])<=1.195)):
            #print True, DFF["[2, 1, 1]"][0]/DFF["[2, 0, 0]"][0]
            ASTMFlagF=True
        else:
            #print False, DFF["[2, 1, 1]"][0]/DFF["[2, 0, 0]"][0]
            ASTMFlagF=False

        if (((DFA["[2, 2, 0]"][0]/ DFA["[2, 0, 0]"][0])>=0.800) and ((DFA["[2, 2, 0]"][0]/ DFA["[2, 0, 0]"][0])<=1.235)):
            #print True, DFA["[2, 2, 0]"][0]/ DFA["[2, 0, 0]"][0]
            ASTMFlagA=True
        else:
            #print False, DFA["[2, 2, 0]"][0]/ DFA["[2, 0, 0]"][0]
            ASTMFlagA=False
        

        DF1.loc[len(DF1)]=(ASTMFlagA and ASTMFlagF)
        DF2.loc[len(DF2)]=(ASTMFlagA and ASTMFlagF)
        DF3.loc[len(DF3)]=(ASTMFlagA and ASTMFlagF)
        
        DF1.loc[len(DF1)]=(ASTMFlagA)
        DF2.loc[len(DF2)]=(ASTMFlagA)
        DF3.loc[len(DF3)]=(ASTMFlagA)
        
        data1 = pd.DataFrame({MixName: DF1})
        data2 = pd.DataFrame({MixName: DF2})
        data3 = pd.DataFrame({MixName: DF3})
        #print data
        
        
        df2pair = pd.concat([df2pair, data1], axis=1)
        df4pair = pd.concat([df4pair, data2], axis=1)
        dfMaxUnique = pd.concat([dfMaxUnique, data3], axis=1)
        
#print df2pair
#df2pair
#TexNames=list(df2pair.columns.values)
#TexNames[1:]
#print list(df[TexNames[0]])
#df1[TexNames[0]]
#print list(df2pair.columns.values)
#df4pair.to_excel('3Pairdata')


# ### Diagnostic checks

# In[40]:


#df2pair
dfMaxUnique
df4pair.to_excel('3Pairdata.xlsx')


# In[41]:


#list(dfMaxUnique.ix[13])[1:]
#list(dfMaxUnique.loc[13])[1:]
list(dfMaxUnique.loc[11])[1:]


# In[42]:


# diagnostic check for searching for ASTM flag
#((15*df2pair.ix[14])[1:]+1).astype(int)
((15*df2pair.loc[14])[1:]+1).astype(int)


# ### Plot all points

# In[43]:


fig, ax = plt.subplots(figsize=(10,5), dpi=600)
#black line for volume fraction
###changed .ix to .loc
#ax.plot([0, max(list(df2pair.ix[13])[1:])], [VF, VF], color='k', linestyle='-', linewidth=.5)

ax.axhline(y=VF,color='k',ls="-",zorder=5)
ax.axhline(y=VF*(1+tolerr),color='.4',ls="--",zorder=1)
ax.axhline(y=VF*(1-tolerr),color='.4',ls="--",zorder=1)
##MCchange, should change from calling "ND single" to 2nd hex scheme
#scatter1=ax.scatter(list(df2pair.loc[13])[1:] ,list(df2pair.loc[0])[1:], color='r',
#                    label="2Pair", alpha=0.33,zorder=4)

scatter1=ax.scatter(list(df2pair.loc[13])[1:] ,list(df2pair.loc[0])[1:], color='r',
                    label="2Pair", alpha=0.33,zorder=4)

ax.plot([1.06,1.06 ], [0, 1], color='k', linestyle='-', linewidth=1, label="")
#ax.scatter(list(df4pair.loc[13])[1:] ,list(df4pair.loc[11])[1:], color='y',
#           label="4Pair", alpha=0.33,zorder=4)
##MCchange, 4 to 3pair
ax.scatter(list(df4pair.loc[13])[1:] ,list(df4pair.loc[0])[1:], color='y',
           label="3Pair", alpha=0.33,zorder=4)

ax.plot([1.43725,1.43725 ], [0, 1], color='#00BFFF', linestyle='-', linewidth=1, label="")
#scatter3=ax.scatter(list(dfMaxUnique.loc[13])[1:] ,list(dfMaxUnique.loc[0])[1:] ,color='b',
scatter3=ax.scatter(list(dfMaxUnique.loc[13])[1:] ,list(dfMaxUnique.loc[0])[1:] ,color='b',
                    label="MaxUnique", alpha=0.33,zorder=4)
##end MCchange
ax.plot([2.2775,2.2775 ], [0, 1], color='m', linestyle='-', linewidth=1, label="")

ax.legend(bbox_to_anchor=(0., 1.12, 1., .102), loc=3, ncol=3,  mode="expand",
          borderaxespad=0.,fontsize=18)


ax.xaxis.set_tick_params(labelsize=16)
ax.yaxis.set_tick_params(labelsize=16)

ax.text(2.4, .9, "ASTM E975 Applicability Limit", fontsize=16)
ax.arrow(2.4, .92, -1.23, 0, head_width=0.05, head_length=0.1, fc='k', ec='k')

ax.text(2.5, .8, "TRIP 700", fontsize=16)
ax.arrow(2.5, .82, -.96, 0, head_width=0.05, head_length=0.1, fc='#00BFFF', ec='#00BFFF')

ax.text(3, .7, "TRIP 780", fontsize=16)
ax.arrow(3, .72, -.63, 0, head_width=0.05, head_length=0.1, fc='m', ec='m')

ax.set_axisbelow(True)
plt.grid()
ax.set_ylabel('Austenite Phase Fraction', fontsize=16)
ax.set_xlabel('Average Texture Index of Austenite and Ferrite', fontsize=16)
#ax.set_title('ND sampling, Effect of Average Texture Index on Phase Fraction')
ax.set_ylim([-0.05,1.05])
ax.set_xlim([0.8,10])
plt.xticks(np.arange(1, 10, 1.0))

# check at 1 the values converge.  Need to comment out the arrow text as well
#ax.set_xlim([0.999,1.01]) 
#ax.set_ylim([0.22,0.28])
#plt.xticks(np.arange(1, 2, 1.0))


# Remove comments to save files
#plt.savefig("TextureIntensity-Draft.pdf", dpi=600,format="pdf",orientation='portrait',bbox_inches='tight')
#plt.savefig("TextureIntensity-Draft.eps", dpi=600,format="eps",orientation='portrait',bbox_inches='tight')
plt.savefig("TextureIntensity-Draft.png", dpi=600,format="png",orientation='portrait',bbox_inches='tight')

plt.show()


# # Redifine Hex Grid Scheme to include subplots

# In[44]:


def HexGrid(name, chi_max, angular_spacing):
    
#chi_max=90.0  #maximum tilt angle in degrees
#angular_spacing=7.0


    d_max=2.0*math.sin(math.radians(chi_max)/2.0)

    N=d_max/math.radians(angular_spacing)

    #print "Max Tilt: ", chi_max
    #print "Angular Spacing: ", angular_spacing

    xaxis=[]
    yaxis=[]

    j=0
    i=0
    y_j=0.0
    x_ij=0.0
    chi_ij=0.0
    phi_ij=0.0

    while np.multiply((math.sqrt(3)*(d_max)/(2*N)), j ) < d_max:

        y_j=np.multiply((math.sqrt(3)*(d_max)/(2*N)), j )
        #print "Current tilt: ", math.degrees(y_j)
        #print "X_ij list: ", x_ij ,"\n"

        #print "j: ",j,"\ty_j: ", y_j 
        #print "\n"    

        #print "x_ij limit: ", (math.sqrt(d_max*d_max-(y_j*y_j)))
        i=0
        x_ij=0.0
        while ((d_max/N)*i) <= (math.sqrt(d_max*d_max-(y_j*y_j))):

            #x.append(1)
            #print "b*j: ",j,"\ti: ", i, "\tx_ij: ", x_ij 

            if (j%2==0):
                #print "Even"
                x_ij=((d_max/N)*i)
                #x[j].append((R/N)*i)
            elif (j%2==1):
                #print "Odd"
                x_ij=((d_max/(2.0*N))+(d_max/N)*i)
                #x[j].append((R/(2*N))+(R/N)*i)
            else:
                pass


            #NextRotation=((R/N)*i)   
            #print "e j: ",j,"\ti: ", i, "\tx_ij: ", x_ij 

            d_ij=math.sqrt(x_ij*x_ij + y_j*y_j)
            chi_ij=2.0*math.asin(d_ij/2.0)
            if math.degrees(chi_ij) <= chi_max:
                if y_j==0:
                    #do once to avoid duplicates
                    #print "do some nothing"

                    phi_ij=((180.0/math.pi)*math.atan2(y_j,x_ij))
                    xaxis.append(math.degrees(chi_ij))
                    yaxis.append(phi_ij)

                    if x_ij!=0:
                        phi_ij=((180.0/math.pi)*math.atan2(y_j,-1.0*x_ij))
                        xaxis.append(math.degrees(chi_ij))
                        yaxis.append(phi_ij)
                    else:
                        #removes reduntant rotation
                        pass

                else:
                    if x_ij>0:
#MCchange
                        #positive i values - quadrant I
                        phi_ij=((180.0/math.pi)*math.atan2(y_j,x_ij))
                        xaxis.append(math.degrees(chi_ij))
                        yaxis.append(phi_ij)

                        #negative i values -quadrant II
                        phi_ij=((180.0/math.pi)*math.atan2(y_j,(-1.0*x_ij)))
                        xaxis.append(math.degrees(chi_ij))
                        yaxis.append(phi_ij)

                        #quadrant III
                        phi_ij=((180.0/math.pi)*math.atan2(y_j,x_ij)+180.0)
                        xaxis.append(math.degrees(chi_ij))
                        yaxis.append(phi_ij)

                        #quadrant IV
                        phi_ij=((180.0/math.pi)*math.atan2(y_j,(-1.0*x_ij))+180.0)
                        xaxis.append(math.degrees(chi_ij))
                        yaxis.append(phi_ij)
#MCchange
                    else:
                        phi_ij=(90.0)
                        xaxis.append(math.degrees(chi_ij))
                        yaxis.append(phi_ij)

                        phi_ij=(270.0)
                        xaxis.append(math.degrees(chi_ij))
                        yaxis.append(phi_ij)
                #print "j: ",j,"\ty_j: ", y_j , "\ti: ", i,  "\tx_ij: ", x_ij ,"\td_ij: ", d_ij ,"\tchi_ij: ",math.degrees(chi_ij), "\tphi_ij: ",phi_ij
                    #anglelist.append([math.degrees(chi_ij) ,phi_ij])

            else:
                #append anyway to debug
            #    xaxis.append(math.degrees(chi_ij))
            #    yaxis.append(phi_ij)
                #print "Excceds bounds at:"
                #print "\tchi_ij: ",chi_ij, "\td_max: ",d_max
                #print "\tchi_ij: ",math.degrees(chi_ij), "\td_max: ",math.degrees(d_max)
                #print "j: ",j,"\ty_j: ", y_j , "\ti: ", i,  "\tx_ij: ", x_ij ,"\td_ij: ", d_ij ,"\tchi_ij: ",math.degrees(chi_ij), "\tphi_ij: ",phi_ij
                pass
            #NextRotation=x_ij
            #print "Next Rotation: ", NextRotation


            ### Iterate
            i=i+1
        #print "X_ij list: ", x_ij ,"\n"
        #NextTilt= np.multiply((math.sqrt(3)*(d_max)/(2*N)), j+1 )
        #print "Next Tilt:", math.degrees(NextTilt)   


        ### Iterate
        j=j+1

        #print "\n"
        #print y_j

    d = {'Tilt' : xaxis, 'Rotation' : yaxis}
    #coordinates=pd.DataFrame(d)
    #MCchange
    coordinates=pd.DataFrame(d)
    Badvals = coordinates[coordinates['Rotation'] > 270.0].index
    coordinates.drop(Badvals , inplace=True)
    Badvals2 = coordinates[coordinates['Rotation'] < 180.0].index
    coordinates.drop(Badvals2 , inplace=True)
    #MCchange

    #print coordsDF.sort_values('Tilt')    
    return name, coordinates


# ### Plot all points for Hex Grid sampling scheme
# 
# Not included in the paper, but shows the bias error is near zero for nearly all texture intensities

# In[45]:


from mpl_toolkits.axes_grid1.inset_locator import inset_axes
fig, ax = plt.subplots(figsize=(10,5), dpi=600)

#black line for volume fraction

##all .ix changed to .loc  ###MCchanged df2 to df4
ax.plot([0, max(list(df4pair.loc[13])[1:])], [VF, VF], color='k', linestyle='-', linewidth=.5)
##MCchange, include error lines
ax.axhline(y=VF*(1+tolerr),color='.4',ls="--",zorder=1)
ax.axhline(y=VF*(1-tolerr),color='.4',ls="--",zorder=1)
# ax.axhline(x=2.597,color='.7',ls="--",zorder=1)
# ax.axhline(x=1.457,color='.7',ls="--",zorder=1)
# ax.axhline(x=1.564,color='.7',ls="--",zorder=1)
# ax.axhline(x=1.530,color='.7',ls="--",zorder=1)

#0 is ND, 10 - 5° fullfine, 11 - 22.5° fullcourse, 12- 5° partial
##MCchange, removed 2pair from plot, changed 4pair to 3pair
#scatter1=ax.scatter(list(df2pair.loc[13])[1:] ,list(df2pair.loc[10])[1:], color='r',label="2Pair", alpha=0.33)

ax.plot([2.28,2.28 ], [0, 1], color='#3487fa', linestyle='-', linewidth=1, label="TRIP 780")
ax.plot([1.44,1.44 ], [0, 1], color='#1c03ff', linestyle='-', linewidth=1, label="TRIP 700")
ax.plot([1.42,1.42 ], [0, 1], color='#ff03ff', linestyle='-', linewidth=1, label="Duplex")
ax.plot([1.48,1.48 ], [0, 1], color='#00BFFF', linestyle='-', linewidth=1, label="TRIP 780b")

ax.scatter(list(df4pair.loc[14])[1:] ,list(df4pair.loc[10])[1:], color='g',label="3Pair: 200, 211, and 310 Ferrite; 200, 220, 222 Austenite", alpha=0.33)



#scatter3=ax.scatter(list(dfMaxUnique.loc[13])[1:] ,list(dfMaxUnique.loc[10])[1:] ,color='b',label="MaxUnique", alpha=0.33)



ax.legend(bbox_to_anchor=(0., 1.12, 1., .102), loc=3, ncol=3,  mode="expand", borderaxespad=0.)


##MCchange ##changed second index in labels to fit within new plot limits

#ax.text(2, .38, "ASTM E975 Applicability Limit", fontsize=12)
#ax.arrow(2, .92, -0.83, 0, head_width=0.05, head_length=0.1, fc='k', ec='k')

#ax.text(2.5, .35, "TRIP 700", fontsize=12)
#ax.arrow(2.5, .82, -.96, 0, head_width=0.05, head_length=0.1, fc='#00BFFF', ec='#00BFFF')

#ax.text(3, .32, "TRIP 780", fontsize=12)
#ax.arrow(3, .72, -.63, 0, head_width=0.05, head_length=0.1, fc='m', ec='m')

### For rescale figures
#ax.text(2, 0.2615, "ASTM E975 Applicability Limit", fontsize=12)
#ax.arrow(2, 0.2615, -0.83, 0, head_width=0.0005, head_length=0.01, fc='k', ec='k')

#ax.text(2.5, 0.2595, "TRIP 700", fontsize=12)
#ax.arrow(2.5, 0.2595, -.96, 0, head_width=0.0005, head_length=0.01, fc='#00BFFF', ec='#00BFFF', lw=1)

#ax.text(3, 0.2575, "TRIP 780", fontsize=12)
#ax.arrow(3, 0.2575, -.63, 0, head_width=0.0005, head_length=0.01, fc='m', ec='m')

ax.set_axisbelow(True)
plt.grid()
ax.set_ylabel('Austenite Phase Fraction')
ax.set_xlabel('Average Texture Index of Austenite and Ferrite')
#ax.set_title('ND sampling, Effect of Average Texture Index on Phase Fraction')
ax.set_ylim([0.15,0.4])
###MCchange  ##ax.set_ylim([-0.05,1.05])
#ax.set_ylim([0.2375,0.2625])
ax.set_xlim([0.8,5])
plt.xticks(np.arange(1, 5, 1.0))

# check at 1 the values converge.  Need to comment out the arrow text as well
#ax.set_xlim([0.999,1.01]) 
#ax.set_ylim([0.22,0.28])
#plt.xticks(np.arange(1, 2, 1.0))


##MCchange, attempt to add subplot
ax1 = fig.add_subplot(339,position=[0.5, 0.65, 0.2, 0.2],projection='stereonet')
ax1.annotate('(i)', xy=(0, 0), xytext=(-3.0,0.1), fontsize=16)
ax1.plane(0.0, 90.0, 'k-', linewidth=1)
ax1.plane(90.0, 90.0, 'k-', linewidth=1)
ax1.set_azimuth_ticks([90,0], labels=['',''])
#SchemeName,Coordinates=HexGrid("HexGrid-90degTilt15degRes",90.0,15)
SchemeName,Coordinates=HexGrid("HexGrid-90degTilt5degRes",90.0,5.0)

for index, row in Coordinates.iterrows():
    #print row['Tilt'], row['Rotation']
    dip, strike =  row['Tilt'], (row['Rotation']-90.0)
    ax1.pole(strike, dip, 'mh', markersize=2, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)
#ax.Position[0.3,0.3]
#ax1.subplots_adjust(left=None, bottom=0.32, right=None, top=0.39, wspace=None, hspace=None)    
#MCchanges

#change name to fit the choice of grid
#plt.savefig("TextureIntensity-FullFineHex.pdf", dpi=600,format="pdf",orientation='portrait',bbox_inches='tight')
#plt.savefig("TextureIntensity-FullFineHex.eps", dpi=600,format="eps",orientation='portrait',bbox_inches='tight')
plt.savefig("TextureIntensity-FullFineHex5degsFull.png", dpi=600,format="png",orientation='portrait',bbox_inches='tight')
plt.show()


####MCchange; attempt to plot all hex schemes
fig, ax = plt.subplots(figsize=(10,5), dpi=600)
#black line for volume fraction

##all .ix changed to .loc
ax.plot([0, max(list(df4pair.loc[14])[1:])], [VF, VF], color='k', linestyle='-', linewidth=.5)
##MCchange, include error lines
ax.axhline(y=VF*(1+tolerr),color='.4',ls="--",zorder=1)
ax.axhline(y=VF*(1-tolerr),color='.4',ls="--",zorder=1)
# ax.axhline(x=2.597,color='.7',ls="--",zorder=1)
# ax.axhline(x=1.457,color='.7',ls="--",zorder=1)
# ax.axhline(x=1.564,color='.7',ls="--",zorder=1)
# ax.axhline(x=1.530,color='.7',ls="--",zorder=1)
ax.plot([2.28,2.28 ], [0, 1], color='#3487fa', linestyle='-', linewidth=1, label="TRIP 780")
ax.plot([1.44,1.44 ], [0, 1], color='#1c03ff', linestyle='-', linewidth=1, label="TRIP 700")
ax.plot([1.42,1.42 ], [0, 1], color='#ff03ff', linestyle='-', linewidth=1, label="Duplex")
ax.plot([1.48,1.48 ], [0, 1], color='#00BFFF', linestyle='-', linewidth=1, label="TRIP 780b")
#0 is ND, 10 - 5° fullfine, 11 - 22.5° fullcourse, 12- 5° partial
##MCchange, removed 2Pair from plot
#scatter1=ax.scatter(list(df2pair.loc[13])[1:] ,list(df2pair.loc[11])[1:], color='r',label="2Pair", alpha=0.33)

##MCchange, ax.plot below produces the vertical marker bars for diff. materials and standars

#ax.plot([1.06,1.06 ], [0, 1], color='k', linestyle='-', linewidth=1, label="")

ax.scatter(list(df4pair.loc[14])[1:] ,list(df4pair.loc[11])[1:], color='g',label="3Pair: 200, 211, and 310 Ferrite; 200, 220, 222 Austenite", alpha=0.33)

#ax.plot([1.43725,1.43725 ], [0, 1], color='#00BFFF', linestyle='-', linewidth=1, label="")

#scatter3=ax.scatter(list(dfMaxUnique.loc[13])[1:] ,list(dfMaxUnique.loc[11])[1:] ,color='b',label="MaxUnique", alpha=0.33)

#ax.plot([2.2775,2.2775 ], [0, 1], color='m', linestyle='-', linewidth=1, label="")

ax.legend(bbox_to_anchor=(0., 1.12, 1., .102), loc=3, ncol=3,  mode="expand", borderaxespad=0.)


##MCchange ##changed second index in labels to fit within new plot limits

#ax.text(2, .38, "ASTM E975 Applicability Limit", fontsize=12)
#ax.arrow(2, .92, -0.83, 0, head_width=0.05, head_length=0.1, fc='k', ec='k')

#ax.text(2.5, .35, "TRIP 700", fontsize=12)
#ax.arrow(2.5, .82, -.96, 0, head_width=0.05, head_length=0.1, fc='#00BFFF', ec='#00BFFF')

#ax.text(3, .32, "TRIP 780", fontsize=12)
#ax.arrow(3, .72, -.63, 0, head_width=0.05, head_length=0.1, fc='m', ec='m')

### For rescale figures
#ax.text(2, 0.2615, "ASTM E975 Applicability Limit", fontsize=12)
#ax.arrow(2, 0.2615, -0.83, 0, head_width=0.0005, head_length=0.01, fc='k', ec='k')

#ax.text(2.5, 0.2595, "TRIP 700", fontsize=12)
#ax.arrow(2.5, 0.2595, -.96, 0, head_width=0.0005, head_length=0.01, fc='#00BFFF', ec='#00BFFF', lw=1)

#ax.text(3, 0.2575, "TRIP 780", fontsize=12)
#ax.arrow(3, 0.2575, -.63, 0, head_width=0.0005, head_length=0.01, fc='m', ec='m')

ax.set_axisbelow(True)
plt.grid()
ax.set_ylabel('Austenite Phase Fraction')
ax.set_xlabel('Average Texture Index of Austenite and Ferrite')
#ax.set_title('ND sampling, Effect of Average Texture Index on Phase Fraction')
ax.set_ylim([0.15,0.4])
#ax.set_ylim([0.2375,0.2625])
ax.set_xlim([0.8,5])
plt.xticks(np.arange(1, 5, 1.0))

# check at 1 the values converge.  Need to comment out the arrow text as well
#ax.set_xlim([0.999,1.01]) 
#ax.set_ylim([0.22,0.28])
#plt.xticks(np.arange(1, 2, 1.0))

##MCchange, attempt to add subplot
ax2 = fig.add_subplot(339, position=[0.5, 0.65, 0.2, 0.2], projection='stereonet')
ax2.annotate('(i)', xy=(0, 0), xytext=(-3.0,0.1), fontsize=16)
ax2.plane(0.0, 90.0, 'k-', linewidth=1)
ax2.plane(90.0, 90.0, 'k-', linewidth=1)
ax2.set_azimuth_ticks([90,0], labels=['',''])
SchemeName,Coordinates=HexGrid("HexGrid-90degTilt15degRes",90.0,15)
#SchemeName,Coordinates=HexGrid("HexGrid-60degTilt5degRes",60.0,5.0)

for index, row in Coordinates.iterrows():
    #print row['Tilt'], row['Rotation']
    dip, strike =  row['Tilt'], (row['Rotation']-90.0)
    ax2.pole(strike, dip, 'mh', markersize=2, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)
#MCchanges


#change name to fit the choice of grid
#plt.savefig("TextureIntensity-FullFineHex.pdf", dpi=600,format="pdf",orientation='portrait',bbox_inches='tight')
#plt.savefig("TextureIntensity-FullFineHex.eps", dpi=600,format="eps",orientation='portrait',bbox_inches='tight')
plt.savefig("TextureIntensity-FullFineHex15degsFull.png", dpi=600,format="png",orientation='portrait',bbox_inches='tight')
plt.show()

####MCchange; attempt to plot all hex schemes
fig, ax = plt.subplots(figsize=(10,5), dpi=600)
#black line for volume fraction

##all .ix changed to .loc
ax.plot([0, max(list(df4pair.loc[14])[1:])], [VF, VF], color='k', linestyle='-', linewidth=.5)
##MCchange, include error lines
ax.axhline(y=VF*(1+tolerr),color='.4',ls="--",zorder=1)
ax.axhline(y=VF*(1-tolerr),color='.4',ls="--",zorder=1)
# ax.axhline(x=2.597,color='.7',ls="--",zorder=1)
# ax.axhline(x=1.457,color='.7',ls="--",zorder=1)
# ax.axhline(x=1.564,color='.7',ls="--",zorder=1)
# ax.axhline(x=1.530,color='.7',ls="--",zorder=1)
ax.plot([2.28,2.28 ], [0, 1], color='#3487fa', linestyle='-', linewidth=1, label="TRIP 780")
ax.plot([1.44,1.44 ], [0, 1], color='#1c03ff', linestyle='-', linewidth=1, label="TRIP 700")
ax.plot([1.42,1.42 ], [0, 1], color='#ff03ff', linestyle='-', linewidth=1, label="Duplex")
ax.plot([1.48,1.48 ], [0, 1], color='#00BFFF', linestyle='-', linewidth=1, label="TRIP 780b")

#0 is ND, 10 - 5° fullfine, 11 - 22.5° fullcourse, 12- 5° partial
##MCchange; removed 2pair from plot
#scatter1=ax.scatter(list(df2pair.loc[13])[1:] ,list(df2pair.loc[12])[1:], color='r',label="2Pair", alpha=0.33)

#ax.plot([1.06,1.06 ], [0, 1], color='k', linestyle='-', linewidth=1, label="")

ax.scatter(list(df4pair.loc[14])[1:] ,list(df4pair.loc[12])[1:], color='g',label="3Pair: 200, 211, and 310 Ferrite; 200, 220, 222 Austenite", alpha=0.33)

#ax.plot([1.43725,1.43725 ], [0, 1], color='#00BFFF', linestyle='-', linewidth=1, label="")

#scatter3=ax.scatter(list(dfMaxUnique.loc[13])[1:] ,list(dfMaxUnique.loc[12])[1:] ,color='b',label="MaxUnique", alpha=0.33)

#ax.plot([2.2775,2.2775 ], [0, 1], color='m', linestyle='-', linewidth=1, label="")

ax.legend(bbox_to_anchor=(0., 1.12, 1., .102), loc=3, ncol=3,  mode="expand", borderaxespad=0.)


##MCchange ##changed second index in labels to fit within new plot limits, commented out arrow and text

#ax.text(2, .38, "ASTM E975 Applicability Limit", fontsize=12)
#ax.arrow(2, .92, -0.83, 0, head_width=0.05, head_length=0.1, fc='k', ec='k')

#ax.text(2.5, .35, "TRIP 700", fontsize=12)
#ax.arrow(2.5, .82, -.96, 0, head_width=0.05, head_length=0.1, fc='#00BFFF', ec='#00BFFF')

#ax.text(3, .32, "TRIP 780", fontsize=12)
#ax.arrow(3, .72, -.63, 0, head_width=0.05, head_length=0.1, fc='m', ec='m')

### For rescale figures
#ax.text(2, 0.2615, "ASTM E975 Applicability Limit", fontsize=12)
#ax.arrow(2, 0.2615, -0.83, 0, head_width=0.0005, head_length=0.01, fc='k', ec='k')

#ax.text(2.5, 0.2595, "TRIP 700", fontsize=12)
#ax.arrow(2.5, 0.2595, -.96, 0, head_width=0.0005, head_length=0.01, fc='#00BFFF', ec='#00BFFF', lw=1)

#ax.text(3, 0.2575, "TRIP 780", fontsize=12)
#ax.arrow(3, 0.2575, -.63, 0, head_width=0.0005, head_length=0.01, fc='m', ec='m')

ax.set_axisbelow(True)
plt.grid()
ax.set_ylabel('Austenite Phase Fraction')
ax.set_xlabel('Average Texture Index of Austenite and Ferrite')
#ax.set_title('ND sampling, Effect of Average Texture Index on Phase Fraction')
ax.set_ylim([0.15,0.4])
#ax.set_ylim([0.2375,0.2625])
ax.set_xlim([0.8,5])
plt.xticks(np.arange(1, 5, 1.0))

# check at 1 the values converge.  Need to comment out the arrow text as well
#ax.set_xlim([0.999,1.01]) 
#ax.set_ylim([0.22,0.28])
#plt.xticks(np.arange(1, 2, 1.0))

##MCchange, attempt to add subplot
ax3 = fig.add_subplot(339, position=[0.5, 0.65, 0.2, 0.2], projection='stereonet')
ax3.annotate('(i)', xy=(0, 0), xytext=(-3.0,0.1), fontsize=16)
ax3.plane(0.0, 90.0, 'k-', linewidth=1)
ax3.plane(90.0, 90.0, 'k-', linewidth=1)
ax3.set_azimuth_ticks([90,0], labels=['',''])
SchemeName,Coordinates=HexGrid("HexGrid-90degTilt22p5degRes",90.0,22.5)
#SchemeName,Coordinates=HexGrid("HexGrid-60degTilt5degRes",60.0,5.0)

for index, row in Coordinates.iterrows():
    #print row['Tilt'], row['Rotation']
    dip, strike =  row['Tilt'], (row['Rotation']-90.0)
    ax3.pole(strike, dip, 'mh', markersize=2, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)
#MCchanges

#change name to fit the choice of grid
#plt.savefig("TextureIntensity-FullFineHex.pdf", dpi=600,format="pdf",orientation='portrait',bbox_inches='tight')
#plt.savefig("TextureIntensity-FullFineHex.eps", dpi=600,format="eps",orientation='portrait',bbox_inches='tight')
plt.savefig("TextureIntensity-FullFineHex22p5degsFull.png", dpi=600,format="png",orientation='portrait',bbox_inches='tight')
plt.show()




fig, ax = plt.subplots(figsize=(10,5), dpi=600)
#black line for volume fraction

##all .ix changed to .loc
ax.plot([0, max(list(df2pair.loc[14])[1:])], [VF, VF], color='k', linestyle='-', linewidth=.5)
##MCchange, include error lines
ax.axhline(y=VF*(1+tolerr),color='.4',ls="--",zorder=1)
ax.axhline(y=VF*(1-tolerr),color='.4',ls="--",zorder=1)

ax.plot([2.28,2.28 ], [0, 1], color='#3487fa', linestyle='-', linewidth=1, label="TRIP 780")
ax.plot([1.44,1.44 ], [0, 1], color='#1c03ff', linestyle='-', linewidth=1, label="TRIP 700")
ax.plot([1.42,1.42 ], [0, 1], color='#ff03ff', linestyle='-', linewidth=1, label="Duplex")
ax.plot([1.48,1.48 ], [0, 1], color='#00BFFF', linestyle='-', linewidth=1, label="TRIP 780b")

#0 is ND, 10 - 5° fullfine, 11 - 22.5° fullcourse, 12- 5° partial
##MCchange; removed 2pair from plot
#scatter1=ax.scatter(list(df2pair.loc[13])[1:] ,list(df2pair.loc[12])[1:], color='r',label="2Pair", alpha=0.33)

#ax.plot([1.06,1.06 ], [0, 1], color='k', linestyle='-', linewidth=1, label="")

ax.scatter(list(df4pair.loc[14])[1:] ,list(df4pair.loc[13])[1:], color='g',label="3Pair: 200, 211, and 310 Ferrite; 200, 220, 222 Austenite", alpha=0.33)

#ax.plot([1.43725,1.43725 ], [0, 1], color='#00BFFF', linestyle='-', linewidth=1, label="")

#scatter3=ax.scatter(list(dfMaxUnique.loc[13])[1:] ,list(dfMaxUnique.loc[12])[1:] ,color='b',label="MaxUnique", alpha=0.33)

#ax.plot([2.2775,2.2775 ], [0, 1], color='m', linestyle='-', linewidth=1, label="")

ax.legend(bbox_to_anchor=(0., 1.12, 1., .102), loc=3, ncol=3,  mode="expand", borderaxespad=0.)


##MCchange ##changed second index in labels to fit within new plot limits, commented out arrow and text

#ax.text(2, .38, "ASTM E975 Applicability Limit", fontsize=12)
#ax.arrow(2, .92, -0.83, 0, head_width=0.05, head_length=0.1, fc='k', ec='k')

#ax.text(2.5, .35, "TRIP 700", fontsize=12)
#ax.arrow(2.5, .82, -.96, 0, head_width=0.05, head_length=0.1, fc='#00BFFF', ec='#00BFFF')

#ax.text(3, .32, "TRIP 780", fontsize=12)
#ax.arrow(3, .72, -.63, 0, head_width=0.05, head_length=0.1, fc='m', ec='m')

### For rescale figures
#ax.text(2, 0.2615, "ASTM E975 Applicability Limit", fontsize=12)
#ax.arrow(2, 0.2615, -0.83, 0, head_width=0.0005, head_length=0.01, fc='k', ec='k')

#ax.text(2.5, 0.2595, "TRIP 700", fontsize=12)
#ax.arrow(2.5, 0.2595, -.96, 0, head_width=0.0005, head_length=0.01, fc='#00BFFF', ec='#00BFFF', lw=1)

#ax.text(3, 0.2575, "TRIP 780", fontsize=12)
#ax.arrow(3, 0.2575, -.63, 0, head_width=0.0005, head_length=0.01, fc='m', ec='m')

ax.set_axisbelow(True)
plt.grid()
ax.set_ylabel('Austenite Phase Fraction')
ax.set_xlabel('Average Texture Index of Austenite and Ferrite')
#ax.set_title('ND sampling, Effect of Average Texture Index on Phase Fraction')
ax.set_ylim([0.15,0.4])
#ax.set_ylim([0.2375,0.2625])
ax.set_xlim([0.8,5])
plt.xticks(np.arange(1, 5, 1.0))

# check at 1 the values converge.  Need to comment out the arrow text as well
#ax.set_xlim([0.999,1.01]) 
#ax.set_ylim([0.22,0.28])
#plt.xticks(np.arange(1, 2, 1.0))

##MCchange, attempt to add subplot
ax4 = fig.add_subplot(339, position=[0.5, 0.65, 0.2, 0.2], projection='stereonet')
ax4.annotate('(i)', xy=(0, 0), xytext=(-3.0,0.1), fontsize=16)
ax4.plane(0.0, 90.0, 'k-', linewidth=1)
ax4.plane(90.0, 90.0, 'k-', linewidth=1)
ax4.set_azimuth_ticks([90,0], labels=['',''])
SchemeName,Coordinates=HexGrid("HexGrid-90degTilt30degRes",90.0,30)
#SchemeName,Coordinates=HexGrid("HexGrid-60degTilt5degRes",60.0,5.0)

for index, row in Coordinates.iterrows():
    #print row['Tilt'], row['Rotation']
    dip, strike =  row['Tilt'], (row['Rotation']-90.0)
    ax4.pole(strike, dip, 'mh', markersize=2, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)
#MCchanges

#change name to fit the choice of grid
#plt.savefig("TextureIntensity-FullFineHex.pdf", dpi=600,format="pdf",orientation='portrait',bbox_inches='tight')
#plt.savefig("TextureIntensity-FullFineHex.eps", dpi=600,format="eps",orientation='portrait',bbox_inches='tight')
plt.savefig("TextureIntensity-FullFineHex30degsFull.png", dpi=600,format="png",orientation='portrait',bbox_inches='tight')
plt.show()


# ##Add the final plot of maximum TI's for each grid
# fig, ax = plt.subplots(figsize=(10,5), dpi=600)
# #black line for volume fraction

# ##all .ix changed to .loc
# ax.plot([0, max(list(df2pair.loc[14])[1:])], [VF, VF], color='k', linestyle='-', linewidth=.5)
# ##MCchange, include error lines
# ax.axhline(y=VF*(1+tolerr),color='.4',ls="--",zorder=1)
# ax.axhline(y=VF*(1-tolerr),color='.4',ls="--",zorder=1)

# ax.plot([2.597,2.597 ], [0, 1], color='k', linestyle='-', linewidth=1, label="TRIP 780")
# ax.plot([1.457,1.457 ], [0, 1], color='#b35e19', linestyle='-', linewidth=1, label="TRIP 700")
# ax.plot([1.564,1.564 ], [0, 1], color='m', linestyle='-', linewidth=1, label="Duplex")
# ax.plot([1.530,1.530 ], [0, 1], color='#00BFFF', linestyle='-', linewidth=1, label="TRIP 780b")

# ##MCchange, include maxTI markers
# ax.plot([1.28465,1.28465 ], [0, 1], color='#f0300e', linestyle='--', linewidth=3, label="5$^\circ$ Partial")
# ax.plot([1.28465,1.28465 ], [0, 1], color='#9be8d4', linestyle='--', linewidth=6, label="15$^\circ$ Partial")
# ax.plot([1.34965,1.34965 ], [0, 1], color='#22787d', linestyle='--', linewidth=3, label="22.5$^\circ$ Partial")
# ax.plot([1.0945,1.0945 ], [0, 1], color='#00BFFF', linestyle='--', linewidth=3, label="30$^\circ$ Full")
# ax.plot([3.27655,3.27655 ], [0, 1], color='#51e8ac', linestyle='--', linewidth=3, label="15$^\circ$ Full")
# ax.plot([1.3504,1.3504 ], [0, 1], color='#1b7df5', linestyle='--', linewidth=3, label="22.5$^\circ$ Full")
# ax.plot([2.01905,2.01905 ], [0, 1], color='#c710e3', linestyle='--', linewidth=3, label="30$^\circ$ Full")


# #0 is ND, 10 - 5° fullfine, 11 - 22.5° fullcourse, 12- 5° partial
# ##MCchange; removed 2pair from plot
# #scatter1=ax.scatter(list(df2pair.loc[13])[1:] ,list(df2pair.loc[12])[1:], color='r',label="2Pair", alpha=0.33)


# ax.legend(bbox_to_anchor=(0., 1.12, 1., .102), loc=3, ncol=3,  mode="expand", borderaxespad=0.)


# ##MCchange ##changed second index in labels to fit within new plot limits, commented out arrow and text

# #ax.text(2.5, .35, "TRIP 700", fontsize=12)
# #ax.arrow(2.5, .82, -.96, 0, head_width=0.05, head_length=0.1, fc='#00BFFF', ec='#00BFFF')

# #ax.text(3, .32, "TRIP 780", fontsize=12)
# #ax.arrow(3, .72, -.63, 0, head_width=0.05, head_length=0.1, fc='m', ec='m')

# ### For rescale figures
# #ax.text(2, 0.2615, "ASTM E975 Applicability Limit", fontsize=12)
# # ax.arrow(1.28465, 0.33, -0.83, 0, head_width=0.0005, head_length=0.01, fc='#9be8d4', ec='#f0300e')
# # ax.arrow(1.34965, 0.30, -0.83, 0, head_width=0.0005, head_length=0.01, fc='#22787d', ec='#22787d')
# # ax.arrow(1.0945, 0.27, -0.83, 0, head_width=0.0005, head_length=0.01, fc='#00BFFF', ec='#00BFFF')
# # ax.arrow(3.27655, 0.24, -0.83, 0, head_width=0.0005, head_length=0.01, fc='#51e8ac', ec='#51e8ac')
# # ax.arrow(1.3504, 0.21, -0.83, 0, head_width=0.0005, head_length=0.01, fc='#1b7df5', ec='#1b7df5')
# # ax.arrow(2.01905, 0.18, -0.83, 0, head_width=0.0005, head_length=0.01, fc='#c710e3', ec='#c710e3')



# #ax.text(2.5, 0.2595, "TRIP 700", fontsize=12)
# #ax.arrow(2.5, 0.2595, -.96, 0, head_width=0.0005, head_length=0.01, fc='#00BFFF', ec='#00BFFF', lw=1)

# #ax.text(3, 0.2575, "TRIP 780", fontsize=12)
# #ax.arrow(3, 0.2575, -.63, 0, head_width=0.0005, head_length=0.01, fc='m', ec='m')

# ax.set_axisbelow(True)
# plt.grid()
# ax.set_ylabel('Austenite Phase Fraction')
# ax.set_xlabel('Average Texture Index of Austenite and Ferrite')
# #ax.set_title('ND sampling, Effect of Average Texture Index on Phase Fraction')
# ax.set_ylim([0.15,0.4])
# #ax.set_ylim([0.2375,0.2625])
# ax.set_xlim([0.8,5])
# plt.xticks(np.arange(1, 5, 1.0))


# plt.savefig("MaxTI.png", dpi=600,format="png",orientation='portrait',bbox_inches='tight')
# plt.show()


# ## Plot which ones would qualify as ASTM E975 (austenite and ferrite)

# ### plot using if statement

# In[19]:


#df2pair.ix[14][1]  ##changed to .loc
df2pair.loc[14][1]


# In[20]:


#searches through datafile, takes some time...


# In[21]:


##changed .ix to .loc

df2pairTI=[]
df2pairVF=[]

df4pairTI=[]
df4pairVF=[]

dfMaxUniqueTI=[]
dfMaxUniqueVF=[]

for counter, value in enumerate(list(df2pair.loc[14])):
    if df2pair.loc[14][counter]==1.0:
        #print df2pair.ix[14][counter], df2pair.ix[13][counter] ,df2pair.ix[0][counter]
        df2pairTI.append((df2pair.loc[13][counter]).tolist())
        df2pairVF.append((df2pair.loc[0][counter]).tolist())
        df4pairTI.append((df4pair.loc[13][counter]).tolist())
        df4pairVF.append((df4pair.loc[0][counter]).tolist())        
        dfMaxUniqueTI.append((dfMaxUnique.loc[13][counter]).tolist())
        dfMaxUniqueVF.append((dfMaxUnique.loc[0][counter]).tolist())       
        
        
        #a=list(df2pair.ix[13][counter])
        #print "Test:", a
        #df2pairTI.extend(a)
    else:
        ()
#print df2pairTI


# ## Plots 10a and 10b
# 
# 

# In[22]:




fig = plt.figure(figsize=(12,5), dpi=600)

ax1 = fig.add_subplot(122) #note the flipped order




#black line for volume fraction
#ax1.plot([0, max(list(df2pair.ix[13])[1:])], [VF, VF], color='k', linestyle='-', linewidth=.5)
ax1.axhline(y=VF,color='k',ls="-",zorder=5)
ax1.axhline(y=VF*(1+tolerr),color='.4',ls="--",zorder=1)
ax1.axhline(y=VF*(1-tolerr),color='.4',ls="--",zorder=1)


#ax.plot([1.43725,1.43725 ], [0, 1], color='m', linestyle='-', linewidth=1, label="TRIP700")
#ax.plot([2.2775,2.2775 ], [0, 1], color='c', linestyle='-', linewidth=1, label="TRIP780")

# using the size of the points as a screen for if the data meets ASTM E975
##MCchange, changed 4pair to 3pair


scatter1=ax1.scatter(df2pairTI ,df2pairVF, s= 20,   color='r',label="2Pair", alpha=0.33,zorder=4)
scatter2=ax1.scatter(df4pairTI ,df4pairVF, s= 20, color='y',label="3Pair", alpha=0.33,zorder=4)
scatter3=ax1.scatter(dfMaxUniqueTI ,dfMaxUniqueVF , s= 20, color='b',label="MaxUnique", alpha=0.33,zorder=4)

ax1.xaxis.set_tick_params(labelsize=16)
ax1.yaxis.set_tick_params(labelsize=16)

plt.grid()
#ax1.annotate('(b)', xy=(1, 0.275), xytext=(1,0.275), fontsize=16)
ax1.annotate('(b)', xy=(1, 0.275), xycoords='data',xytext=(-.15,1.1),
             textcoords='axes fraction',fontsize=20)
#ax1.set_ylabel('Austenite Phase Fraction')
ax1.set_xlabel('Average Texture Index of \nAustenite and Ferrite', fontsize=16)
#ax.set_title('ND sampling, Austenite and Ferrite Limits (ASTM E975)')
#ax.set_ylim([-0.05,1.05])
#ax.set_xlim([0,10])

#plt.xticks(np.arange(0, 10, 1.0))

ax1.set_ylim([0.22,.28])
ax1.set_xlim([.998,1.06])

ax2 = fig.add_subplot(121)


ax2.xaxis.set_tick_params(labelsize=16)
ax2.yaxis.set_tick_params(labelsize=16)

ax2.axhline(y=VF,color='k',ls="-",zorder=5)
ax2.axhline(y=VF*(1+tolerr),color='.4',ls="--",zorder=1)
ax2.axhline(y=VF*(1-tolerr),color='.4',ls="--",zorder=1)

scatter1=ax2.scatter(df2pairTI ,df2pairVF, s= 20,   color='r',label="2Pair", alpha=0.33,zorder=4)
scatter2=ax2.scatter(df4pairTI ,df4pairVF, s= 20, color='y',label="3Pair", alpha=0.33,zorder=4)
scatter3=ax2.scatter(dfMaxUniqueTI ,dfMaxUniqueVF , s= 20, color='b',label="MaxUnique", alpha=0.33,zorder=4)


#ax2.annotate('(a)', xy=(.999796, 0.275), xytext=(1,0.275), fontsize=16)
ax2.annotate('(a)', xy=(1, 0.275), xycoords='data',xytext=(-.15,1.1),
             textcoords='axes fraction',fontsize=20)


ax2.set_ylabel('Austenite Phase Fraction', fontsize=16)
ax2.set_xlabel('Average Texture Index of \nAustenite and Ferrite', fontsize=16)

ax2.xaxis.set_major_formatter(FormatStrFormatter('%.4f'))

ax2.set_ylim([0.22,.28])
ax2.set_xlim([.9998,1.002])
plt.grid()

ax2.legend(bbox_to_anchor=(.22, 1.2, 1., .102), loc=3, ncol=3, 
           borderaxespad=0., fontsize=18)

fig.subplots_adjust(wspace=.3)

# Remove comments to save files
#plt.savefig("MeetsASTME975Zoom-Draft.pdf", dpi=300,format="pdf",orientation='portrait',bbox_inches='tight')
#plt.savefig("MeetsASTME975Zoom-Draft.eps", dpi=600,format="eps",orientation='portrait',bbox_inches='tight')
plt.savefig("MeetsASTME975Zoom-Draft.png", dpi=600,format="png",orientation='portrait',bbox_inches='tight')

plt.show()


# # Unused plots

# ### Original plot, just minimized the size of the points, but still plotted them

# In[23]:


fig, ax = plt.subplots(figsize=(10,5), dpi=600)
#black line for volume fraction
##adjusted .ix to .loc due to deprecation
ax.plot([0, max(list(df2pair.loc[13])[1:])], [VF, VF], color='k', linestyle='-', linewidth=.5)

#ax.plot([1.43725,1.43725 ], [0, 1], color='m', linestyle='-', linewidth=1, label="TRIP700")
#ax.plot([2.2775,2.2775 ], [0, 1], color='c', linestyle='-', linewidth=1, label="TRIP780")

# using the size of the points as a screen for if the data meets ASTM E975
##MCchange, changed 4pair to 3pair
scatter1=ax.scatter(list(df2pair.loc[13])[1:] ,list(df2pair.loc[0])[1:], s= ((20*df2pair.loc[14])[1:]).astype(int),
                            color='r',label="2Pair", alpha=0.33)
scatter2=ax.scatter(list(df4pair.loc[13])[1:] ,list(df4pair.loc[0])[1:],  s= ((20*df2pair.loc[14])[1:]).astype(int),
                            color='y',label="3Pair", alpha=0.33)
scatter3=ax.scatter(list(dfMaxUnique.loc[13])[1:] ,list(dfMaxUnique.loc[0])[1:] , s= ((20*df2pair.loc[14])[1:]).astype(int),
                            color='b',label="MaxUnique", alpha=0.33)

ax.legend(bbox_to_anchor=(0., 1.12, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.)

plt.grid()
ax.set_ylabel('Phase Fraction')
ax.set_xlabel('Average Texture Index of Austenite and Ferrite')
#ax.set_title('ND sampling, Austenite and Ferrite Limits (ASTM E975)')
#ax.set_ylim([-0.05,1.05])
#ax.set_xlim([0,10])

#plt.xticks(np.arange(0, 10, 1.0))

ax.set_ylim([0.22,.28])
ax.set_xlim([.99,1.07])

plt.savefig("MeetsASTME975-Draft.pdf", dpi=300,format="pdf",orientation='portrait',bbox_inches='tight')
plt.savefig("MeetsASTME975-Draft.eps", dpi=600,format="eps",orientation='portrait',bbox_inches='tight')
plt.show()


# ## Austenite Only Limit

# In[24]:


fig, ax = plt.subplots(figsize=(10,5), dpi=600)
#black line for volume fraction
##replaced .ix with .loc to match module update
ax.plot([0, max(list(df2pair.loc[13])[1:])], [VF, VF], color='k', linestyle='-', linewidth=.5)

ax.plot([1.43725,1.43725 ], [0, 1], color='m', linestyle='-', linewidth=1, label="TRIP700")
ax.plot([2.2775,2.2775 ], [0, 1], color='c', linestyle='-', linewidth=1, label="TRIP780")

##MCchange, changed 4pair to 3pair
scatter1=ax.scatter(list(df2pair.loc[13])[1:] ,list(df2pair.loc[0])[1:], s= ((10*df2pair.loc[15])[1:]).astype(int),
                            color='r',label="2Pair")
scatter2=ax.scatter(list(df4pair.loc[13])[1:] ,list(df4pair.loc[0])[1:],  s= ((10*df4pair.loc[15])[1:]).astype(int),
                            color='y',label="3Pair")
scatter3=ax.scatter(list(dfMaxUnique.loc[13])[1:] ,list(dfMaxUnique.loc[0])[1:] , s= ((10*dfMaxUnique.loc[15])[1:]).astype(int),
                            color='b',label="MaxUnique")

ax.legend(bbox_to_anchor=(0., 1.12, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.)

plt.grid()
ax.set_ylabel('Phase Fraction')
ax.set_xlabel('Average Texture Index of Austenite and Ferrite')
ax.set_title('ND sampling, Austenite Only Limit (Jatczak, Hinton)')
ax.set_ylim([-0.05,1.05])
ax.set_xlim([0,10])

#plt.xticks(np.arange(0, 10, 1.0))

#ax.set_ylim([0.15,.35])
#ax.set_xlim([.99,1.5])

plt.savefig("MeetsAusteniteOnly.pdf", dpi=300,format="pdf",orientation='portrait',bbox_inches='tight')
#plt.show()


# ### Other plot attempts

# In[25]:


fig, ax = plt.subplots(figsize=(8,4), dpi=600)

#list of textures
TexNames=list(df2pair.columns.values)

#list of sampling methods
SamplingMethods=list(df2pair[TexNames[0]])

elem=np.arange(len(TexNames[1:]))

#for each sampling method
for i, item in enumerate(df2pair[TexNames[0]]):
    #print i
    #print item
    #Row slice from dataframe, without name
    #print list(df.ix[i])[1:]

    
    #print len(elem),len(list(df1.ix[i])[1:])
    #elem.fill(i)
    if i==0:
        scatter1=ax.scatter(list(df2pair.loc[i])[1:] ,  list(df2pair.loc[i])[1:] ,
                            color='r',label="2Pair")
        #scatter2=ax.scatter(np.full(elem.shape, i) ,  list(df4pair.ix[i])[1:] ,
        #                    color='y',label="4Pair")
        #scatter3=ax.scatter(np.full(elem.shape, i+0.15) ,  list(dfMaxUnique.ix[i])[1:] ,
        #                    color='b',label="MaxUnique")
    else:
        scatter1=ax.scatter(np.full(elem.shape, i-0.15) ,  list(df2pair.loc[i])[1:] ,
                            color='r',label="")
        #scatter2=ax.scatter(np.full(elem.shape, i) ,  list(df4pair.ix[i])[1:] ,
        #                    color='y',label="")
        #scatter3=ax.scatter(np.full(elem.shape, i+0.15) ,  list(dfMaxUnique.ix[i])[1:] ,
        #                    color='b',label="")
        #elem.fill(i)   
    
    

ax.legend(bbox_to_anchor=(0., 1.12, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.)

plt.xticks(np.arange(len(SamplingMethods)),SamplingMethods, rotation=90)
ax.plot([0, 13], [VF, VF], color='k', linestyle='-', linewidth=.5)
plt.grid()
ax.set_ylabel('Phase Fraction')
ax.set_title('Effect of Texture Combinations')

#plt.savefig("TextureCombinations.pdf", dpi=600,format="pdf")
plt.show()


# 
# 
# ## Attempt violin plot of texture combinations
# ### I was able to create the plot, but as this type of plotting is not typically seen, we have chosen to abandon it.

# In[26]:


fig, ax = plt.subplots(figsize=(8,4), dpi=600)

#list of textures
TexNames=list(df2pair.columns.values)
#TexNamestor=TexNames[0]
#TexNames[1:]


#list of sampling methods
SamplingMethods=list(df2pair[TexNames[0]])[:13]
print (SamplingMethods)
elem=np.arange(len(TexNames[1:]))[:13]


#for each sampling method
for i, item in enumerate(SamplingMethods):
    #print i
        print (item)
    #Row slice from dataframe, without name
    #print list(df2pair.ix[i])[1:]
    #print np.arange(len(list(df2pair)))
    

        scatter1=ax.violinplot(list(df2pair.loc[i])[1:],[i-.25] , showmeans=True, showextrema=False, widths=.15)
        for pc in scatter1['bodies']:
            pc.set_facecolor('red')
            pc.set_edgecolor('black')
    
    
        scatter2=ax.violinplot(list(df4pair.loc[i])[1:],[i] , showmeans=True, showextrema=False, widths=.15)
        for pc in scatter2['bodies']:
            pc.set_facecolor('yellow')
            pc.set_edgecolor('black')
    
        scatter3=ax.violinplot(list(dfMaxUnique.loc[i])[1:],[i+.25] , showmeans=True, showextrema=False, widths=.15)
        for pc in scatter3['bodies']:
            pc.set_facecolor('blue')
            pc.set_edgecolor('black')
    
    


ax.legend(bbox_to_anchor=(0., 1.12, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.)

plt.xticks(np.arange(len(SamplingMethods)),SamplingMethods, rotation=90)
ax.plot([0, 13], [VF, VF], color='k', linestyle='-', linewidth=.5)
plt.grid()
ax.set_ylabel('Phase Fraction')
ax.set_title('Effect of Texture Combinations')

plt.savefig("ViolinPlt.png",dpi=600,format="png",orientation='portrait',bbox_inches='tight')

plt.show()


# In[ ]:




