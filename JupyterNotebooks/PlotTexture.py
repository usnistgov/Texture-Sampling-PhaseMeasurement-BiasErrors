#####################################
# ODF Key
#####################################
def ODFKey(save=False,cmd=False,savename='test.png'):
    """
    Figure 2: Phi2 = 45° ODF Key for Cubic-Orthotropic Symmetry

    Figure 2 from J App Cryst paper
    """
    import matplotlib.pyplot as plt

    fig=plt.figure(figsize=(4,3), dpi=600)
    ax1 = fig.add_subplot(111)


    labels=[0,15,30,45,60,75,90]
    ax1.set_xlim(0,90)
    ax1.set_xlabel(r'$\phi_1$', fontsize=20)
    ax1.xaxis.set_tick_params(labelsize=16)
    plt.xticks(labels, labels)

    ax1.set_ylim(-0,90)
    ax1.set_ylabel(r'$\Phi$', fontsize=20)
    ax1.set_ylim(ax1.get_ylim()[::-1]) #mirror for conventional plot
    ax1.yaxis.set_tick_params(labelsize=16)
    plt.yticks(labels, labels)


    #legend=ax1.legend(loc='upper right',fontsize=10)


    ax1.scatter([45],[0], label='Cube',color='k', marker='s', clip_on=False,s=100.0,zorder=10)
    ax1.scatter([0,90],[0,0], label='Shear',color='k', marker='D', clip_on=False,s=100.0,zorder=10)

    ax1.scatter([90],[90], label='Goss',color='#9400D3', marker='o', clip_on=False,s=150.0,zorder=10)
    ax1.scatter([0],[90], label='Rot. Goss', facecolors='none', edgecolors='#9400D3', marker='o', clip_on=False,s=150.0,zorder=10)
    ax1.scatter([54.74],[90], label='Brass',color='g', marker='^', clip_on=False,s=150.0,zorder=10)
    ax1.scatter([90],[35.26], label='Copper',color='g', marker='v', clip_on=False,s=150.0,zorder=10)

    ax1.scatter([0],[15.79], label='alpha1',color='r', marker='1', clip_on=False,s=200.0,zorder=10)
    ax1.scatter([0],[25.24], label='alpha2',color='r', marker='2', clip_on=False,s=200.0,zorder=10)
    ax1.scatter([0],[35.26], label='alpha3',color='r', marker='3', clip_on=False,s=200.0,zorder=10)
    ax1.scatter([0],[43.31], label='alpha4',color='r', marker='4', clip_on=False,s=200.0,zorder=10)

    ax1.scatter([0,60],[54.7,54.7], label='gamma1',color='b', marker='x', clip_on=False,s=150.0,zorder=10)
    ax1.scatter([30,90],[54.7,54.7], label='gamma2',color='b', marker='+', clip_on=False,s=200.0,zorder=10)
    ax1.scatter([90],[60.5], label='554',color='#20B2AA', marker='p', clip_on=False,s=200.0,zorder=10)

    ax1.plot([0,90],[54.7,54.7],label='gamma fiber',color='b', clip_on=False,linewidth=1, zorder=2)
    ax1.plot([0,0],[0,90],label='alpha fiber',color='r', clip_on=False,linewidth=1, zorder=9)
    #ax1.plot(TRIP700Data[:,0],TRIP700Data[:,1], label='',color='m')

    ax1.set_axisbelow(True)
    ax1.grid(color='#9F9F9F', linestyle='-', linewidth=0.5)

    #ax1.set_xlim([40, 145])

    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0, box.width * 0.75, box.height])

    # Put a legend to the right of the current axis
    ax1.legend(loc='center left', bbox_to_anchor=(1.1, 0.5),fontsize=16,ncol=2)


    if save==True:
        plt.savefig(savename, bbox_inches='tight')

    if cmd==False:
        plt.show()

#####################################
# Simple Test Plot
#####################################
def SimpleTestPlot(Name, Coordinates, save=False,cmd=False,savename='test.png'):
    """
    Test plot to work out conventions
    """
    import pandas as pd
    import mplstereonet
    from matplotlib import pyplot as plt
    
    fig = plt.figure(figsize=(8,8), dpi=600)

    ax2 = fig.add_subplot(111, projection='stereonet')

    #SingleOrientation - Name, Tilt, Rotation
    #SchemeName,Coordinates=SingleOrientation("RD Single", 90.0,180.0)
    #SchemeName,Coordinates=SingleOrientation("TD Single", 90.0,270.0)

    dip, strike =Coordinates['Tilt'], Coordinates['Rotation']-90.0
    l1=ax2.pole(strike, dip, 'bD', markersize=10, clip_on=False)
    #dip tilts about the 0 axis (RD), left handed
    #strike tilts about the normal axis (ND), left handed
    #Looking at the bottom of the sphere, not the top
    #this is just convention of mplstereonet, does not affect averaging methods

    if save==True:
        plt.savefig(savename, bbox_inches='tight')

    if cmd==False:
        plt.show()


#####################################
# Plot a single sampling scheme
#####################################
def SingleSchemePlot(Name, Coordinates, Marker, Markersize, save=False, cmd=False, savename='test.png'):
    """
    Test plot to work out conventions
    """
    import pandas as pd
    import mplstereonet
    from matplotlib import pyplot as plt
    
    fig = plt.figure(figsize=(8,9), dpi=600)


    #key
    ax1 = fig.add_subplot(111, projection='stereonet')
    ax1.set_azimuth_ticks([0,90], labels=['RD','-- TD'],fontsize=14)
    ax1.plane(0.0, 90.0, 'k-', linewidth=1)
    ax1.plane(90.0, 90.0, 'k-', linewidth=1)
    ax1.annotate('ND', xy=(0, 0), xytext=(0.1,0.1), fontsize=14)
        
    
    #SingleOrientation - Name, Tilt, Rotation
    #SchemeName,Coordinates=SingleOrientation("RD Single", 90.0,180.0)
    #SchemeName,Coordinates=SingleOrientation("TD Single", 90.0,270.0)

    dip, strike =Coordinates['Tilt'], Coordinates['Rotation']-90.0
    l1=ax1.pole(strike, dip, Marker, markersize=10, clip_on=False)
    #dip tilts about the 0 axis (RD), left handed
    #strike tilts about the normal axis (ND), left handed
    #Looking at the bottom of the sphere, not the top
    #this is just convention of mplstereonet, does not affect averaging methods

    if save==True:
        plt.savefig(savename, bbox_inches='tight')

    if cmd==False:
        plt.show()





###################################
# Texture Component Heatmap
###################################
def PlotHeatmap(hw,PeakCombo,Scheme):
    #def PlotHeatmap(Directory, HW, PeakCombination, SamplingScheme, Range):
    """
    Read list of files from directory and then operate on them
    """

    import numpy as np
    import pandas as pd
    import scipy as scipy
    from scipy import interpolate
    from scipy import signal
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter
    #import mplstereonet
    import math
    import os
    import glob
    import statistics



    # Get the current working directory path
    cwd=os.getcwd()
    #print cwd
    xpcdatapath=os.path.abspath(os.path.join(os.path.dirname( cwd)))
    #print xpcdatapath
    # Use join instead, it's better cross platform
    Folder=os.path.join(cwd,'AveragedIntensites')



    import fnmatch
    HW=str(hw)

    VF=0.25 #Austenite Volume Fraction


    #Read in Texture files matching HW of 20 degrees
    AusteniteTextures=[]
    FerriteTextures=[]  


    for file in os.listdir(Folder):
    
        if (file.find("A-HW"+HW)>0):
            AusteniteTextures.append(file)
        elif (file.find("F-HW"+HW)>0):
            FerriteTextures.append(file)
        else: ()
        
        #print FerriteTextures
        #print AusteniteTextures

    # create a dataframe shape from existing data
    DFA=pd.read_excel((os.path.join(Folder, FerriteTextures[0])),header=1,skip_footer=0)
    

    #copy the HKL reflection indexes
    df2pair=DFA["HKL"]
    df4pair=DFA["HKL"]
    dfMaxUnique=DFA["HKL"]

    for AustOrient in AusteniteTextures:
        for FerrOrient in FerriteTextures:
        

            DFF=pd.read_excel(os.path.join(Folder,FerrOrient),header=1,skip_footer=0)
            DFA=pd.read_excel(os.path.join(Folder,AustOrient),header=1,skip_footer=0)

            #Switch position of Tilt -row 7 and Rotate -row 8 data rows, matches Figure 4 better
            #http://stackoverflow.com/questions/32929927/pandas-swap-rows-between-dataframes
            tempF=DFF.loc[8]
            DFF.loc[8,:]=DFF.loc[7,:].values
            DFF.loc[7,:]=tempF.values

            tempA=DFA.loc[8]
            DFA.loc[8,:]=DFA.loc[7,:].values
            DFA.loc[7,:]=tempA.values
        
            #print (DFF)
            #print (DFA)
        
            # add minus VF for plotting
            DF1=(VF*DFA["2Pairs"]/(VF*DFA["2Pairs"]+((1.0-VF)*DFF["2Pairs"])))#-VF
            DF2=(VF*DFA["4Pairs"]/(VF*DFA["4Pairs"]+((1.0-VF)*DFF["4Pairs"])))#-VF
            DF3=(VF*DFA["MaxUnique"]/(VF*DFA["MaxUnique"]+((1.0-VF)*DFF["MaxUnique"])))#-VF
        
            Fname=FerrOrient.split(".")
            Aname=AustOrient.split(".")
        
            MixName= Fname[0]+"-"+Aname[0]
        
            
        
            data1 = pd.DataFrame({MixName: DF1})
            data2 = pd.DataFrame({MixName: DF2})
            data3 = pd.DataFrame({MixName: DF3})
        
            
        
            df2pair = pd.concat([df2pair, data1], axis=1)
            df4pair = pd.concat([df4pair, data2], axis=1)
            dfMaxUnique = pd.concat([dfMaxUnique, data3], axis=1)
        
    if(PeakCombo.lower()=='df2'):
        data=df2pair
    elif (PeakCombo.lower()=='df4'):
        data=df4pair
    else:
        data=dfMaxUnique
    
    
    i=0
    for i in range (len(data.iloc[:,0])):
        if (data.iloc[:,0][i].lower()==Scheme.lower()):
            break
    
    relevantdata=data.iloc[i]
    print(relevantdata) #relevantdata is the dataframe of interest here



    print(PlotHeatmap(20,"df4","ND Single"))  #must type df2 for df2pair and df4 for df4pair; type anything else you want
    #and you get dfMaxUnique!
       
    
    ## Reshape list to a matrix

    ## plot matrix as a heatmap
    
    ## need to end the definition
    return






    # The following code is referenced from the following website: https://towardsdatascience.com/better-heatmaps-and-correlation-matrix-plots-in-python-41445d0f2bec
    # Multiple spaces between code segments shows these are carried out in different cells in Jupyter Notebook
    # Cell 1- Import all necessary libraries and packages;

#    import pandas as pd
#    from matplotlib import pyplot as plt
#    from pylab import rcParams
#    rcParams['figure.figsize'] = 7,7
#    import seaborn as sns
#    import numpy as np
#    sns.set(color_codes=True, font_scale=0.8)

#    #AC - not sure if these belong here
#    %matplotlib inline
#    %config InlineBackend.figure_format = 'retina'
#    %load_ext autoreload
#    %autoreload 2
#
#    #Cell 2- The most important import to carry out this function!
#
#    #AC - pip commands belong elsewhere. We should work out a conda env.
#    !pip install heatmapz
#
#
#    #Cell 3- Import specific methods from Heatmap library
#
#    # Import the two methods from heatmap library
#    from heatmap import heatmap, corrplot
#
#    #Cell 4- Read CSV File (COULD BE CHANGED)
#
#    data = pd.read_csv(#https://raw.githubusercontent.com/drazenz/heatmap/master/autos.clean.csv'#) #need to replace the CSV to be generated
#
#
#    #Debating whether or not to accept the databases simply as parameters and go from there?
#
#
#
#    #Cell 5- Final configurations for plot size and display
#
#    plt.figure(figsize=(8, 8))
#    corrplot(data.corr(), size_scale=500);

    #Can change the size of the figure (right now it is 8 inches by 8 inches), or can change how the squares in the heatmap are sized relative to the heatmap (size scale)
