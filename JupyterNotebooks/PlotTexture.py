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
def SingleSchemePlot(Name, Coordinates, Marker, MarkerSize,RD=True,RDup=True,  save=False, cmd=False,savename='test.png'):
    """
    Test plot to work out conventions
    
    Note - mplstereonet plots on the bottom of the sphere
    """
    import pandas as pd
    import mplstereonet
    from matplotlib import pyplot as plt
    
    fig = plt.figure(figsize=(8,9), dpi=600)


    #key
    if RDup==True:
        if RD==True:
            ax1 = fig.add_subplot(111, projection='stereonet')
            
            ax1.set_azimuth_ticks([0,90], labels=['RD','-- TD'],fontsize=14)

            ax1.plane(0.0, 90.0, 'k-', linewidth=1)

            ax1.plane(90.0, 90.0, 'k-', linewidth=1)

            ax1.annotate('ND', xy=(0, 0), xytext=(0.1,0.1), fontsize=14)
        else:
            ax1 = fig.add_subplot(111, projection='stereonet')
            ax1.set_azimuth_ticks([0,90], labels=['X','-- Y'],fontsize=14)
            ax1.plane(0.0, 90.0, 'k-', linewidth=1)
            ax1.plane(90.0, 90.0, 'k-', linewidth=1)
            ax1.annotate('Z', xy=(0, 0), xytext=(0.1,0.1), fontsize=14)
    else:
        if RD==True:
            ax1 = fig.add_subplot(111, projection='stereonet')
            ax1.set_azimuth_ticks([0,270], labels=['RD','TD'],fontsize=14)
            ax1.plane(0.0, 90.0, 'k-', linewidth=1)
            ax1.plane(90.0, 90.0, 'k-', linewidth=1)
            ax1.annotate('ND', xy=(0, 0), xytext=(0.1,0.1), fontsize=14)
            ax1.rotation=-90
            
        else:
            ax1 = fig.add_subplot(111, projection='stereonet')
            ax1.set_azimuth_ticks([0,270], labels=['X','Y'],fontsize=14)
            ax1.plane(0.0, 90.0, 'k-', linewidth=1)
            ax1.plane(90.0, 90.0, 'k-', linewidth=1)
            ax1.annotate('Z', xy=(0, 0), xytext=(0.1,0.1), fontsize=14)
            ax1.rotation=-90
    #SingleOrientation - Name, Tilt, Rotation
    #SchemeName,Coordinates=SingleOrientation("RD Single", 90.0,180.0)
    #SchemeName,Coordinates=SingleOrientation("TD Single", 90.0,270.0)

    # The -90 in rotation is needed to make the X/RD axes consistent.  Otherwise 0° rotation is parallel to Y/TD...
    # NEED to check if this an issue with other plots/summations!
    dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)
    l1=ax1.pole(strike, dip, Marker, markersize=MarkerSize, clip_on=False)
    #dip tilts about the 0 axis (RD), left handed
    #strike tilts about the normal axis (ND), left handed
    #Looking at the bottom of the sphere, not the top
    #this is just convention of mplstereonet, does not affect averaging methods

    if save==True:
        plt.savefig(savename, bbox_inches='tight')

    if cmd==False:
        plt.show()

#####################################
# Plot a Filled Contour with Density of points
#####################################
def DensityContourPlot(Name, Coordinates,RDup=True, Weights=False, save=False, cmd=False, savename='test.png'):
    """
    Uses the density_contourf function of mplstereonet to display how a particular scheme over- or under-samples locations on the pole figure.
    mplstereonet has a hardcoded 1° sampling spacing, the resulting plot may look a little noisy.
    """

    import pandas as pd
    import mplstereonet
    from matplotlib import pyplot as plt
    import numpy as np
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    fig = plt.figure(figsize=(6,6), dpi=300)

    ### Using same colormap as the Matlab/Colormap4, adding alpha in the vector
    n=4
    v = np.linspace(0, n, n*10+1)

    newcolors2=np.array([[0, 0, 0, 1],    [.1, .1, .1, 1], [.2, .2, .2, 1], [.3, .3, .3, 1 ],  [.4, .4, .4, 1], [.5, .5, .5, 1],
             [.6, .6, .6, 1], [.7, .7, .7, 1], [.8, .8, .8, 1], [.9, .9, .9, 1],
             [ 0, 0, 1, 1],   [ 0, .1, .9, 1], [0, .2, .8, 1],  [0, .3, .7, 1],  [0, .4, .6, 1], [0, .5, .5, 1],
             [0, .6, .4, 1],  [0, .7, .3, 1],  [0, .8, .2, 1],  [0, .9, .1, 1],  [0, 1, 0, 1],   [ .1, 1, 0, 1],
             [.2, 1, 0, 1],   [.3, 1, 0, 1],   [.4, 1, 0, 1],   [.5, 1, 0, 1],   [.6, 1, 0, 1],  [.7, 1, 0, 1],
             [.8, 1, 0, 1], [.9, 1, 0, 1],   [1, 1, 0, 1],    [1, .9, 0, 1],   [1, .8, 0, 1],   [1, .7, 0, 1],
             [1, .6, 0, 1], [1, .5, 0, 1],
             [1, .4, 0, 1],   [1, .3, 0, 1],   [1, .2, 0, 1],   [1, .1, 0, 1],   [1, 0, 0, 1]] )
             
    newcmp2 = ListedColormap(newcolors2)


    if RDup==True:
        ax09 = fig.add_subplot(111, projection='stereonet')
        # FIX location
        #ax09.set_azimuth_ticks([0,90], labels=['RD','TD'],fontsize=14)
        ax09.set_azimuth_ticks([90,0], labels=['',''])
    else:
        ax09 = fig.add_subplot(111, projection='stereonet')
        # FIX location
        #ax09.set_azimuth_ticks([0,90], labels=['RD','TD'],fontsize=14)
        ax09.set_azimuth_ticks([90,0], labels=['',''])
        ax09.rotation=-90
    dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)
    
    #Gridsize option? gridsize=500,
    if Weights==False:
        cax = ax09.density_contourf(strike, dip,levels=v, measurement='poles',method='schmidt',cmap=newcmp2 )
    else:
        cax = ax09.density_contourf(strike, dip,levels=v, measurement='poles', method='schmidt',
                                   cmap=newcmp2, weights=list(Coordinates['Weights']) )
    cbar=fig.colorbar(cax,orientation='vertical')
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel("Oversampling Multiples", rotation=270)
    
    if save==True:
        plt.savefig(savename, bbox_inches='tight')

    if cmd==False:
        plt.show()



###################################
# Texture Component Heatmap
###################################
def PlotHeatmap(hw, PeakCombo,Scheme, Folder, Scheme2=None, VF=0.25, cbarMap=False, cbarRange=[0,0.5], save=False, cmd=False, savename='test.png'):
    """
    A method that accepts a desired HalfWidth value, Peak Combination, and a specific Sampling Scheme, and outputs a heatmap of the 
    resulting Austenite Phase Fraction Values. In order for the function to work as intended, please make sure you note the format 
    of the parameters and their specifications!
    
    Parameters
    ----------
    
    hw: int
    Desired HalfWidth Value. Halfwidths should range from 5 to 50, increasing by multiples of 5 (check me on 
    this).
    
    
    PeakCombo: String
    Desired Peak Combination to study. Peak Combinations fall under DF2, DF4, and DFUnique, each of which refer to the specific combinations of XRD peaks used in calculations. Must be defined name from TextureSampling.py
    
    Scheme: String
    Desired Sampling Scheme. Make sure the format of the inputted scheme exactly matches that which is prescribed by the original
    data. For example, to view "HexGrid-60degTilt5degRes" sampling scheme, you need to type that exactly for the code to work.
    
    Folder: string
    Folder to save the data
    
    Scheme2: String (default None)
    When using different two theta values, include data from a second scheme.  Enter the Austenite phase for Scheme2, Ferrite phase for Scheme
    
    
    """
    import fnmatch
    import os
    import pandas as pd
    import numpy as np
    import seaborn as sns
    from matplotlib import pyplot as plt
    HW=str(hw)

    
    AusteniteTextures=[]
    FerriteTextures=[]  


    for file in os.listdir(Folder):
    
        if (file.find("A-HW"+HW)>0):
            AusteniteTextures.append(file)
        elif (file.find("F-HW"+HW)>0):
            FerriteTextures.append(file)
        else: ()
        
        #print (FerriteTextures)
        #print AusteniteTextures
    
    # Create as dictionary, easier for seaborn later, and single valued entries
    AustList=[]
    FerrList=[]
    VFList=[]
    
    for AustOrient in AusteniteTextures:
        for FerrOrient in FerriteTextures:
        
            
            DFF=pd.read_excel(os.path.join(Folder,FerrOrient),header=1,skipfooter=0)
            DFA=pd.read_excel(os.path.join(Folder,AustOrient),header=1,skipfooter=0)

            if Scheme2 != None:
            
                F_val=float(DFF.loc[DFF['HKL'] == Scheme][PeakCombo])
                A_val=float(DFA.loc[DFA['HKL'] == Scheme2][PeakCombo])
            
            else :
                # Revised, some issues with what's returned
                # Pandas returns a larger array instead of just one value - AC 2021 Mar 29
                F_val=float(DFF.loc[DFF['HKL'] == Scheme][PeakCombo])
                A_val=float(DFA.loc[DFA['HKL'] == Scheme][PeakCombo])
                #print (DFF)
                #print (DFA)
                
            #split at file type and halfwidth for names
            Fname=FerrOrient.split(".")[0].split("-")[0]
            Aname=AustOrient.split(".")[0].split("-")[0]
            
            # Append to lists
            AustList.append(Aname)
            FerrList.append(Fname)
            VFList.append(VF*A_val/(VF*A_val+((1.0-VF)*F_val)))
        
        
    dataDict={"Austenite":AustList, "Ferrite": FerrList, "VF":VFList}
    SaveTable=pd.pivot_table(pd.DataFrame.from_dict(dataDict),index="Austenite", columns="Ferrite")
    
    #used for debugging
    #return dataDict
    #return EmptyDF
    
    
    # Need to rename columns since pivot_tables return a multiindex data frame
    FerriteNames=[col[1] for col in SaveTable.columns.values]
 
    # I think this has been depricated- AC 29 Mar 2021
#    values[values == ''] = 0.0
#    Values = values.astype(np.float)
#    labels=[]
#    for val in Values:
#        if(val>=0.255):
#            labels.append("+")
#        elif(val<=0.245):
#            labels.append("-")
#        else:
#            labels.append("O")
#    labels=np.asarray(labels)
    #labels.resize(len(yaxis),len(xaxis)) 
    #labels=pd.DataFrame(labels)

    # plotting
    
    if cbarMap=='grey':
        color= sns.diverging_palette(359, 359, 99, l=0, sep=1, n=50, center='light', as_cmap=True)
    else:
        color= sns.color_palette("coolwarm", 25)
    plt.figure(figsize = (13,7))
    figure=sns.heatmap(SaveTable, vmin=cbarRange[0], vmax=cbarRange[1], cmap=color, center=VF, annot=True, fmt="3.3f", linewidths=0.5,square=True,cbar_kws={"shrink": .80}, xticklabels=FerriteNames)
    figure.set_xlabel('Ferrite')
    
    #figure=sns.heatmap(df_wide,vmin=0.0, vmax=0.50, cmap=color,center=0.25,annot=dw, annot_kws={"size": 18},fmt='',linewidths=0.5,square=True,cbar_kws={"shrink": .80})
    plt.title("Halfwidth of "+HW+" , "+ Scheme + " Sampling Scheme, "+ PeakCombo.upper()+ " Peak Combination" ,fontsize =18)
    bottom, top = figure.get_ylim()
    figure.set_ylim(bottom + 0.5, top - 0.5)
    
    #return figure
    if save==True:
        plt.savefig(savename, bbox_inches='tight')

    if cmd==False:
        plt.show()



