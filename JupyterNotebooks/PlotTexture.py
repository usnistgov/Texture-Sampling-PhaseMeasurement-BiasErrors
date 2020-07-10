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
def PlotHeatmap(Directory, HW, PeakCombination, SamplingScheme, Range):
    """
    Read list of files from directory and then operate on them
    """


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

    ### Example pseudo code ###
    
    ## Use pandas mask or matching to subset the data file to only specific HW
    
    ## Further select a specific peak combination
    
    ## Further selection a specific sampling scheme
    
    ## Reshape list to a matrix

    ## plot matrix as a heatmap
    
    ## need to end the definition
    return
