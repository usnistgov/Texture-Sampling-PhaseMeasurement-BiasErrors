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
