#####################################
#####################################
# Begin Operations on Intensities
#####################################
#####################################

#####################################
# Read in Polefigures in .xpc format
#####################################
import fortranformat as ff
import numpy as np
import pandas as pd
import math
def xpcformat(mode=None, filename=None):
    
    """
    
    FIX: Docstring too long and complicated
    
    NOTE: .xpc files seem to start (rotation=0) along the Y (TD) direction. AC-2021 Apr 14
    
    MAUD uses an .xpc format for pole figures, likely derived from BearTex [1].  This format is similar to the General Intensity File Format in POPLA [2, appendix B2], with a slightly different header

    [1] http://eps.berkeley.edu/~wenk/TexturePage/beartex.htm

    [2] Popla Manual http://pajarito.materials.cmu.edu/rollett/27750/popLA_Manual.pdf


    #        Adapted .xpc format parser, from
    #        https://github.com/usnistgov/texture
    #        commit 9c0ac85
    #        Program upf.py
    # direct link: https://github.com/usnistgov/texture/blob/9c0ac8531276a5d8d27c0e895074ca66fda76608/src/upf.py

    # .xpc format
    # line with a phase, reflection and trailing '#'
    # 4 blank lines
    # 1 line with lattice parameters
    # 1 line with reflection and pole figure explaination
    # 4 lines, with 18 terms each, integers - fixed width (fortran format)
        
    
    Experimental pole figure format controller
    mode:
  ->    "steglich"
  ->    "bruker"  *.uxd file
  ->    "epf" (2011-Oct-6) epf popLA experimental pole figure format
    Returns the pole figure data as
  ->  the standard format (m x n numpy array)
  ->  each of axes stands for rotating (phi) and tilting (khi)
  ->  angle in the laboratory space.
    These angles will be angle and radius in the
    space onto which a pole figure is projected.
    conventions:
     tilting angle: khi
     rotating angle: phi
     dk: incremental khi angle
     dp: incremental phi angle
     nk: number of points along khi axis
     np: number of points along phi axis
     angles are converted into radian whenever possible
    Arguments
    =========
    mode     = None
    filename = None
    """
   
    
    print ('Pole Figure Parsing')

    if mode=='xpc':
        """
        Adapted .xpc format parser, from
        https://github.com/usnistgov/texture
        commit 9c0ac85
        based on ready made popLA epf format parser
        consider the possibility of multiple number of polefigure
        ## phi must be 0~355, khi must be 0~90 with 5 as an ang
        resolution
        """
        print ('You are now reading experimental pole figure(s): \n\t%s'%filename)
        blocks = open(filename, 'r').read().split('\n\n\n\n')[1:]
        #blocks = open(filename, 'rU').read().split('\n\n\n\n')[1:]   #changed, U mode deprecated
        print ('There are %s blocks of data found'%len(blocks))
        if len(blocks)==0:
            msg1 = 'xpc parser in upf assumes that pole figures are separated by 4 new lines'
            msg2 = ' searching %s finds no set of 4 new lines in '%filename
            msg  = '%s \n %s'%(msg1,msg2)
            raise IOError (msg)
            # blocks = parse_epf(filename)
        npf = len(blocks)
        if npf==0: raise IOError ('No pf block found.')

        datasets = []; max_khi = []
        if  npf>1: hkls=["HKL"] ## multiple number of pole figures in a file
        
        for part in blocks:
            line=part.split('\n')
            #print len(line)

            structureline=ff.FortranRecordReader('(6f10.4,1x,i4,1x,i4)')
            [a,b,c,alpha,beta,gamma,crystalclass,something]=structureline.read(line[1])
            pfDefline=ff.FortranRecordReader('(1x,3i3,6f5.1,2i2)')
            [h,k,l,unknown1,tilt,tiltinc,unknown2,rotation,rotationinc,unknown3,unknown4]=pfDefline.read(line[2])
            
            #for the rest of the lines, do the following
            dataline=ff.FortranRecordReader('(1x,18i4)')
            
            # Pretty ugly code, but works...
            grouping=[[3,4,5,6],[7,8,9,10],[11,12,13,14],[15,16,17,18],[19,20,21,22],[23,24,25,26],
                      [27,28,29,30],[31,32,33,34],[35,36,37,38],[39,40,41,42],[43,44,45,46],[47,48,49,50],
                      [51,52,53,54],[55,56,57,58],[59,60,61,62],[63,64,65,66],[67,68,69,70],[71,72,73,74],
                     [75,76,77,78]]
            
            
            dataset=[]
            for item in grouping:
                #print item[0],item[1],item[2],item[3]
                # Changed the order to deal with 90 rotation in xpc files AC-2021 Apr 14
                parsed=dataline.read(line[item[1]])
                parsed.extend(dataline.read(line[item[2]]))
                parsed.extend(dataline.read(line[item[3]]))
                parsed.extend(dataline.read(line[item[0]]))
                dataset.append(parsed)
            #print (dataset)
            
            #Saves as a Pandas dataframe, and maps the 360 degree phi data from the 0 degree phi data
            #row and column indexes are by degrees
            
            # FIX - changed in new version of pandas
            df=pd.DataFrame(dataset, index=np.arange(0,91,5))
            df.columns=[np.arange(0,360,5)]
            df[360]=df.iloc[:,0]  #tried changing .loc to .iloc
            
            # Save the hkl value
            hkl = [h,k,l] #hkl
            #print hkl
            hkls.append(hkl)

            datasets.append(df)
            
        #print hkls
        print("number of pole figures:", len(datasets))

        return datasets, hkls
    else: raise IOError ('Unexpected mode is given')
    #return data


#####################################
# Read in Polefigures in mtex format
#####################################
def mtexPFformat(filepath=None):
    """
    Read in pole figures exported from mtex using the 'export' command
    Each pole figure is a separate file with the hkl given in the filename inside paratheses
    This reader assumes that the collection of pole figures are the only .txt files in the filepath given
    
    
    Note: have not checked rotation issue noted in .xpc files... AC-2021 Apr 14
    
    conventions:
    
    Input Arguments
    =========
    filepath : str
       File path to the folder with the pole figures

  
    Output Arguments
    =========
    PFs
    hkllist
    """

    import fortranformat as ff
    import numpy as np
    import pandas as pd
    import math
    import os
    import glob
    
    print ('Pole Figure Parsing')
    # read list of files, assumes they are saved as .txt and are the only files in directory
    PFFilenamesList = glob.glob(os.path.join(filepath, '*.txt'))
    
    # dataset is a list of pandas dataframes
    #pf data format is a pandas array with tilt and rotate as the indexes


    # Initialize list to add to the dataframes
    hkls=["HKL"]
    datasets=[]
    for PFFile in PFFilenamesList:
        (head,tail)=os.path.split(PFFile)
        print('You are now reading experimental pole figure(s): \n\t%s'%tail)
        hkl_name=tail.split("(")[1].split(")")[0]
        #print("HKL: ",hkl)

        # mtex data is 3 columns: tilt, rotate, intensity
        #CHECK - tilt convention
        OriginalData=pd.read_csv(PFFile,sep='\s+',names=["Tilt", "Rotate","Intensity"])
        #Convert to Ints
        #2021 June 08 ??? Why are limiting to Ints
        OriginalData["Tilt"]=OriginalData["Tilt"].astype(int)
        OriginalData["Rotate"]=OriginalData["Rotate"].astype(int)
        #print(OriginalData)
        
        
        MatrixData=OriginalData.pivot_table(index="Tilt", columns="Rotate", values="Intensity")
        # replace the 0 tilt row with 180 row, since mtex keeps the rotation all at zero
        MatrixData=MatrixData.fillna(value=MatrixData[0].iloc[0])
        # add a 360 point for interpolation
        MatrixData[360]=MatrixData[0]
        
        #print(MatrixData)
        #Adjust intensities for similarity to Beartex format
        MatrixData=MatrixData*100

        # Save the hkl value
        hkl = [int(hkl_name[0]),int(hkl_name[1]),int(hkl_name[2])] #hkl
        #print hkl
        hkls.append(hkl)
#
        datasets.append(MatrixData)
#
#    #print hkls
#    print("number of pole figures:", len(datasets))

    # return the same structure as the .xpc reader
    return datasets, hkls




#####################################
# Calculate Intensity from Pole Figure Coordinates
#####################################
def pfIntensitySum(name, PoleFigures, Coordinates, Weights=False):

    """
    Define function to take a pole figure (or series of pole figures) and series of coordinate pairs, returning an average intensity for all of the pole figures and coordinates
    For each coordinate pair passed to this function, interpolate the pole figure to find intensity
    Average intenisity for all coordinate pairs
    Average intensity for all pole figures passed
    
    conventions:

    Input Arguments
    =========
    name: passed through program to help mark columns
    PoleFigures - series of
    -> FORMAT
    Coordinates
    -> FORMAT
    array of coordinates (rotation,tilt)
  
    Output Arguments
    =========
    AverageIntensity Array for each pole figure with average intensity
    """
    import numpy as np
    import pandas as pd
    import scipy as scipy
    from scipy import interpolate
    
    #print "Averaging Intensity from Pole Figures"
    AverageIntensity=[name]

    ## For each pole figure:
    for pf in PoleFigures:
        IntensityValues=[]
        #print "Pole Figure Data:"
        #print len(list(pf.columns.values))
        #print len(list(pf.index.values))
        #print pf

        ## Interpolate the pole figure
        #set kx,ky=1 for linear interpolation, otherwise got odd behavior near zero on edges
        InterpPF=scipy.interpolate.RectBivariateSpline(list(pf.index.values),list(pf.columns.values), pf, kx=1, ky=1)

        ## For each coordinate:
        ## Read the value from the pole figure, append to new array
        #print Coordinates
        
        if Weights==True:
            for index, row in Coordinates.iterrows():
                #print row['Tilt'], row['Rotation'], InterpPF.ev(row['Tilt'],row['Rotation'])
                IntensityValues.append(InterpPF.ev(row['Tilt'],row['Rotation'])*row['Weights'])
        else:
            for index, row in Coordinates.iterrows():
                #print row['Tilt'], row['Rotation'], InterpPF.ev(row['Tilt'],row['Rotation'])
                IntensityValues.append(InterpPF.ev(row['Tilt'],row['Rotation']))
        
        #print IntensityValues
        #Factor of 100 is divided to convert from POPLA style format of 100 = 1 MRD/MUD
        if Weights==True:
            AverageIntensity.append(sum(IntensityValues)/(100))
            # Length of intensity values included in the weighting factor, so it doesn't need to be included here
        else:
            AverageIntensity.append(sum(IntensityValues)/(100*len(IntensityValues)))

    #print AverageIntensity

    ## return average value
    return AverageIntensity


def GenerateAveIntesity(SchemesListDF, pftype, DataFolder, SaveFolder):
    """
    Generate average intensity based on pole figures and coordinates

        This section calculates the average intensity and saves to file
        Looks for the list of XPC files and calculates a table (.xlsx) for each

        SchemesList: Pandas dataframe of schemes
        pftype: mtex (.txt) or Beartex (.xpc) format
        DataFolder: Where is the data
        SaveFolder: Where to save the .xlsx files
        
        TO ADD - list of peak combinations, citations
        
        TODO - allow for subsets instead of all of them in a folder?
        TODO - reference by HKL not column name
    """
    import os
    import pandas as pd
    import numpy as np
    import glob
    #SchemeName,Coordinates=

    # Get the current working directory path
    #cwd=os.getcwd()
    #print cwd
    #xpcdatapath=os.path.abspath(os.path.join(os.path.dirname( cwd)))
    #print xpcdatapath
    #Folder=xpcdatapath+'/MAUD/XPCFiles'
    #Folder='C:\Research\Texture-Sampling-PhaseMeasurement-BiasErrors-master\MAUD\XPCFiles'
    #SaveFolder="AveragedIntensites"  Temporary fix
    #SaveFolder='C:\Research\Texture-Sampling-PhaseMeasurement-BiasErrors-master\JupyterNotebooks\AveragedIntensites'

    if not os.path.isdir(SaveFolder):
        os.makedirs(SaveFolder)
        
    os.chdir(SaveFolder)

    # XPC files have all the pole figures in single files

    #for file in os.listdir(DataFolder):
    for file in glob.glob(os.path.join(DataFolder, '*')):
        print("File from listdir: ",file)
        print("pftype: ", pftype)
        if pftype=="xpc":
            if file.endswith(".xpc"):
                XPCfile=(os.path.join(DataFolder, file))
                (head,tail)=os.path.split(file)
                #Split for HW
                if "-" in tail:
                    orientation, hw=tail.split('-')
                else:
                    orientation, ext=tail.split('.')
                    
                PhaseType= orientation[-1:]
                
                #for XPCfile in listoffiles:
                (pfs,hkllist)=xpcformat('xpc',XPCfile)

        # Mtex files have all the pole figures as separate files
        elif pftype=="mtex":
            MtexFolder=(os.path.join(DataFolder, file))
            (head,tail)=os.path.split(MtexFolder)
            print("Mtex data from: ",tail)
            # Look into sub directory
            if "-" in tail:
                orientation, hw=tail.split('-')
            else:
                orientation, ext=tail.split('.')
                
            PhaseType= orientation[-1:]
            
            #for XPCfile in listoffiles:
            (pfs,hkllist)=mtexPFformat(MtexFolder)
                
        else:
            print("Not a suppored type")


            #create subsets for phase fractions
        hkllist.append('1Pair-A')
        hkllist.append('1Pair-B')
        hkllist.append('2Pairs-A')
        hkllist.append('2Pairs-B')
        hkllist.append('3Pairs-A')
        hkllist.append('3Pairs-B')
        hkllist.append('3Pairs-C')
        hkllist.append('4Pairs')
        hkllist.append('5A4F')
        hkllist.append('MaxUnique')

        OutputList=[hkllist]
        #print(OutputList)
        for Scheme in SchemesListDF["SchemeName"]:
            #print(Scheme)
            
            #referencing a dataframe in a dataframe seems difficult in pandas.
            # this is a workaround
            DFsubset=SchemesListDF["Coordinates"][SchemesListDF["SchemeName"] == Scheme]
            
            if "Weight" in Scheme:
                PfIS=pfIntensitySum(Scheme,pfs, DFsubset.iloc[0], Weights=True)
                print("Weighting the sampling scheme")
            else:
                # df.loc[df['column_name'] == some_value]
                PfIS=pfIntensitySum(Scheme,pfs, DFsubset.iloc[0] )
            
            # TO DO - fix what columns are read in
            
            # need to check order of hkl and adjust
            # Glob order is kind of random when returning a file list...

            if PhaseType=="A":
                #Austenite Reflections
                #[4, 2, 0]	[4, 0, 0]	[2, 0, 0]	[2, 2, 0]	[3, 3, 3]	[1, 1, 1]	[2, 2, 2]	[3, 3, 1]	[3, 1, 1]	[4, 2, 2]	[5, 1, 1]
                Aindex420=hkllist.index([4, 2, 0])
                Aindex400=hkllist.index([4, 0, 0])
                Aindex200=hkllist.index([2, 0, 0])
                Aindex220=hkllist.index([2, 2, 0])
                Aindex333=hkllist.index([3, 3, 3])
                Aindex111=hkllist.index([1, 1, 1])
                Aindex222=hkllist.index([2, 2, 2])
                Aindex331=hkllist.index([3, 3, 1])
                Aindex311=hkllist.index([3, 1, 1])
                Aindex422=hkllist.index([4, 2, 2])
                Aindex511=hkllist.index([5, 1, 1])
            elif PhaseType=="F" or "M":
                #Ferrite/Martensite Reflections
                #[2, 1, 1]	[4, 0, 0]	[2, 0, 0]	[1, 1, 0]	[2, 2, 2]	[3, 1, 0]	[3, 2, 1] [2, 0, 0]
                Mindex211=hkllist.index([2, 1, 1])
                Mindex400=hkllist.index([4, 0, 0])
                Mindex200=hkllist.index([2, 0, 0])
                Mindex110=hkllist.index([1, 1, 0])
                Mindex220=hkllist.index([2, 2, 0])
                Mindex222=hkllist.index([2, 2, 2])
                Mindex310=hkllist.index([3, 1, 0])
                Mindex321=hkllist.index([3, 2, 1])
            else:
                print ("Unrecognized Phase")
            #####################
            # 1 Pair A:  Austenite (111);Ferrite (110)
            # Used in Chess data https://doi.org/10.1016/j.msea.2019.05.017
            if PhaseType=="A":
                PfIS.append(PfIS[Aindex111])
            elif PhaseType=="F" or "M":
                PfIS.append(PfIS[Mindex110])
            else:
                print ("Unrecognized Phase")

            #####################
            # 1 Pair B:  Austenite (220);Ferrite (211)
            # Back reflection from a Cr source
            if PhaseType=="A":
                PfIS.append(PfIS[Aindex220])
            elif PhaseType=="F" or "M":
                PfIS.append(PfIS[Mindex211])
            else:
                print ("Unrecognized Phase")
                
                
            # 2 Pairs A:  Austenite (200), (220);Ferrite (200), (211)
            # Used in JAC paper, matches ASTM E975 with Chromium radiation, Jacques 2009 Round Robin (XRD3)
            if PhaseType=="A":
                PfIS.append(np.mean([PfIS[Aindex200],PfIS[Aindex220]]))
            elif PhaseType=="F" or "M":
                PfIS.append(np.mean([PfIS[Mindex200],PfIS[Mindex211]]))
            else:
                print ("Unrecognized Phase")
            
            # 2 Pairs B:  Austenite (220), (311);Ferrite (200), (211)
            # Jacques 2009 Round Robin (XRD5)
            if PhaseType=="A":
                PfIS.append(np.mean([PfIS[Aindex220],PfIS[Aindex311]]))
            elif PhaseType=="F" or "M":
                PfIS.append(np.mean([PfIS[Mindex200],PfIS[Mindex211]]))
            else:
                print ("Unrecognized Phase")
                
            #####################

            # 3 Pairs A:  Austenite (200), (220), (222); Ferrite (200), (211), (310)
            # Used in DXC 2020 presentation
            if PhaseType=="A":
                PfIS.append(np.mean([PfIS[Aindex200],PfIS[Aindex220],PfIS[Aindex222]]))
            elif PhaseType=="F" or "M":
                PfIS.append(np.mean([PfIS[Mindex200],PfIS[Mindex211],PfIS[Mindex310]]))
            else:
                print ("Unrecognized Phase")
            
            # 3 Pairs B:  Austenite (200), (220), (311); Ferrite (200), (211), (220)
            # Skip A111/F110 and A222 (weak peak)
            if PhaseType=="A":
                PfIS.append(np.mean([PfIS[Aindex200],PfIS[Aindex220],PfIS[Aindex311]]))
            elif PhaseType=="F" or "M":
                PfIS.append(np.mean([PfIS[Mindex200],PfIS[Mindex211],PfIS[Mindex220]]))
                #used to include 4, exclude 5, but that is incorrect
            else:
                print ("Unrecognized Phase")

            # 3 Pairs C:   Austenite (111), (200), (220); Ferrite (110), (200), (211)
            # Jacques 2009 Round Robin (XRD2)
            if PhaseType=="A":
                PfIS.append(np.mean([PfIS[Aindex111],PfIS[Aindex200],PfIS[Aindex220]]))
            elif PhaseType=="F" or "M":
                PfIS.append(np.mean([PfIS[Mindex110],PfIS[Mindex200],PfIS[Mindex211]]))
                #used to include 4, exclude 5, but that is incorrect
            else:
                print ("Unrecognized Phase")

            #####################

            # 4 Pairs A:  Austenite (111), (200), (220), (311); Ferrite (110), (200), (211), (220)
            # Used in JAC paper, Jacques 2009 Round Robin (XRD1 and XRD 4)
            if PhaseType=="A":
                PfIS.append(np.mean([PfIS[Aindex111],PfIS[Aindex200],PfIS[Aindex220],PfIS[Aindex311]]))
            elif PhaseType=="F" or "M":
                PfIS.append(np.mean([PfIS[Mindex110],PfIS[Mindex200],PfIS[Mindex211],PfIS[Mindex220]]))
                #used to include 4, exclude 5, but that is incorrect
            else:
                print ("Unrecognized Phase")
            #####################

            # 5A4F Pairs:  Austenite (111), (200), (220), (311), (222); Ferrite (110), (200), (211), (220)
            # Jacques 2009 Round Robin (XRD6)
            if PhaseType=="A":
                PfIS.append(np.mean([PfIS[Aindex111],PfIS[Aindex200],PfIS[Aindex220],PfIS[Aindex311],PfIS[Aindex222]]))
            elif PhaseType=="F" or "M":
                PfIS.append(np.mean([PfIS[Mindex110],PfIS[Mindex200],PfIS[Mindex211],PfIS[Mindex220]]))
                #used to include 4, exclude 5, but that is incorrect
            else:
                print ("Unrecognized Phase")

            ####################

            # Max unique
            # used in JAC paper
            if PhaseType=="A":
                PfIS.append(np.mean([PfIS[Aindex111],PfIS[Aindex200],PfIS[Aindex220],PfIS[Aindex311],PfIS[Aindex331],PfIS[Aindex420],PfIS[Aindex422],PfIS[Aindex511]]))
            elif PhaseType=="F" or "M":
                PfIS.append(np.mean([PfIS[Mindex110],PfIS[Mindex200],PfIS[Mindex211],PfIS[Mindex310],PfIS[Mindex222],PfIS[Mindex321]]))
                #used to include 4, exclude 5, but that is incorrect
            else:
                print ("Unrecognized Phase")

            OutputList.append(PfIS)

            #print q,SchemeName
            #print "List of average pole Figure Intensities:\n", PfIS
            #print "Average of all pole figures listed: ", sum(PfIS)/len(PfIS)
            #print ""

            #np.mean()
        #print (XPCfile.rsplit('/',1))
        
        if pftype=="xpc":
        #check robustness on other systems
            if os.name=='posix':
                directory, outfile=XPCfile.rsplit('/', 1)  #use this version for MAC
            else:
                directory, outfile=XPCfile.rsplit('\\', 1)  #use this version for Windows
            print (outfile)
        elif pftype=="mtex":
            (head,tail)=os.path.split(MtexFolder)
            if os.name=='posix':
                directory, outfile=MtexFolder.rsplit('/', 1)  #use this version for MAC
            else:
                directory, outfile=MtexFolder.rsplit('\\', 1)  #use this version for Windows
            print (outfile)
        # Output the XPC file name, not sure why it's needed?
        #OutputList.append([XPCfile])
        IntensitiesDF=pd.DataFrame(OutputList)
        IntensitiesDF

        # save to excel

        s=""

        writer = pd.ExcelWriter(s.join([outfile.split('.')[0], ".xlsx"]))
        IntensitiesDF.to_excel(writer,outfile)
        writer.save()
    os.chdir("..")




###################################
# ????
###################################
#Plot pole figures of sampling positions:
#
#
#
#Test functions:
#
#    SchemeName,Coordinates=SingleOrientation("Morris", 60.0,90.0)
#    Coordinates
#
#
#
#Simple plot to work out mplsteronet conventions:
#    """This function is to help us undersrand the conventions of mplstereonet. It accepts no parameters and outputs a plot that translates the 3 Dimensional representations
#    of the sampling sphere into a stereographic 2-Dimensional projection, with the help of creating a sample orientation from the "SingleOrientation" method.
#
#    Parameters
#    ----------
#    None
#        """
#
#
#
#fig = plt.figure(figsize=(8,8), dpi=600)
#
#ax2 = fig.add_subplot(111, projection='stereonet')
#
##SingleOrientation - Name, Tilt, Rotation
##SchemeName,Coordinates=SingleOrientation("RD Single", 90.0,180.0)
#SchemeName,Coordinates=SingleOrientation("TD Single", 90.0,270.0)
#
#dip, strike =Coordinates['Tilt'], Coordinates['Rotation']-90.0
#l1=ax2.pole(strike, dip, 'bD', markersize=10, clip_on=False)
##dip tilts about the 0 axis (RD), left handed
##strike tilts about the normal axis (ND), left handed
##Looking at the bottom of the sphere, not the top
##this is just convention of mplstereonet, does not affect averaging methods
#
#plt.show()



#####################################
#####################################
# End Operations on Intensities
#####################################
#####################################


#####################################
#####################################
# Begin Sampling Schemes
#####################################
#####################################



#####################################
# convert rotations per minute to radians per second
#####################################
def rpm2radpsec(rpm):
    """
    TO ADD: docstring
    """
    import math
    radpsec=(rpm*math.pi*2.0)/60.0
    return radpsec

#####################################
# Define function to tilt and rotate
#####################################
def TiltRotate(name, time, datapoints, rpm,maxtilt,tiltcpm):
    #time=120  # in seconds
    #datapoints = 5000.0
    """
    A function used to represent the tilting angles and rotational angles
    of a sample for a given amount of datapoints, and outputs arrays of
    rotationposition and tiltposition that correspond to certain "time
    values"
    
    Citation: Adapted from C. F. Jatczak, J. A. Larson, and S. W. Shin, Retained austenite and its measurements by X-ray diffraction: an information manual. Warrendale, PA: Society of Automotive Engineers, 1980.
    
    Parameters
    ----------
    name : str
       The name of the sample used
    time : int
       The value representing, in seconds, the total amount of time the sample will
       be subjected to measurements
    datapoints : int
       The number of times we would like to measure the tilt angles and rotation positions of the sample
    rpm : float
       The revolutions per minute expected that our apparatus will revolve; to be converted
       to a more appropriate radians/second measure
    maxtilt : float
       The maximum value the cycle will be tilted(?)
    tiltcpm : float
       The tilt cycles per minute, operating almost as a frequency of the number of times
       the apparatus reaches the other side of the diameter of the pole figure and back to
       original position
    """
    
    import numpy as np
    import pandas as pd
    import math
    from scipy import signal
    rotationspeed=(rpm2radpsec(rpm)) # radians per second
    #maxtilt=60.0  # in degrees
    tiltspeed=(tiltcpm/60.0)  #cycles per second

    #print time,  datapoints, rpm,maxtilt,tiltcpm
    
    timelist=np.ndarray.tolist(np.arange(0, time, time/datapoints))
    #print timelist

    rotationposition=[]
    tiltposition=[]


    #print rotationspeed

    #create the list of tilts and rotations
    for i,item in enumerate(timelist):
        #rotationposition.append(math.degrees(item*rotationspeed % 2.0*math.pi))
        rotationposition.append((math.degrees(rotationspeed*item) % 360.0))

        tiltposition.append( maxtilt* signal.sawtooth(2 * np.pi * tiltspeed * item +np.pi/2, 0.5))

        #print "time:", item, "\tRotation: ", rotationposition[-1],  "\tTilt: ", tiltposition[i]


    #print rotationposition
    # function for tilt
    #56 cycles/minute


    #plt.plot(timelist, rotationposition, 'r',timelist, tiltposition, 'b')
    #plt.show()

    d = {'Tilt' : tiltposition, 'Rotation' : rotationposition}
    coordinates=pd.DataFrame(d)
    #help(TiltRotate)

        #print coordsDF.sort_values('Tilt')    
    return name, coordinates

#############################################
# Define function to create a hexagonal grid:
#############################################
def HexGrid(name, chi_max, angular_spacing, CoverageType="full", IncludeND=True):
    """
    A function used to create a hexagonal grid scheme, to be used for displaying the crystallographic texture of, in this case, various
    steels. The Hexagonal Grid Scheme has experimentally proven to be effective in the reduction of error associated with measuring Austenite
    Phase Fraction, compared to other potential schemes.   This function accepts the sample name, maximum "chi-angle" for measurements
    (elaborated upon in the "Parameters" description), and the desired angle amount to increment each data point by. Outputs arrays of
    "Tilt" and "Rotation" angle values.
    
    Citation: Adapted from A. C. Rizzie, “Elaboration on the Hexagonal Grid and Spiral Method for Data Collection Via Pole Figures,” Spring 2008 [Online]. Available: http://www.bsu.edu/libraries/beneficencepress/mathexchange/05-01/rizzie.pdf. [Accessed: 13-Dec-2016]
    
    Parameters
    ----------
    name : str
       The name of the sample used
       
    chi_max : float
       A float value that represents the maximum possible "chi" degree value (In a pole figure, this represents the maximum possible tilt angle
       made with the normal of the plane and the location of the point on the plane translated upwards onto the surface of the reference sphere
       (an imaginary point).

    angular_spacing : float
       A float value that represents the desired incrementing angle to be applied in measurements; correspondingly increases the
       spacing of the grid.
    
    CoverageType : "full" or "quad"
        generate either "full" pole figure coverage, or only for a "quad" or quadrant of the pole figure.
        Default is "full"
    """
    import numpy as np
    import pandas as pd
    import math
    #chi_max=90.0  #maximum tilt angle in degrees
    #angular_spacing=7.0
    d_max=2.0*math.sin(math.radians(chi_max)/2.0)

    N=d_max/math.radians(angular_spacing)

    #print "Max Tilt: ", chi_max
    #print "Angular Spacing: ", angular_spacing

    xaxis=[] # tilt
    yaxis=[] # rotation

    # the following variable names are taken from Rizzle
    # Note the quadrants aren't quite the same in the mplstereonet covention
    j=0
    i=0
    y_j=0.0
    x_ij=0.0
    chi_ij=0.0 #tilt
    phi_ij=0.0 #rotation

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
            
                # points parallel to ND and -RD-ND plane
                if y_j==0:
                    #do once to avoid duplicates
                    #print "do some nothing"
                    

                    phi_ij=((180.0/math.pi)*math.atan2(y_j,x_ij))
                    #print(phi_ij,chi_ij)
                    
                    if CoverageType=="full":
                        if chi_ij==0:
                            if IncludeND==True:
                                # Just append all points
                                xaxis.append(math.degrees(chi_ij))
                                yaxis.append(phi_ij)
                                #print("Full: ND added")
                            # ND point
                            else:
                                # don't plot the point at chi_ij=0
                                #print("Full: ND skipped")
                                pass
                        else:
                            xaxis.append(math.degrees(chi_ij))
                            yaxis.append(phi_ij)
                            #print("Full: other RD-ND points added")

                    elif CoverageType=="quad":
                        # Is this the right place for this?
                        if chi_ij==0:
                            if IncludeND==True:
                                xaxis.append(math.degrees(chi_ij))
                                yaxis.append(phi_ij)
                                #print("Quad: ND added")
                            else:
                                #print("Quad: ND skipped")
                                pass
                        
                        else:
                            if phi_ij>=-0.2 and phi_ij<=0.2:
                                #print("Quad: -RD skipped")
                                pass
                            # ND point
                            else:
                                xaxis.append(math.degrees(chi_ij))
                                yaxis.append(phi_ij)
                                #print("Quad: +RD added")

                    else:
                        print(CoverageType, "is not a supported option")
                        
                    # Points along the RD-ND plane?
                    if x_ij!=0:
                        phi_ij=((180.0/math.pi)*math.atan2(y_j,-1.0*x_ij))
                        #print(phi_ij,chi_ij)
                        if CoverageType=="full":
                            xaxis.append(math.degrees(chi_ij))
                            yaxis.append(phi_ij)
                            #print("Full: +/-RD added")
                        # FIX somethings' not right here
                        elif CoverageType=="quad":
                            if phi_ij>=-0.2 and phi_ij<=0.2:
                                #print("Quad: -RD skipped RD-ND section")
                                pass
                            else:
                                xaxis.append(math.degrees(chi_ij))
                                yaxis.append(phi_ij)
                                #print("Quad: +/-RD added RD-ND section")
                        else:
                            print(CoverageType, "is not a supported option")
                    else:
                        #removes reduntant rotation
                        pass
                        
                # Points not along y_j=0
                else:
                    if x_ij>0:
                        #quadrant III
                        phi_ij=((180.0/math.pi)*math.atan2(y_j,x_ij)+180.0)
                        xaxis.append(math.degrees(chi_ij))
                        yaxis.append(phi_ij)

                        if CoverageType=="full":
                            #negative i values -quadrant II
                            phi_ij=((180.0/math.pi)*math.atan2(y_j,(-1.0*x_ij)))
                            xaxis.append(math.degrees(chi_ij))
                            yaxis.append(phi_ij)
                       
                            #positive i values - quadrant I
                            phi_ij=((180.0/math.pi)*math.atan2(y_j,x_ij))
                            xaxis.append(math.degrees(chi_ij))
                            yaxis.append(phi_ij)


                            #quadrant IV
                            phi_ij=((180.0/math.pi)*math.atan2(y_j,(-1.0*x_ij))+180.0)
                            xaxis.append(math.degrees(chi_ij))
                            yaxis.append(phi_ij)
                            
                        elif CoverageType=="quad":
                            pass
                        
                        else:
                            print(CoverageType, "is not a supported option")

                    # Points along the TD-ND plane
                    else:
                        if CoverageType=="full":
                            # poisnt along -TD-ND plane
                            phi_ij=(90.0)
                            xaxis.append(math.degrees(chi_ij))
                            yaxis.append(phi_ij)
                            
                            # points along +TD-ND plane
                            phi_ij=(270.0)
                            xaxis.append(math.degrees(chi_ij))
                            yaxis.append(phi_ij)
                        
                        elif CoverageType=="quad":
                            # points along along +TD-ND plane
                            phi_ij=(270.0)
                            xaxis.append(math.degrees(chi_ij))
                            yaxis.append(phi_ij)
                        
                        else:
                            print(CoverageType, "is not a supported option")
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
    coordinates=pd.DataFrame(d)

    #print coordsDF.sort_values('Tilt')    
    return name, coordinates


###################################

###################################
#Define function to create a spiral grid scheme (INCOMPLETE):
###################################
def SpiralScheme(name, chi_max, phi_max):
    """
    A function used to create a spiral grid scheme (a polar graphical representation of data in certain atomic planes, to be used in analysis of material texture).
    This function accepts three parameters: the name of the sample to be used, the maximum chi degree value (elaborated upon below in the parameters description), as well
    as the maximum phi degree value (similarly focused upon in the parameters section of this docstring. Like other methods in this family of grid schemes, it outputs arrays
    of tilt angle (represented by chi) and rotation angle (represented by phi).

    This function relies on a logarithmic spiral pattern, which expands outwards at a quicker rate than an Archimedean Spiral. 

    Citation: Adapted from A. C. Rizzie, “Elaboration on the Hexagonal Grid and Spiral Method for Data Collection Via Pole Figures,” Spring 2008 [Online]. Available: http://www.bsu.edu/libraries/beneficencepress/mathexchange/05-01/rizzie.pdf. [Accessed: 24-Feb-2021]

    Parameters
    ----------
    name : str
       The name of the sample used
       
    chi_max : float
       A float value that represents the maximum possible "chi" degree value (In a pole figure, this represents the maximum possible tilt angle
       made with the normal of the plane and the location of the point on the plane translated upwards onto the surface of the reference sphere
       (an imaginary point).

    phi_max : float
       A float value that represents the maximum possible "phi" degree value (In a pole figure and using our spiral grid scheme, this is related to the amount of times our cycle
       will revolve around the pole figure as it moves radially outwards and collects points throughout the plane.
    """
    xaxis=[] #Chi Array of values to be outputted
    yaxis=[] #Phi Array of values to be outputted
    i=1 # Chi Counter, important for the determination of Phi 
    chi_N=0 
    phi_N=0
    
    xaxis.append(0) #The paper sees we must manually input the coordinate values (0,0) before we enter the
    
    yaxis.append(0)# while-loop.
    
    xaxis.append(5) #Since this is also a one-time occurence, I decided to take care of it beforehand...
    
    yaxis.append(0)
    
    chi_N=5
    
    b=np.log(chi_max)/phi_max # A constant that is to be used in the values for Chi_n
    while (phi_N<=phi_max): # might change this condition to while (phi...)
        phi_N=phi_N+ (5*chi_max/xaxis[i]) 
        yaxis.append(phi_N)
        i=i+1
        if (phi_N<phi_max):
            chi_N=np.exp(b*phi_N)
        elif (phi_N==phi_max):
            chi_N=chi_max
            xaxis.append(chi_N)
            break
        
        xaxis.append(chi_N)
               
            # Need to perform the translations (arctan, scale down by 1/45 as the research paper, but would like to
            #understand why this needs to be done before implementing it...)
        
    
    
    del(xaxis[len(xaxis)-1])
    del(yaxis[len(yaxis)-1])
    xaxis.append(chi_max)
    yaxis.append(phi_max)
    values = {'Tilt' : xaxis, 'Rotation' : yaxis}
    coordinates=pd.DataFrame(values)   
    DroppedVals = coordinates[(coordinates['Tilt'] > 0.0) & (coordinates['Tilt'] < 5.0)].index #drops the chi values less than 5.0 as prescribed by the document
    coordinates.drop(DroppedVals , inplace= True)                      
    return name, coordinates



###################################
# Rings Perpendicular to ND
###################################
def RingPerpND(res):
    """
    This method takes a look at graphically exploring specific rings perpendicular to the normal direction of the sample. The parameters to be accepted are very simple:
    only the desired angluar distance to be incremented as we measure the rotation degree value along the circular plane. Outputs the directional relationship between the ring and 
    the sample direction, and a set of coordinates, with "Tilt" being 90 degrees for each of the corresponding "Rotation values" as we work around the circular shaped plane. 
    Tilt is always 90 degrees, of course, because the ring is perpendicular to the normal direction.
    
   
    Parameters
    ----------
    res : float
        A float value that determines the step size of the rotation array base as it goes from 0 to 360 degrees (around the circular plane).
    """
    import numpy as np
    import pandas as pd
    import math
    
    # Perpendicular to ND

    name="Ring Perpendicular to ND"
    yaxis=np.ndarray.tolist(np.arange(0.0, 360.0001, res))#rotation
    xaxis=[90.0] * len(yaxis) #tilt 
    d = {'Tilt' : xaxis, 'Rotation' : yaxis}
    coordsDF=pd.DataFrame(d)
    return name, coordsDF
   
###################################
# Rings Perpendicular to RD
###################################
def RingPerpRD(res):
    """This method takes a look at graphically exploring specific rings perpendicular to the rolling direction of the sample. The parameters to be accepted are very simple:
    only the desired angluar distance to be incremented as we measure the tilt degree value along the circular plane. Outputs the directional relationship between the ring and 
    the sample direction, and a set of coordinates, with "Rotation" values being a constant 90 or 270 degrees, and tilt varying until we hit 90 degrees (no need for the full 180 
    degrees because the two sections are symmetric and identical to each other).
   
    Parameters
    ----------
    res : float
        A float value that determines the step size of the tilt array base as it goes from 0 to 90 degrees (around the circular plane).
        """
    import numpy as np
    import pandas as pd
    import math
    
   # Perpendicular to RD

    name="Ring Perpendicular to RD"
    xaxis=np.ndarray.tolist(np.arange(0.0, 90.001, res))#tilt
    #yaxis=[90.0] * len(xaxis) #rotation
    yaxis=[270.0] * len(xaxis) #rotation
    d = {'Tilt' : xaxis, 'Rotation' : yaxis}
    coordsDF=pd.DataFrame(d)
    return name, coordsDF
   # Perpendicular to RD

    
###################################
# Rings Perpendicular to RD
###################################
def RingPerpTD(res):
    """This method takes a look at graphically exploring specific rings perpendicular to the transverse direction of the sample. The parameters to be accepted are very simple:
    only the desired angluar distance to be incremented as we measure the tilt degree value along the circular plane. Outputs the directional relationship between the ring and 
    the sample direction, and a set of coordinates, with "Rotation" values being a constant 0 or 180 degrees, and tilt varying until we hit 90 degrees (no need for the full 180 
    degrees because the two sections are symmetric and identical to each other).
    
   
    Parameters
    ----------
    res : float
        A float value that determines the step size of the tilt array base as it goes from 0 to 90 degrees (around the circular plane).
        """
    
    import numpy as np
    import pandas as pd
    import math
    
    # Perpendicular to TD

    name="Ring Perpendicular to TD"
    xaxis=np.ndarray.tolist(np.arange(0.0, 90.001, res))#tilt Note for me (Surya) Ask Dr. Creuziger why this doesn't go upto 180.001
    #yaxis=[0.0] * len(xaxis) #rotation
    yaxis=[180.0] * len(xaxis) #rotation
    d = {'Tilt' : xaxis, 'Rotation' : yaxis}
    coordsDF=pd.DataFrame(d)
    return name, coordsDF


###################################
# Single Orientation
###################################
def SingleOrientation(name, tilt, rotation):
    """
    This function defines a singular, specific orientation for sampling, as it works in concert with other methods that output tilt and rotation arrays. It accepts the
    sample name, "tilt" array and "rotation" array that have been provided by other methods, and simply outputs a neater version of the data, giving it a name and creating
    a Pandas dataframe that records the tilt and rotation of each relevant data point.
    
   
    Parameters
    ----------
    name : str
        The name of the sample orientation to be used
     
    tilt : array of float values
        An array of float values representing the corresponding tilt angles (or chi), to be used as reference for pole figure geometry
    
    rotation : array of floar values
        An array of floar values representing the corresponding rotation angles (or phi), to be used as reference for pole figure geometry.
        """
    import numpy as np
    import pandas as pd
    import math
    
    coordslist=[[tilt,rotation]]
    xaxis=[]
    yaxis=[]
    for item in coordslist: 
        xaxis.append(item[0])
        yaxis.append(item[1])
    #coordslist[0][:]
    d = {'Tilt' : xaxis, 'Rotation' : yaxis}
    coordsDF=pd.DataFrame(d)
    return name, coordsDF

#######################################
#Equal Angle Sampling Scheme
#######################################

# Perpendicular to RD
def RingRot(res,theta,omega_list,rotaxis, Weight=False): #add omega
    """

    res: Resolution about the rings
    theta: Theta value (1/2 2-theta), deviation from vertical due to finite energy
    omega_list: List of the angular rotations about the axis
    rotaxis: Which axis to rotate the rings about

    FIX - Y rotaxis has an extra line in X direction
    """
    import numpy as np
    import pandas as pd
    import math
    
    if Weight==True:
        name="RotRing Axis-%s Res-%s Theta-%s OmegaMax-%s Weighted" % (rotaxis, res,theta,max(omega_list))
    else:
        name="RotRing Axis-%s Res-%s Theta-%s OmegaMax-%s" % (rotaxis, res,theta,max(omega_list))
    
    #name="Ring Perpendicular to ND"
    #print rotaxis
    
    #Generate the ring of rotations
    yaxis=np.ndarray.tolist(np.arange(0.0, 360.0, res))#rotation
    xaxis=[90.0-theta] * len(yaxis) #tilt
    # Don't save, or this gets repeated
    #d = {'Tilt' : xaxis, 'Rotation' : yaxis,'Weights' : np.ones(len(xaxis))}
    #coordsDF=pd.DataFrame(d)
    
    ## Additional rotations
    #omega_list=[5]
    
    for omega_counter, omega_step in enumerate(omega_list):
        #print "Generating omega", omega_counter
        #print "Additonal omegas: ", omega_step
        r_list=[]
        t_list=[]
        weight_list=[]
        rotation_index=[]
        omega_step_list=[]
        for i, value in enumerate(yaxis):
            #print "Creating Point on ring", i
            if rotaxis=='X':
                #print "Rotating in X axis"
                (r,t)=cart2sphDeg(np.einsum('ij,j', RotateXMatrix(omega_step), sph2cartDeg(value,xaxis[i])))
                
                # Weight by:
                # Area of a zone in a sphere (also called spherical segment)
                # Area of a zone = 2*pi*R*h.  h is the x1-x2 distance
                # Area of sphere = 4*pi*r^2
                # Area of zone normalized by unit sphere = 0.5*h (r=1)
                
                #check if close to pole
                if (value!=0.0):
                    x1,y,z=sph2cartDeg(value+res/2.0,xaxis[i])
                    x2,y,z=sph2cartDeg(value-res/2.0,xaxis[i])
                    weight=0.5*abs(x1-x2)
                else:
                    x1,y,z=sph2cartDeg(value,xaxis[i]) # the segment curves around when close to the pole
                    x2,y,z=sph2cartDeg(value-res/2.0,xaxis[i])
                    weight=0.5*abs(x1-x2)

            elif rotaxis=='Y':
                #print "Rotating in Y axis"
                (r,t)=cart2sphDeg(np.einsum('ij,j', RotateYMatrix(omega_step), sph2cartDeg(value,xaxis[i])))
                
                # Weight by:
                # Area of a zone in a sphere
                if (value!=90.0):
                    x,y1,z=sph2cartDeg(value+res/2.0,xaxis[i])
                    x,y2,z=sph2cartDeg(value-res/2.0,xaxis[i])
                    weight=0.5*abs(y1-y2)
                else:
                    x,y1,z=sph2cartDeg(value,xaxis[i]) # the segment curves around
                    x,y2,z=sph2cartDeg(value-res/2.0,xaxis[i])
                    weight=0.5*abs(y1-y2)


                
                #using this as an index on the rotation
                # also fix omega step below
                
                
            else:
                print("Not a supported rotation axis")
            
            #(r,t)=StepRotateRing(value,xaxis[i],omega_step)
            #print value, xaxis[i], omega_step,r,t
            
            #print "Adjusting angles"
            # Add some logic to keep 0<=r<360 and 0<=t<90
            if r<0:
                r=r+360
            elif r>=360:
                r=r-360

            # correct for tilt below the sphere.  Don't need to correct rotation here
            if t<0:
                t=-t
            if t>90:
                t=180-t

            
            #print "Appending lists"
            r_list.append(r)
            t_list.append(t)
            rotation_index.append(i)
            omega_step_list.append(omega_step)
            #print weight, np.sin(np.radians(t))
            
            weight_list.append(weight)
            
            # correct for duplicate points at ±90°
#             if abs(omega_step)==90.0:
#                 #weight_list.append(weight/2.0)
#                 weight_list.append(weight) #only for indexing
#             else:
#                 weight_list.append(weight)
#                 #weight_list.append(weight)
            #22 Oct 2019 need to change weight here, after correcting t - doesn't work
            #weight_list.append(weight*np.sin(np.radians(t)))
        #print "Creating Data Frames"
        #print "Counter value",omega_counter
        
        # normalize weights by number of data points and weight values

        #print(weight_list)
        # still inside a single ring, so just need number of points
        weight_nparray=np.array(weight_list)
        normfactor=sum(weight_nparray)*len(omega_list) #includes the number of points
        weight_nparray=weight_nparray/normfactor
        
        # Flag weight to reduce weigtht for +/- 90 omegas due to overlap in the +ND/-ND plane
        if omega_step==90:
            weight_nparray=weight_nparray/2.0
            #print("omega 90")
        elif omega_step==-90:
            weight_nparray=weight_nparray/2.0
            #print("omega -90")
        else:
            pass
        weight_list=list(weight_nparray)
        #print(weight_list)

        if omega_counter==0:
            # initialize dataframe
            #print "Intialize Data Frame at x=",omega_counter
            d2 = {'Tilt' : t_list, 'Rotation' : r_list,
                  'RotationIndex' : rotation_index,'Omega' : omega_step_list,'Weights' : weight_list}
            coordsDF=pd.DataFrame(d2)
            
        else:
            d2 = {'Tilt' : t_list, 'Rotation' : r_list,
                  'RotationIndex' : rotation_index,'Omega' : omega_step_list, 'Weights' : weight_list}
            d2DF=pd.DataFrame(d2)
            coordsDF=coordsDF.append(d2DF, ignore_index=True)
    
    #normalization by number of points, needs to be all points since the ±90 rings are weighted differently
    coordsDF['Weights']=coordsDF['Weights']/sum(coordsDF['Weights'])
        
        
        
    # range and step of rotations
    
    # append
    #df2 = pd.DataFrame([[5, 6], [7, 8]], columns=list('AB'))
    #df.append(df2)
    
    return name, coordsDF


# In[12]:


# adapted from Steronet
### Helper function for RotRing?
def sph2cartDeg(rotation_deg, tilt_deg):
    """
    Converts a longitude and latitude (or sequence of lons and lats) given in
    _radians_ to cartesian coordinates, `x`, `y`, `z`, where x=0, y=0, z=0 is
    """
    import numpy as np
    import math
    elevation_deg=90-tilt_deg
    elevation=elevation_deg*math.pi/180
    rotation=rotation_deg*math.pi/180
    
    x = np.cos(elevation) * np.cos(rotation)
    y = np.cos(elevation) * np.sin(rotation)
    z = np.sin(elevation)

    return [x, y, z]
    


### Helper function for RotRing?
def cart2sphDeg(a):
    """

    """
    import numpy as np
    import math
    x=a[0]
    y=a[1]
    z=a[2]
    #https://www.mathworks.com/help/matlab/ref/cart2sph.html
    azimuth = np.arctan2(y,x)
    elevation = np.arctan2(z,np.sqrt(x**2 + y**2))
    #tilt=
    
    
    #check on value of r?
    r = np.sqrt(x**2 + y**2 + z**2)
    #lat = np.arcsin(z/r) # original
    #lat = np.arccos(z/r) # changed for tilt/rotate convention - didn't work out
    #lat = -np.arccos(z/r) # changed to invert with RotateTilt

    #lon = np.arctan2(y, x)
    return np.degrees(azimuth), 90-np.degrees(elevation)


### Helper function for RotRing?
def RotateYMatrix(omega_deg):
    """
    Helper function for RotRing, hard coded rotations about the Y axis
    """

    import numpy as np
    import math
    omega=omega_deg*math.pi/180
    R=[[math.cos(omega),0, -math.sin(omega)  ],
        [0, 1,0],
        [math.sin(omega), 0, math.cos(omega)]]
    return np.transpose(np.array(R)) # transpose to change from active to passive



### Helper function for RotRing?
def RotateXMatrix(omega_deg):
    """
    Helper function for RotRing, hard coded rotations about the Y axis
    """
    import numpy as np
    import math
    omega=np.radians(omega_deg)
    R=[[1,0, 0  ],
        [0, math.cos(omega),math.sin(omega)],
        [0,-math.sin(omega), math.cos(omega)]]
    return np.transpose(np.array(R)) # transpose to change from active to passive




#######################################
#Equal Angle Sampling Scheme
#######################################

def EqualAngleGrid(name,chi_max, stepsize, CoverageType ='full'):
    """
    The Equal Angle Grid is currently how most XRD systems tend to operate, evenly stepping 5 degrees in rotation until one revolution is
    complete, then incrementing tilt by that same 5 degrees and repeating till maximum tilt is obtained. 
    While easy to understand, this scheme tends to oversample around the ND axis, because 360/5= 72 sampling points are taken for each tilt "ring", 
    regardless of tilt value. The Equal Angle scheme also requires nearly 1400 sampling points, whereas novel hex schemes can obtain even sampling of the
    pole grid while requiring fewer sampling points. Thus, the Equal Angle scheme can be compared with traditional hex schemes to demonstrate that fewer
    sampling without compromising measurement accuracy can be obtained.

    Parameters
    ----------
    name : string
        The name of the sampling scheme.
     
    chi_max : float
        Maximum tilt of the pole figure, in degrees. Typically, 60 degrees is the practical limit due to defocusing phenomena.
    
    stepsize : float
        Value to evenly increment both tilt and rotation steps by, in degrees. 5 degrees is the recommended parameter to account for smoothly-varying textured materials.
    
    CoverageType : string
        Determines full (the default option) or partial coverage of the pole figure. Entries of 'full' or 'quad' designate 'full' coverage or 'single-quadrant' coverage, respectively. 



    

    """
    tilt=[]
    rotation=[]
    i=0
    j=0
    tilt.append(i)
    rotation.append(j)
    if CoverageType.lower()=='full':
        phi_max=360.0

    elif CoverageType.lower()=='quad':
        phi_max=90.0

    else:
        raise Exception("Invalid Coverage Type: Select between \"full\" and \"quad\".")

    while(round(i)<=chi_max):
        while(round(j)<=phi_max):
            rotation.append(j)
            j=j+stepsize
            tilt.append(i)
        j=0
        i=i+stepsize
    d = {'Tilt' : tilt, 'Rotation' : rotation}
    coordsDF=pd.DataFrame(d)
    return name, coordsDF
        
def CLRGrid(name,chi_max, stepsize, CoverageType='full'):
    """
    In response to the oversampling observed along the ND axis in the equal-angle scheme, this CLR scheme attempts a more even sampling of the pole grid
    by implementing a "constant local resolution," which is what the "CLR" stands for. For a constant increasing tilt, the rotation increment decreases asymptotically,
    such that the rotation increment equals the tilt increment when tilt is maximum. 
    Due to some slight issues with the formula used, some locations of the CLR scheme are imprecise, but some approximations have been calculated to mostly 
    allow for a smoothly-sampled pole grid.
     

    Parameters
    ----------
    name : string
        The name of the sampling scheme.
     
    chi_max : float
        Maximum tilt of the pole figure, in degrees. Typically, 60 degrees is the practical limit due to defocusing phenomena.
    
    stepsize : float
        Value to evenly increment both tilt and rotation steps by, in degrees. 5 degrees is the recommended parameter to account for smoothly-varying textured materials.
    
    CoverageType : string
        Determines full (the default option) or partial coverage of the pole figure. Entries of 'full' or 'quad' designate 'full' coverage or 'single-quadrant' coverage, respectively. 


    """
    tilt=[]
    rotation=[]
    error=set() #debug purposes
    step=math.radians(stepsize)
    j=0
    k=1.0 
    tilt.append(0)
    rotation.append(0)

    if CoverageType.lower()=='full':
        phi_max=360.0

    elif CoverageType.lower()=='quad':
        phi_max=90.0

    else:
        raise Exception("Invalid Coverage Type: Select between \"full\" and \"quad\".")
    

    while(k<=chi_max/stepsize):
        while(j<phi_max):
            rotation.append(j)
            phi_step=math.degrees(math.acos(((math.cos(step)-math.pow(math.cos(k*step),2)))/(math.pow(math.sin(k*step),2))))
            if (360%phi_step!=0):
                numsteps=round(360.0/phi_step)
                phi_step=360.0/numsteps
                #tup=phi_step,k,numsteps
                #error.add(tup)
            j=j+phi_step
            tilt.append(k*stepsize)
        k=k+1 
        j=0


    #m=sorted(error)
    #m.reverse()    
    #return pd.DataFrame(data=m,columns=['Rotation Step', 'K-value (ring number)', 'Number of points per rotation']) 

    d = {'Tilt' : tilt, 'Rotation' : rotation}
    coordsDF=pd.DataFrame(d)
    return name, coordsDF

def BT8_HexGrid(name, chi_max, stepsize, CoverageType="full"):
    """
    Hex schemes are considered one to be one of the most promising sampling schemes, given their combination of fewer sampling points and ability to
    mitigate phase fraction measurement bias due to crystallographic texture. This variation of the hex scheme samples half as many points as a traditional
    hex scheme, so it is being tested for its ability to reduce texture effects on phase fraction calculations.
     

    Parameters
    ----------
    name : string
        The name of the sampling scheme.
     
    chi_max : float
        Maximum tilt of the pole figure, in degrees. Typically, 60 degrees is the practical limit due to defocusing phenomena.
    
    stepsize : float
        Value to evenly increment both tilt and rotation steps by, in degrees. 5 degrees is the recommended parameter to account for smoothly-varying textured materials.
    
    CoverageType : string
        Determines full (the default option) or partial coverage of the pole figure. Entries of 'full' or 'quad' designate 'full' coverage or 'single-quadrant' coverage, respectively. 


    """
    i=1
    j=0
    k=0
    rotation=[]
    tilt=[]
    nphi=round(360/stepsize)
    n=nphi
    t=0.0
    if CoverageType.lower()=="full":

        while(6*i<=nphi):
            n=n+nphi-i*6
            i=i+1
        
        for i in np.arange(0,round(nphi/6)+1,step=1):
            for j in range(1,(nphi-6*i)+1):
                chi=2*(np.arcsin((math.sqrt(2)/2/(nphi/6)*((nphi/6) - i))))*180/np.pi
                if round(chi)<=chi_max:
                    rotation.append((360/(nphi-6*i))*(j-1))
                    if chi>90.0:
                        tilt.append(180-chi)
                    else:
                        tilt.append(chi)

    

    elif CoverageType.lower()=="quad":

        while(6*i<=nphi):
            n=n+nphi-i*6
            i=i+1
        
        for i in np.arange(0,round(nphi/6)+1,step=1):
            for j in range(1,(nphi-6*i)+1):
                phi=(360/(nphi-6*i))*(j-1)
                chi=2*(np.arcsin((math.sqrt(2)/2/(nphi/6)*((nphi/6) - i))))*180/np.pi
                if round(chi)<=chi_max:
                    if phi<=90.0:
                        rotation.append(phi)
                        if chi>90.0:
                            tilt.append(180-chi)
                        else:
                            tilt.append(chi)

     
    
    
    
    rotation.append(0)
    tilt.append(0) 
    
    d = {'Tilt' : tilt, 'Rotation' : rotation}
    coordsDF=pd.DataFrame(d)
    return name, coordsDF


def SpiralGrid(name, resolution):
    """
    The initial inspirations behind the spiral
    grid were to continuously and smoothly iterate through the
    range of tilt and rotation angles. In order to have an even sampling from this spiral grid, 
    the number of sampling points at certain tilt rings was set directly proportional to the magnitude of the
    sine of the tilt angle at a given sampling ring. Then, the discrete nature of the tilt rings was removed and the tilt angle
    was also dynamically changed, thereby reducing the total number of sampling points used in the sampling scheme. The final results can be 
    better appreciated with the plotting options seen in the code (via matplotlib), plotting the number of sampling points in a tilt regions 
    as a function of instantaneous coordinate angle value.


    Parameters
    ----------
    name : string
        The name of the sampling scheme.
     
    resolution: float
        The spacing between sampling points in the pole figure grid. A higher resolution value implies more sampling points used, or a denser 
        distribution of sampling points, and vice versa for lower resolution values. 5 degrees is assumed to be the conventional resultion used 
        in this project.



    """
    import matplotlib.pyplot as plt
    tilt=[]
    rotation=[]
    rings=[]
    chi=90.0
    phi=0.0
    count=0
    while(chi>0):
        while(phi<360):
            points=(360/resolution)*math.sin(math.radians(chi))
            phi_step=float(360/points)
            chi_step=float(resolution/points)
            rotation.append(phi)
            tilt.append(chi)
            phi=phi+phi_step
            chi=chi-chi_step
            count=count+1    
        phi=phi-360
        rings.append(count)
        count=0
        
    

    tilt.append(0)
    rotation.append(0)

    #Plot functions for tilt/rotation over time
    
    width = 0.35       # the width of the bars: can also be len(x) sequence

    labels=[]
    i=0
    for i in range(len(rings)):
        labels.append(str(i+1))
    
    
    plt.bar(labels,rings)


    
    plt.xlabel('Tilt Range (Starting from higher Tilt)')
    plt.ylabel('Number of Sampling Points')
    plt.title('Sampling points of Spiral Grid as a Function of Tilt')

    plt.show()
    
    






    d = {'Tilt' : tilt, 'Rotation' : rotation}
    coordsDF=pd.DataFrame(d)
    return name, coordsDF
    

def OffsetRing(name, res, theta, tilt_center , rotation_center):
    """
    The main motivation behind the Offset ring is to account for 1 diffraction peak in each of the ferrite and austenite phases. This would
    be done by 2 seperate sampling rings, which would be shifted by a certain value to account for the fact that X-ray energies are finite 
    as opposed to infinite (along the ND axis).
    
    This Offset single ring implementation was largely adapted from the 'Ringrot' function and modified appropriately for these purposes, but there 
    are still some improvements that need to be made in order for this scheme to be feasible for phase fraction measurement of textured 
    steels.

    Parameters
    ----------
    name : string
        The name of the sampling scheme.
        
    res:
        Resolution about the rings
    theta:
        Theta value (1/2 2-theta), deviation from vertical due to finite energy
    tilt_center : float (degrees)
        Tilt angle for the center of the ring
    rotation_center: float (degrees)
        rotation angle for the center of the ring

    """

    import numpy as np
    import pandas as pd
    import math
    
    yaxis=np.ndarray.tolist(np.arange(0.0, 360.0, res))
    xaxis=[90.0-theta] * len(yaxis) 
    
    r_list=[]
    t_list=[]
    for i, value in enumerate(yaxis):
        (r,t)=cart2sphDeg(np.einsum('ij,j', RotateTiltRotation(tilt_center , rotation_center), sph2cartDeg(value,xaxis[i])))
            
        
        if r<0:
            r=r+360
        elif r>=360:
            r=r-360

        if t<0:
            t=-t
        if t>90:
            t=180-t

        
        
        r_list.append(r)
        t_list.append(t)

    df = {'Tilt' : t_list, 'Rotation' : r_list}
    coordsDF=pd.DataFrame(df)
    
    
    return name, coordsDF 


### Helper function for RotRing?
def RotateTiltRotation(psi_deg, phi_deg):
    """
    Helper function for OffsetRing,
    Uses Rotation matrix from "Residual Stress" by Noyan & Cohen (1987)
    """

    import numpy as np
    import math
    psi=psi_deg*math.pi/180
    phi=phi_deg*math.pi/180
    R=[[math.cos(phi)*math.cos(psi),math.sin(phi)*math.cos(psi), -math.sin(psi)  ],
        [-math.sin(phi), math.cos(phi),0],
        [math.cos(phi)*math.sin(psi), math.sin(phi)*math.sin(psi), math.cos(psi)]]
    return np.transpose(np.array(R)) # transpose to change from active to passive



#print(EqualAngleSampling("Austenite",30)) (Debug)

#####################################
#####################################
# End Sampling Schemes
#####################################
#####################################
