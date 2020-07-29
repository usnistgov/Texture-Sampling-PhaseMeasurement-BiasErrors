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

# ## Packages used in this Jupyter Notebook

# In[1]:


import fortranformat as ff
import numpy as np 
import pandas as pd
import scipy as scipy
from scipy import interpolate 
from scipy import signal
import matplotlib.pyplot as plt
import mplstereonet
import math
import os


# ## Function: Read in Polefigure in .xpc format
# MAUD uses an .xpc format for pole figures, likely derived from BearTex [1].  This format is similar to the General Intensity File Format in POPLA [2, appendix B2], with a slightly different header
# 
# [1] http://eps.berkeley.edu/~wenk/TexturePage/beartex.htm
# 
# [2] Popla Manual http://pajarito.materials.cmu.edu/rollett/27750/popLA_Manual.pdf
# 

# In[ ]:


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


def xpcformat(mode=None, filename=None):
    """
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

    print "Pole Figure Parsing"

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
        print 'You are now reading experimental pole figure(s) :%s'%filename
        blocks = open(filename, 'rU').read().split('\n\n\n\n')[1:]
        print 'There are %s blocks of data found'%len(blocks)
        if len(blocks)==0:
            msg1 = 'xpc parser in upf assumes that pole figures are separated by 4 new lines'
            msg2 = ' searching %s finds no set of 4 new lines in '%filename
            msg  = '%s \n %s'%(msg1,msg2)
            raise IOError, msg
            # blocks = parse_epf(filename)
        npf = len(blocks)
        if npf==0: raise IOError, 'No pf block found.'

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
                parsed=dataline.read(line[item[0]])
                parsed.extend(dataline.read(line[item[1]]))
                parsed.extend(dataline.read(line[item[2]]))
                parsed.extend(dataline.read(line[item[3]]))
                dataset.append(parsed)
            #print dataset
            
            #Saves as a Pandas dataframe, and maps the 360 degree phi data from the 0 degree phi data
            #row and column indexes are by degrees
            df=pd.DataFrame(dataset, index=np.arange(0,91,5))
            df.columns=[np.arange(0,360,5)]
            df[360]=df.ix[:,0]
            
            # Save the hkl value
            hkl = [h,k,l] #hkl
            #print hkl
            hkls.append(hkl)

            datasets.append(df)
            
        #print hkls
        print "number of pole figures:", len(datasets)

        return datasets, hkls
    else: raise IOError, 'Unexpected mode is given'
    #return data


# ## Define function to take a pole figure (or series of pole figures) and series of coordinate pairs, returning an average intensity for all of the pole figures and coordinates

# In[ ]:


def pfIntensitySum(name, PoleFigures, Coordinates):

    """
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
        for index, row in Coordinates.iterrows():
            #print row['Tilt'], row['Rotation'], InterpPF.ev(row['Tilt'],row['Rotation'])
            IntensityValues.append(InterpPF.ev(row['Tilt'],row['Rotation']))

        #print IntensityValues
        #Factor of 100 is divided to convert from POPLA style format of 100 = 1 MRD/MUD
        AverageIntensity.append(sum(IntensityValues)/(100*len(IntensityValues)))

    #print AverageIntensity

    ## return average value
    return AverageIntensity


# ## Define function to create a hexagonal grid
# 
# Adapted from A. C. Rizzie, “Elaboration on the Hexagonal Grid and Spiral Method for Data Collection Via Pole Figures,” Spring 2008 [Online]. Available: http://www.bsu.edu/libraries/beneficencepress/mathexchange/05-01/rizzie.pdf. [Accessed: 13-Dec-2016]
# 

# In[2]:


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
    coordinates=pd.DataFrame(d)

    #print coordsDF.sort_values('Tilt')    
    return name, coordinates


# ## Define function to tilt and rotate
# 
# Adapted from C. F. Jatczak, J. A. Larson, and S. W. Shin, Retained austenite and its measurements by X-ray diffraction: an information manual. Warrendale, PA: Society of Automotive Engineers, 1980. 
# 

# In[3]:


# convert rotations per minute to radians per second
def rpm2radpsec(rpm):
    radpsec=(rpm*math.pi*2.0)/60.0
    return radpsec


# In[4]:


def TiltRotate(name, time, datapoints, rpm,maxtilt,tiltcpm):
    #time=120  # in seconds
    #datapoints = 5000.0
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

        tiltposition.append( maxtilt * signal.sawtooth(2 * np.pi * tiltspeed * item +np.pi/2, 0.5))

        #print "time:", item, "\tRotation: ", rotationposition[-1],  "\tTilt: ", tiltposition[i]


    #print rotationposition
    # function for tilt
    #56 cycles/minute


    #plt.plot(timelist, rotationposition, 'r',timelist, tiltposition, 'b')
    #plt.show()

    d = {'Tilt' : tiltposition, 'Rotation' : rotationposition}
    coordinates=pd.DataFrame(d)

        #print coordsDF.sort_values('Tilt')    
    return name, coordinates


# # Define function for Spiral

# In[5]:


def Spiral(name, maxtilt, numrotation, spacing):
    """Name refers to the specific scheme. maxtilt is the maximum tilt angle measured in degrees.
    numrotation is the number of times the spiral rotates. spacing is the relative spacing between datapoints; a default of 5 is recommended."""
    rotationposition=[0.0,0.0]
    tiltposition=[0.0,5.0]
    
    spacingcoefficient = np.log(maxtilt)/(math.degrees(numrotation * 2 * np.pi))
    
    while rotationposition[len(rotationposition)-1]<(math.degrees(numrotation * 2 * np.pi)) and tiltposition[len(tiltposition)-1]<maxtilt:
        rotationposition.append(rotationposition[len(rotationposition)-1] + spacing*maxtilt/tiltposition[len(tiltposition)-1])
        tiltposition.append(np.exp(spacingcoefficient * rotationposition[len(rotationposition)-1]))

    del(rotationposition[len(rotationposition)-1])
    del(tiltposition[len(tiltposition)-1])
        
    rotationposition.append(math.degrees(numrotation * 2 * np.pi))
    tiltposition.append(maxtilt)
    
    d = {'Tilt' : tiltposition, 'Rotation' : rotationposition}
    coordinates=pd.DataFrame(d,columns=['Tilt','Rotation'])
    coordinates = coordinates.drop(coordinates[(coordinates.Tilt>0.0) & (coordinates.Tilt<5.0)].index)
    return name, coordinates    


# # Define function for Gaussian Quadrature

# In[6]:


def GaussQuad(name):
    rotationposition=[0.0,0.0,45.0,90.0,135.0,0.0,45.0,90.0,135.0,0.0,45.0,90.0,135.0,0.0,45.0,90.0,135.0]
    tiltposition=[0.0,0.533296*180/math.pi,0.533296*180/math.pi,0.533296*180/math.pi,0.533296*180/math.pi,1.2239*180/math.pi,1.2239*180/math.pi,1.2239*180/math.pi,1.2239*180/math.pi,1.91769*180/math.pi,1.91769*180/math.pi,1.91769*180/math.pi,1.91769*180/math.pi,2.6083*180/math.pi,2.6083*180/math.pi,2.6083*180/math.pi,2.6083*180/math.pi]
    d = {'Tilt' : tiltposition, 'Rotation' : rotationposition}
    coordinates = pd.DataFrame(d)
    return name, coordinates


# # Create list of positions to evaluate

# ## Single Sample Orientation

# In[7]:


def SingleOrientation(name, tilt, rotation):
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


# ## Rings Perpendicular to a sample direction

# In[8]:


# Perpendicular to ND
def RingPerpND(res):
    name="Ring Perpendicular to ND"
    yaxis=np.ndarray.tolist(np.arange(0.0, 360.0001, res))#rotation
    xaxis=[90.0] * len(yaxis) #tilt 
    d = {'Tilt' : xaxis, 'Rotation' : yaxis}
    coordsDF=pd.DataFrame(d)
    return name, coordsDF


# In[9]:


# Perpendicular to RD
def RingPerpRD(res):
    name="Ring Perpendicular to RD"
    xaxis=np.ndarray.tolist(np.arange(0.0, 90.001, res))#tilt 
    #yaxis=[90.0] * len(xaxis) #rotation
    yaxis=[270.0] * len(xaxis) #rotation
    d = {'Tilt' : xaxis, 'Rotation' : yaxis}
    coordsDF=pd.DataFrame(d)
    return name, coordsDF


# In[10]:


# Perpendicular to TD
def RingPerpTD(res):
    name="Ring Perpendicular to TD"
    xaxis=np.ndarray.tolist(np.arange(0.0, 90.001, res))#tilt 
    #yaxis=[0.0] * len(xaxis) #rotation
    yaxis=[180.0] * len(xaxis) #rotation
    d = {'Tilt' : xaxis, 'Rotation' : yaxis}
    coordsDF=pd.DataFrame(d)
    return name, coordsDF


# ## Series of rotated rings 

# In[34]:


# Perpendicular to RD
def RingRot(res,theta,omega_list,rotaxis): #add omega
    name="Ring Perpendicular to "+rotaxis
    #name="Ring Perpendicular to ND"
    #print rotaxis
    
    #Generate the ring of rotations
    yaxis=np.ndarray.tolist(np.arange(0.0, 360.0001, res))#rotation
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
        for i, value in enumerate(yaxis):
            #print "Creating Point on ring", i
            if rotaxis=='X':
                #print "Rotating in X axis"
                (r,t)=cart2sphDeg(np.einsum('ij,j', RotateXMatrix(omega_step), sph2cartDeg(value,xaxis[i])))
                # 1- cos**6 works well, but I don't know why...
                weight=(1.0-(np.einsum('i,i',
                                 (np.einsum('ij,j', RotateXMatrix(omega_step), sph2cartDeg(value,xaxis[i]))),
                                 [1,0,0]))**6)
                
                #y=(1/(2*(1-(np.cos(x))**4))) - doesn't work either               
                #weight=(1.0/(2*(1.0-(np.einsum('i,i',
                #                 (np.einsum('ij,j', RotateXMatrix(omega_step), sph2cartDeg(value,xaxis[i]))),
                #                 [1,0,0]))**4)))               

                
                # weight function 22 Oct 2019 - doesn't work
                #weight=(1.0-abs(np.einsum('i,i',
                #                 (np.einsum('ij,j', RotateXMatrix(omega_step), sph2cartDeg(value,xaxis[i]))),
                #                 [1,0,0])))           
                
                
                # Try a different function  - doesn't work well
                #weight=(1.0-2*(np.einsum('i,i',
                #                 (np.einsum('ij,j', RotateXMatrix(omega_step), sph2cartDeg(value,xaxis[i]))),
                #                 [1,0,0])) -
                #        abs(np.einsum('i,i',
                #                 (np.einsum('ij,j', RotateXMatrix(omega_step), sph2cartDeg(value,xaxis[i]))),
                #                 [1,0,0]))**2)                  
            elif rotaxis=='Y':     
                #print "Rotating in Y axis"
                (r,t)=cart2sphDeg(np.einsum('ij,j', RotateYMatrix(omega_step), sph2cartDeg(value,xaxis[i])))
                
                # Try area of a zone in a sphere
                if (value!=90.0):
                    x,y1,z=sph2cartDeg(value+res/2.0,xaxis[i])
                    x,y2,z=sph2cartDeg(value-res/2.0,xaxis[i])
                    weight=0.5*abs(y1-y2)
                else:
                    x,y1,z=sph2cartDeg(value,xaxis[i]) # the segment curves around                 
                    x,y2,z=sph2cartDeg(value-res/2.0,xaxis[i])
                    weight=0.5*abs(y1-y2)                   

                # 1- cos**6 works well, but I don't know why...
#                 weight=(1.0-(np.einsum('i,i',
#                                  (np.einsum('ij,j', RotateXMatrix(omega_step), sph2cartDeg(value,xaxis[i]))),
#                                  [0,1,0]))**6)   

                #y=(1/(2*(1-(np.cos(x))**4))) - doesn't work               
                #weight=(1.0/(2*(1.0-(np.einsum('i,i',
                #                 (np.einsum('ij,j', RotateXMatrix(omega_step), sph2cartDeg(value,xaxis[i]))),
                #                 [0,1,0]))**4)))                     
                
                # weight function 22 Oct 2019 - doesn't work
                #weight=(1.0-abs(np.einsum('i,i',
                #                 (np.einsum('ij,j', RotateXMatrix(omega_step), sph2cartDeg(value,xaxis[i]))),
                #                 [0,1,0])))  
                
                #weight=(1.0-2*(np.einsum('i,i',
                #                 (np.einsum('ij,j', RotateXMatrix(omega_step), sph2cartDeg(value,xaxis[i]))),
                #                 [0,1,0])) -
                #        abs(np.einsum('i,i',
                #                 (np.einsum('ij,j', RotateXMatrix(omega_step), sph2cartDeg(value,xaxis[i]))),
                #                 [0,1,0]))**2)s
                
                #using this as an index on the rotation
                # also fix omega step below                 
                
                
            else:
                print "Not a supported rotation axis"
            
            #(r,t)=StepRotateRing(value,xaxis[i],omega_step)
            #print value, xaxis[i], omega_step,r,t
            
            #print "Adjusting angles"
            # Add some logic to keep 0<=r<360 and 0<=t<90
            if r<0:
                r=r+360
            elif r>=360:
                r=r-360

            # correct for tilt below the sphere.  Don't need to correct rotation here
            if t>90:
                t=180-t
            
            #print "Appending lists"
            r_list.append(r)
            t_list.append(t)
            rotation_index.append(i)   
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
        if omega_counter==0:
            # initialize dataframe
            #print "Intialize Data Frame at x=",omega_counter
            d2 = {'Tilt' : t_list, 'Rotation' : r_list,
                  'RotationIndex' : rotation_index,'Weights' : weight_list}
            coordsDF=pd.DataFrame(d2)
            
        else:    
            d2 = {'Tilt' : t_list, 'Rotation' : r_list,
                  'RotationIndex' : rotation_index, 'Weights' : weight_list}
            d2DF=pd.DataFrame(d2)
            coordsDF=coordsDF.append(d2DF, ignore_index=True)
        
        
        
        
    # range and step of rotations
    
    # append
    #df2 = pd.DataFrame([[5, 6], [7, 8]], columns=list('AB'))
    #df.append(df2)
    
    return name, coordsDF


# In[12]:


# adapted from Steronet
def sph2cartDeg(rotation_deg, tilt_deg):
    """
    Converts a longitude and latitude (or sequence of lons and lats) given in
    _radians_ to cartesian coordinates, `x`, `y`, `z`, where x=0, y=0, z=0 is
    """
    elevation_deg=90-tilt_deg
    elevation=elevation_deg*math.pi/180
    rotation=rotation_deg*math.pi/180
    
    x = np.cos(elevation) * np.cos(rotation)
    y = np.cos(elevation) * np.sin(rotation)
    z = np.sin(elevation)

    return [x, y, z]
    


# In[13]:


# adjust the range of cart2sph, currenlty getting negaitve values of rotate and tilt greater than 90


# In[14]:


def cart2sphDeg(a):
    """

    """

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


# In[15]:


#### Now add rotation around X axis between these


# In[16]:


def RotateYMatrix(omega_deg):
    omega=omega_deg*math.pi/180
    R=[[math.cos(omega),0, -math.sin(omega)  ],
        [0, 1,0],
        [math.sin(omega), 0, math.cos(omega)]]
    return np.transpose(np.array(R)) # transpose to change from active to passive


# In[17]:


def RotateXMatrix(omega_deg):
    omega=np.radians(omega_deg)
    R=[[1,0, 0  ],
        [0, math.cos(omega),math.sin(omega)],
        [0,-math.sin(omega), math.cos(omega)]]
    return np.transpose(np.array(R)) # transpose to change from active to passive


# In[18]:


a=sph2cartDeg(5,1)  #something goes awry when tilt =0


print a
B=RotateYMatrix(1)

c= np.einsum('ij,j', B, a)


cart2sphDeg(c)


# In[ ]:


### Summ rotations


# ## Generate average intensity based on pole figures and coordinates
# - This section calculates the average intensity and saves to file
# - Looks for the list of XPC files and calculates a table (.xlsx) for each
# 
# _Skip if you only wish to plot the sampling schemes_

# In[ ]:


#SchemeName,Coordinates=

# Get the current working directory path
cwd=os.getcwd()
#print cwd
xpcdatapath=os.path.abspath(os.path.join(os.path.dirname(cwd)))
#print xpcdatapath

Folder=os.path.join(os.path.dirname(cwd),'MAUD','XPCFiles')
print Folder

SaveFolder="AveragedIntensites"

if not os.path.isdir(SaveFolder):
    os.makedirs(SaveFolder)
    
os.chdir(SaveFolder)


for file in os.listdir(Folder):
    print file
    if file.endswith(".xpc"):
        XPCfile=(os.path.join(Folder, file))
        
        #
        if "-" in file: 
            orientation, hw=file.split('-')
        else:
            orientation, ext=file.split('.')
            
        PhaseType= orientation[-1:]


    #for XPCfile in listoffiles:

        (pfs,hkllist)=xpcformat('xpc',XPCfile)

        #create subsets for phase fractions

        hkllist.append('2Pairs')
        hkllist.append('4Pairs')
        hkllist.append('MaxUnique')

        OutputList=[hkllist]

        for q in list([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]):

            if q==1: SchemeName,Coordinates=SingleOrientation("ND Single", 0.0,0.0)
            elif q==2: SchemeName,Coordinates=SingleOrientation("RD Single", 90.0,0.0)
            elif q==3: SchemeName,Coordinates=SingleOrientation("TD Single", 90.0,90.0)
            elif q==4: SchemeName,Coordinates=SingleOrientation("Morris", 60.0,90.0)
                ## Add other orientation

            elif q==5: SchemeName,Coordinates=RingPerpND(5.0)
            elif q==6: SchemeName,Coordinates=RingPerpRD(5.0)
            elif q==7: SchemeName,Coordinates=RingPerpTD(5.0)
            #Gaussian Quadrature
            elif q==8: SchemeName,Coordinates=GaussQuad("Gaussian-Quadrature")
            #Spiral Scheme
            elif q==9: SchemeName,Coordinates=Spiral("Spiral-90deg-10rot-5space",90,10,5)
            elif q==10: SchemeName,Coordinates=Spiral("Spiral-90deg-15rot-5space",90,15,5)
            elif q==11: SchemeName,Coordinates=Spiral("Spiral-90deg-10rot-3space",90,10,3)
                # re-arranged to match plotting
            elif q==12: SchemeName,Coordinates=TiltRotate("NoRotation-tilt60deg", 120.0, 1600.0, 0.0,60.0,56.0)
            elif q==13: SchemeName,Coordinates=TiltRotate("Rotation-NoTilt", 120.0, 1600.0, 30.0,0.0,56.0)
            elif q==14: SchemeName,Coordinates=TiltRotate("Rotation-60detTilt", 120.0, 5000.0, 30.0,60.0,56.0)
                
            elif q==15: SchemeName,Coordinates=HexGrid("HexGrid-90degTilt5degRes",90.0,5.0)
            elif q==16: SchemeName,Coordinates=HexGrid("HexGrid-90degTilt22p5degRes",90.0,22.5)
            elif q==17: SchemeName,Coordinates=HexGrid("HexGrid-60degTilt5degRes",60.0,5.0)
                #New testing Schemes start here
                
            else: 
                print "No Schemes of that index"
                break

            PfIS=pfIntensitySum(SchemeName,pfs,Coordinates)
            PfIS.append(np.mean([PfIS[2],PfIS[3]]))
            PfIS.append(np.mean([PfIS[1],PfIS[2],PfIS[3],PfIS[4]]))

            if PhaseType=="A":
                PfIS.append(np.mean([PfIS[1],PfIS[2],PfIS[3],PfIS[4],PfIS[7],PfIS[8],PfIS[9],PfIS[10]]))
            elif PhaseType=="F":
                PfIS.append(np.mean([PfIS[1],PfIS[2],PfIS[3],PfIS[5],PfIS[6],PfIS[7]])) 
                #used to include 4, exclude 5, but that is incorrect
            else:
                print "Unrecognized Phase"

            OutputList.append(PfIS)  

            #print q,SchemeName
            #print "List of average pole Figure Intensities:\n", PfIS
            #print "Average of all pole figures listed: ", sum(PfIS)/len(PfIS)
            #print ""

            #np.mean()

        directory, outfile=XPCfile.rsplit('\\', 1)
        outfile2 = outfile.split('-')
        print outfile
        #OutputList.append([XPCfile]) 
        IntensitiesDF=pd.DataFrame(OutputList)
        IntensitiesDF

        # save to excel

        s=""     

        writer = pd.ExcelWriter(s.join([outfile.split('.')[0], ".xlsx"]))
        IntensitiesDF.to_excel(writer,outfile2[0])
        writer.save()
    else:
        print "Not an .xpc file"
os.chdir("..")


# # Plots

# # Plot pole figures of sampling positions

# ### Test functions

# In[ ]:


SchemeName,Coordinates=SingleOrientation("Morris", 60.0,90.0)
Coordinates


# ### Simple plot to work out mplsteronet conventions

# In[ ]:


fig = plt.figure(figsize=(8,8), dpi=600)

ax2 = fig.add_subplot(111, projection='stereonet')

#SingleOrientation - Name, Tilt, Rotation
#SchemeName,Coordinates=SingleOrientation("RD Single", 90.0,180.0)
SchemeName,Coordinates=SingleOrientation("TD Single", 90.0,270.0)

dip, strike =Coordinates['Tilt'], Coordinates['Rotation']-90.0
l1=ax2.pole(strike, dip, 'bD', markersize=10, clip_on=False)
#dip tilts about the 0 axis (RD), left handed
#strike tilts about the normal axis (ND), left handed
#Looking at the bottom of the sphere, not the top
#this is just convention of mplstereonet, does not affect averaging methods

plt.show()


# ### Simple plot to test new sampling schemes

# In[ ]:


fig = plt.figure(figsize=(8,8), dpi=600)

ax10 = fig.add_subplot(111, projection='stereonet')

ax10.plane(0.0, 90.0, 'k-', linewidth=1)
ax10.plane(90.0, 90.0, 'k-', linewidth=1)

SchemeName,Coordinates=Spiral("Spiral",90,10,3)
dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)
l1=ax10.pole(strike, dip, marker ='o',mfc='purple', markersize=10, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)

plt.show()


# In[ ]:


fig = plt.figure(figsize=(8,8), dpi=600)

ax10 = fig.add_subplot(111, projection='stereonet')

ax10.plane(0.0, 90.0, 'k-', linewidth=1)
ax10.plane(90.0, 90.0, 'k-', linewidth=1)

SchemeName,Coordinates=HexGrid("HexGrid-90degTilt5degRes",90.0,10.0)
dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)
l1=ax10.pole(strike, dip, marker ='o',mfc='purple', markersize=10, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)

plt.show()


# ## Rotated ring plotting

# In[32]:


fig = plt.figure(figsize=(8,8), dpi=600)

ax10 = fig.add_subplot(111, projection='stereonet')

ax10.plane(0.0, 90.0, 'k-', linewidth=1)
ax10.plane(90.0, 90.0, 'k-', linewidth=1)

#RingRotRD
res=10
theta=2.5

#omega_list=np.arange(-90, 91, 5 )
#omega_list=np.arange(-90, 1, 22.5 )
omega_list=[0,10]

#l1=ax10.pole(0, 0, marker ="s",mfc='k', markersize=10, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)

SchemeName,Coordinates=RingRot(res,theta,omega_list,'Y')
#print Coordinates
dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)
l1=ax10.pole(strike, dip, marker ='o',mfc='magenta', markersize=5, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)

omega_list=[]
#SchemeName,Coordinates=RingRot(res,theta,omega_list,'X')
#print Coordinates
#dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)
#l1=ax10.pole(strike, dip, marker ='o',mfc='k', markersize=5, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)
#l1=ax10.plane(-90, 10, marker ="+",mfc='purple', markersize=10, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)
#l1=ax10.pole(-90, 10, marker ="x",mfc='purple', markersize=10, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)

#plt.savefig("SingleRing-Theta2p5-90deg.png", dpi=600,format="png")
plt.show()


# In[31]:


sph2cartDeg(0,0)


# In[ ]:





# In[36]:


#RingRotRD
res=5
theta=5

#omega_list=np.arange(-90, 91, 5 )
#omega_list=np.arange(-90, 1, 22.5 )
omega_list=[0]

#l1=ax10.pole(0, 0, marker ="s",mfc='k', markersize=10, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)

SchemeName,Coordinates=RingRot(res,theta,omega_list,'Y')
#df['col_3'] = df.apply(lambda x: f(x.col_1, x.col_2), axis=1)
#Coordinates['Cart']=Coordinates.apply(lambda row: sph2cartDeg(row['Rotation'],rwo['Tilt']), axis=1)
#(Coordinates['X'],Coordinates['Y'],Coordinates['Z'])=Coordinates.apply(lambda row: sph2cartDeg(row['Rotation'],
#                                                                                             row['Tilt']), axis=1)

X=Coordinates.apply(lambda row: pd.Series(sph2cartDeg(row['Rotation'],row['Tilt']),index=['x','y','z']), axis=1)
#X
Coordinates=pd.concat([Coordinates, X], axis=1)

print Coordinates

# Just print the coordinates of a specific rotation
#print Coordinates.loc[Coordinates['Weights'] == 3]
#sph2cartDeg

#plot these
fig = plt.figure(figsize=(8,8), dpi=600)

ax10 = fig.add_subplot(111, projection='stereonet')

ax10.plane(0.0, 90.0, 'k-', linewidth=1)
ax10.plane(90.0, 90.0, 'k-', linewidth=1)
dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)
l1=ax10.pole(strike, dip, marker ='o',mfc='magenta', markersize=5, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)
plt.show()


fig = plt.figure(figsize=(8,8), dpi=600)
plt.scatter(Coordinates['RotationIndex'],Coordinates['y'])
plt.show()


# In[ ]:





# In[ ]:


Coordinates['Weights']


# In[38]:


from matplotlib.colors import ListedColormap, LinearSegmentedColormap

n=6

jet = plt.cm.get_cmap('jet', 100)
newcolors = jet(np.linspace(0, 1, n+1))

gray = np.array([166.0/256.0, 166.0/256.0, 166.0/256.0, 1])
newcolors[:1, :] = gray
print newcolors
newcmp = ListedColormap(newcolors)

v = np.linspace(0, n, n+1)
#plt.contourf(phi1list, PHIlist, Z,v, cmap='jet')
#print v


# In[39]:


### Using Colormap4, adding alpha
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
         


# In[42]:


fig = plt.figure(figsize=(12,6), dpi=600)



#RingRotRD
res=5
theta=10

omega_list=np.arange(-180, 180, 5 )
#omega_list=np.arange(-90, 91, 5 )
#omega_list=np.arange(-90, 86, 5 )
#omega_list=np.arange(-85, 86, 5 )

#l1=ax10.pole(0, 0, marker ="s",mfc='k', markersize=10, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)

ax09 = fig.add_subplot(231, projection='stereonet')
ax09.set_azimuth_ticks([90,0], labels=['','']) 
ax09.annotate('Ring Rotate', xy=(0, 0), xytext=(-4.0,-1.0), fontsize=16)
#Ring Rotation
SchemeName,Coordinates=RingRot(res,theta,omega_list,'X')
dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)

cax = ax09.density_contourf(strike, dip,levels=v, measurement='poles',method='schmidt',cmap=newcmp2 )


ax10 = fig.add_subplot(232, projection='stereonet')
ax10.set_azimuth_ticks([90,0], labels=['','']) 
ax10.annotate('Ring Rotate \n Weighted', xy=(0, 0), xytext=(-4.0,-1.0), fontsize=16)
#Ring Rotation
SchemeName,Coordinates=RingRot(res,theta,omega_list,'X')
dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)

cax = ax10.density_contourf(strike, dip,levels=v, measurement='poles',
                            method='schmidt',gridsize=500,cmap=newcmp2, weights=list(Coordinates['Weights']) )

ax13 = fig.add_subplot(233, projection='stereonet')
ax13.set_azimuth_ticks([90,0], labels=['','']) 
#ax13.plane(0.0, 90.0, 'k-', linewidth=1)
#ax13.plane(90.0, 90.0, 'k-', linewidth=1)
fig.colorbar(cax)

ax11 = fig.add_subplot(234, projection='stereonet')
ax11.set_azimuth_ticks([90,0], labels=['','']) 
ax11.annotate('Hex Grid', xy=(0, 0), xytext=(-4.0,-1.0), fontsize=16)
#Hex Grid
SchemeName,Coordinates=HexGrid("HexGrid-90degTilt5degRes",90.0,5.0)
dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)
cax = ax11.density_contourf(strike, dip,levels=v, measurement='poles',method='schmidt',gridsize=100,cmap=newcmp2)
#l1=ax11.pole(strike, dip, marker ='o',mfc='k', markersize=2, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)


ax12 = fig.add_subplot(235, projection='stereonet')
ax12.set_azimuth_ticks([90,0], labels=['',''])
ax12.annotate('Tilt & Rotate', xy=(0, 0), xytext=(-4.0,-1.0), fontsize=16)
# Tilt & Rotate
SchemeName,Coordinates=TiltRotate("Rotation-60detTilt", 120.0, 5000.0, 30.0,60.0,56.0)
dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)
cax = ax12.density_contourf(strike, dip,levels=v, measurement='poles',method='schmidt',cmap=newcmp2)


#l1=ax10.pole(strike, dip, marker ='o',mfc='k', markersize=2, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)

plt.show()


# In[36]:


print len(v)


# ## Plot the schemes with diffraction vectors and density of points

# In[68]:


fig = plt.figure(figsize=(12,6), dpi=600)

# Tilt & Rotate
fig.suptitle("Tilt & Rotate", fontsize=16)


SchemeName,Coordinates=TiltRotate("Rotation-60detTilt", 120.0, 5000.0, 30.0,60.0,56.0)
dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)

ax1 = fig.add_subplot(121, projection='stereonet')
ax1.set_title('Diffraction Vectors')
ax1.set_azimuth_ticks([90,0], labels=['','']) 
l1=ax1.pole(strike, dip, marker ='o',mfc='k', markersize=2, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)
#fig.colorbar(cax,orientation='horizontal')


ax2 = fig.add_subplot(122, projection='stereonet')
ax2.set_title('Density of points')
ax2.set_azimuth_ticks([90,0], labels=['','']) 
#ax2.annotate('Hex Grid', xy=(0, 0), xytext=(-4.0,-1.0), fontsize=16)
#cax = ax11.density_contourf(strike, dip,levels=v, measurement='poles',method='schmidt',gridsize=100,cmap=newcmp2)
cax = ax2.density_contourf(strike, dip,levels=v, measurement='poles',method='schmidt',cmap=newcmp2 )
cbar=fig.colorbar(cax,orientation='horizontal')
cbar.ax.set_xlabel("Oversampling Multiples")


#plt.savefig("TiltRotate_Diff_Density.png", dpi=600,format="png")
plt.show()


# In[66]:


fig = plt.figure(figsize=(12,6), dpi=600)

#RingRotRD
res=5
theta=10

omega_list=np.arange(-90, 91, 5 )
#omega_list=np.arange(-90, 86, 5 )
#omega_list=np.arange(-85, 86, 5 )

#ax09 = fig.add_subplot(111, projection='stereonet')
ax09.set_azimuth_ticks([90,0], labels=['','']) 

#Ring Rotation
#ax09.annotate('Ring Rotate', xy=(0, 0), xytext=(-4.0,-1.0), fontsize=16)
#SchemeName,Coordinates=RingRot(res,theta,omega_list,'X')
#dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)

#Hex Grid
fig.suptitle("Hex Grid", fontsize=16)


SchemeName,Coordinates=HexGrid("HexGrid-90degTilt5degRes",90.0,5.0)
dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)

ax1 = fig.add_subplot(121, projection='stereonet')
ax1.set_title('Diffraction Vectors')
ax1.set_azimuth_ticks([90,0], labels=['','']) 
l1=ax1.pole(strike, dip, marker ='o',mfc='k', markersize=2, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)
#fig.colorbar(cax,orientation='horizontal')


ax2 = fig.add_subplot(122, projection='stereonet')
ax2.set_title('Density of points')
ax2.set_azimuth_ticks([90,0], labels=['','']) 
#ax2.annotate('Hex Grid', xy=(0, 0), xytext=(-4.0,-1.0), fontsize=16)
#cax = ax11.density_contourf(strike, dip,levels=v, measurement='poles',method='schmidt',gridsize=100,cmap=newcmp2)
cax = ax2.density_contourf(strike, dip,levels=v, measurement='poles',method='schmidt',cmap=newcmp2 )
cbar=fig.colorbar(cax,orientation='horizontal')
cbar.ax.set_xlabel("Oversampling Multiples")


#plt.savefig("HexGrid_Diff_Density.png", dpi=600,format="png")
plt.show()


# In[70]:


fig = plt.figure(figsize=(12,6), dpi=600)

#RingRotRD

#Hex Grid
fig.suptitle("Ring Rotate", fontsize=16)

res=5
theta=5
omega_list=np.arange(-90, 91, 5 )

SchemeName,Coordinates=RingRot(res,theta,omega_list,'X')
dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)

ax1 = fig.add_subplot(121, projection='stereonet')
ax1.set_title('Diffraction Vectors')
ax1.set_azimuth_ticks([90,0], labels=['','']) 
l1=ax1.pole(strike, dip, marker ='o',mfc='k', markersize=2, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)
#fig.colorbar(cax,orientation='horizontal')


ax2 = fig.add_subplot(122, projection='stereonet')
ax2.set_title('Density of points')
ax2.set_azimuth_ticks([90,0], labels=['','']) 
#ax2.annotate('Hex Grid', xy=(0, 0), xytext=(-4.0,-1.0), fontsize=16)
#cax = ax11.density_contourf(strike, dip,levels=v, measurement='poles',method='schmidt',gridsize=100,cmap=newcmp2)
cax = ax2.density_contourf(strike, dip,levels=v, measurement='poles',method='schmidt',cmap=newcmp2 )
cbar=fig.colorbar(cax,orientation='horizontal')
cbar.ax.set_xlabel("Oversampling Multiples")


#plt.savefig("RingRotateX_Diff_Density.png", dpi=600,format="png")
plt.show()


# In[72]:


## Interpolate and plot the directions


res=5
theta=10
omega_list=np.arange(-90, 91, 5 )
SchemeName,Coordinates=RingRot(res,theta,omega_list,'X')
dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)

#SchemeName,Coordinates=HexGrid("HexGrid-90degTilt5degRes",90.0,5.0)
#dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)

print len(dip)

MPLgridsize=100

# always returns 100x100 - set by grid size, default 100
a=mplstereonet.density_grid(strike, dip, method='schmidt', gridsize=MPLgridsize)

######
### Need to add interpolation step
######

aDF=pd.DataFrame({'strike': np.degrees(a[0].flatten()), 'dip': np.degrees(a[1].flatten()),'density': a[2].flatten()}  )

aDF[['CartX','CartY','CartZ']]=aDF.apply(lambda x: sph2cartDeg(x['strike'],x['dip']),axis=1)

#for i,value in aDF['dip']:
#    aDF['Cart'][i]=sph2cartDeg(value,aDF['strike'][i])

#print aDF

# reminder: dip=tilt, strike=rotation

# plotting the filled contour
fig = plt.figure(figsize=(12,4), dpi=600)

ax1 = fig.add_subplot(131, projection='stereonet')
cax = ax1.density_contourf(strike, dip,levels=v, measurement='poles',
                           method='schmidt',gridsize=MPLgridsize,cmap=newcmp2 )
fig.colorbar(cax)

# plot the grid
ax2 = fig.add_subplot(132, projection='stereonet')
l1=ax2.pole(strike, dip, marker ='o',mfc='k', markersize=2, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)
l1=ax2.pole(np.degrees(a[0].flatten()), np.degrees(a[1].flatten()),
            marker ='+',mfc='c', markersize=10, clip_on=False)

# plot the weighting factor
ax3 = fig.add_subplot(133)
ax3.set_ylim([0,10])
#ax3.scatter((aDF['CartX']**2+aDF['CartY']**2+aDF['CartZ']**2), aDF['density'])
#ax3.scatter(np.degrees(np.cos(aDF['CartX'])), aDF['density'])
ax3.scatter(aDF['CartZ'], aDF['density']) # plots what I expect to see
#ax3.scatter((1-aDF['CartZ']**2-aDF['CartY']**2), aDF['density']) # odd plot

x=np.arange(0, 1, 0.01 )
#y=np.ones(len(x))
y=(1/(2*(1-(np.cos(x))**4))) # matches curve fairly well at x>0.3

#y=1/(np.arcsin(x))  # also works fairly well at x>0.3
#y=1/(np.arc(x))


ax3.plot(x,y,c='r')


# Doesn look like a log norm
#from scipy.stats import lognorm
#stddev = 0.1
#mean = -.50
#dist=lognorm([stddev],loc=mean)

#ax3.plot(x,dist.pdf(x))


ax3.axhline(1.0, c='y')
ax3.axvline(np.sin(np.radians(theta)), c='c') # at theta

plt.show()


# In[79]:


fig = plt.figure(figsize=(12,6), dpi=600)

#RingRotRD

fig.suptitle("Ring Rotate", fontsize=16)

res=5
theta=5
omega_list=np.arange(-90, 91, 5 )

SchemeName,Coordinates=RingRot(res,theta,omega_list,'X')
dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)

#ax1 = fig.add_subplot(121, projection='stereonet')
#ax1.set_title('Diffraction Vectors')
#ax1.set_azimuth_ticks([90,0], labels=['','']) 
#l1=ax1.pole(strike, dip, marker ='o',mfc='k', markersize=2, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)
#fig.colorbar(cax,orientation='horizontal')


ax2 = fig.add_subplot(121, projection='stereonet')
ax2.set_title('Density of points')
ax2.set_azimuth_ticks([90,0], labels=['','']) 
#ax2.annotate('Hex Grid', xy=(0, 0), xytext=(-4.0,-1.0), fontsize=16)
#cax = ax11.density_contourf(strike, dip,levels=v, measurement='poles',method='schmidt',gridsize=100,cmap=newcmp2)
cax = ax2.density_contourf(strike, dip,levels=v, measurement='poles',method='schmidt',cmap=newcmp2 )
cbar=fig.colorbar(cax,orientation='horizontal')
cbar.ax.set_xlabel("Oversampling Multiples")

ax3 = fig.add_subplot(122, projection='stereonet')
ax3.set_title('Weighted Density of points')
ax3.set_azimuth_ticks([90,0], labels=['','']) 
#ax3.annotate('Ring Rotate \n Weighted', xy=(0, 0), xytext=(-4.0,-1.0), fontsize=16)
#Ring Rotation

cax3 = ax3.density_contourf(strike, dip,levels=v, measurement='poles',
                            method='schmidt',gridsize=500,cmap=newcmp2, weights=list(Coordinates['Weights']) )
cbar=fig.colorbar(cax3,orientation='horizontal')
cbar.ax.set_xlabel("Oversampling Multiples")



plt.savefig("RingRotateX_Density_Weighted.png", dpi=600,format="png")
plt.show()


# In[ ]:





# In[ ]:




#interpolation function
f=scipy.interpolate.Rbf(a[0], a[1],a[2])

#sph2cartDeg(value,xaxis[i])

#print np.array([[strike], [dip],[f(strike, dip)]])
print f(strike, dip)
#print strike.shape
#print dip.shape
#print f(strike, dip).shape

fig = plt.figure(figsize=(12,6), dpi=600)
ax09 = fig.add_subplot(111, projection='stereonet')
ax09.set_azimuth_ticks([90,0], labels=['','']) 
ax09.annotate('Test', xy=(0, 0), xytext=(-4.0,-1.0), fontsize=16)
#Ring Rotation
#SchemeName,Coordinates=RingRot(res,theta,omega_list,'X')
#dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)

cax = ax09.density_contourf(strike, dip,levels=v, measurement='poles',method='schmidt',cmap=newcmp )
fig.colorbar(cax)
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


def fullprint(*args, **kwargs):
    from pprint import pprint
    import numpy
    opt = numpy.get_printoptions()
    numpy.set_printoptions(threshold='nan')
    numpy.set_printoptions(precision=3)
    np.set_printoptions(suppress=True)
    pprint(*args, **kwargs)
    numpy.set_printoptions(**opt)


# In[ ]:


SchemeName,Coordinates=HexGrid("HexGrid-90degTilt5degRes",90.0,1.0)
dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)

#SchemeName,Coordinates=RingRot(res,theta,omega_list,'X')
#dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)

a=mplstereonet.density_grid(strike, dip,method='schmidt')
fullprint(a[2])


# In[ ]:


# how much overcoverage is there
# shows I need to adjust the range of tilt & rotate


# In[ ]:


#with pd.option_context('display.max_rows', None, 'display.max_columns', None,'display.float_format', lambda x: '%.3f' % x):  # more options can be specified also
#    print Coordinates


# In[ ]:


fig = plt.figure(figsize=(12,3), dpi=600)
res=5
theta=10
omega_list=np.arange(-90, 91, 5 )
SchemeName,Coordinates=RingRot(res,theta,omega_list,'X')

h=plt.hist2d(Coordinates['Rotation'], Coordinates['Tilt'], bins=[72, 18], range=[[0, 360.1], [0, 90.1]])

plt.colorbar(h[3])

plt.axhline(90,c='white')

plt.show()


# In[ ]:


print h


# In[ ]:


## expanded range to show additional sampling at tilt =90 


# In[ ]:


res=5
theta=10
omega_list=np.arange(-90, 91, 5 )
SchemeName,Coordinates=RingRot(res,theta,omega_list,'X')

h=plt.hist2d(Coordinates['Rotation'], Coordinates['Tilt'], bins=[54, 18], range=[[-180, 360], [0, 180]])

plt.colorbar(h[3])

plt.show()


# In[ ]:





# ## Plot all the Sampling methods detailed in the Paper

# In[ ]:


# Note - MPLStereonet uses a different angle convention (geosciences) than is used for crystallography
# A rotation of 90° is needed to align the coordiantes
# TD, Morris rotated by 270° to be in the same quadrant as Matlab
#dip = (Coordinates['Rotation']-90.0)


# 3 x 3 subplot grid
# single, ring, key
# Tilt/Rotate
# Hex Grids

fig = plt.figure(figsize=(8,9), dpi=600)


#key
ax1 = fig.add_subplot(331, projection='stereonet')
ax1.set_azimuth_ticks([0,90], labels=['RD','-- TD'],fontsize=14) 
ax1.plane(0.0, 90.0, 'k-', linewidth=1)
ax1.plane(90.0, 90.0, 'k-', linewidth=1)
ax1.annotate('(a)', xy=(0, 0), xytext=(-3.0,0.1), fontsize=16)
ax1.annotate('ND', xy=(0, 0), xytext=(0.1,0.1), fontsize=14)


#Single orientations

ax2 = fig.add_subplot(332, projection='stereonet')

ax2.plane(0.0, 90.0, 'k-', linewidth=1)
ax2.plane(90.0, 90.0, 'k-', linewidth=1)
            #if q==1: SchemeName,Coordinates=SingleOrientation("ND Single", 0.0,0.0)
            #elif q==2: SchemeName,Coordinates=SingleOrientation("RD Single", 90.0,0.0)
            #elif q==3: SchemeName,Coordinates=SingleOrientation("TD Single", 90.0,90.0)
            #elif q==4: SchemeName,Coordinates=SingleOrientation("Morris", 60.0,90.0)
                
SchemeName,Coordinates=SingleOrientation("ND Single", 0.0,0.0)
dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)
l1=ax2.pole(strike, dip, 'bs', markersize=10, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)

SchemeName,Coordinates=SingleOrientation("RD Single", 90.0,180.0)
dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)
l2=ax2.pole(strike, dip, 'rs', markersize=10, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)

SchemeName,Coordinates=SingleOrientation("TD Single", 90.0,270.0)
dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)
l3=ax2.pole(strike, dip, 'gs', markersize=10, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)

SchemeName,Coordinates=SingleOrientation("Morris", 60.0,270.0)
dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)
l4=ax2.pole(strike, dip, 'ys', markersize=10, clip_on=False, markeredgecolor='black', markeredgewidth=0.5)

ax2.set_azimuth_ticks([90,0], labels=['',''])
ax2.annotate('(b)', xy=(0, 0), xytext=(-3.0,0.1), fontsize=16)

#ax1.grid()
#for index, row in Coordinates.iterrows():
    #print row['Tilt'], row['Rotation']
    #dip, strike =  row['Tilt'], row['Rotation']
    #ax1.pole(strike, dip, 'g^', markersize=5)
    
    
    #ax.plane(strike, dip, 'g-', linewidth=2)
    #ax.rake(strike, dip, -70)
    
# Ring Orientations
ax3 = fig.add_subplot(333, projection='stereonet')
ax3.annotate('(c)', xy=(0, 0), xytext=(-3.0,0.1), fontsize=16)
ax3.set_azimuth_ticks([90,0], labels=['','']) 


SchemeName,Coordinates=RingPerpND(5.0)
ax3.plane(0.0, 0.0, 'b-', linewidth=3)
ax3.plane(0.0, 180.0, 'b-', linewidth=3)
for index, row in Coordinates.iterrows():
    #print row['Tilt'], row['Rotation']
    dip, strike =  row['Tilt'], (row['Rotation']-90.0)
    ax3.pole(strike, dip, 'k+', markersize=10, clip_on=False)



SchemeName,Coordinates=RingPerpRD(5.0)
ax3.plane(-90.0, 90.0, 'r-', linewidth=3) #0.0 (RD)-90.0 (Strike Converntion ) = -90.0Dip
for index, row in Coordinates.iterrows():
    #print row['Tilt'], row['Rotation']
    dip, strike =  row['Tilt'], (row['Rotation']-90.0)
    ax3.pole(strike, dip, 'r+', markersize=10, clip_on=False)


SchemeName,Coordinates=RingPerpTD(5.0)
ax3.plane(0.0, 90.0, 'g-', linewidth=3)#90.0 (TD)-90.0 (Strike Converntion ) = 0.0Dip
for index, row in Coordinates.iterrows():
    #print row['Tilt'], row['Rotation']
    dip, strike =  row['Tilt'], (row['Rotation']-90.0)
    ax3.pole(strike, dip, 'g+', markersize=10, clip_on=False)






ax4 = fig.add_subplot(334, projection='stereonet')
ax4.annotate('(d)', xy=(0, 0), xytext=(-3.0,0.1), fontsize=16)
ax4.set_azimuth_ticks([90,0], labels=['','']) 
ax4.plane(0.0, 90.0, 'k-', linewidth=1)
ax4.plane(90.0, 90.0, 'k-', linewidth=1)
SchemeName,Coordinates=TiltRotate("Rotation-NoTilt", 120.0, 1600.0, 30.0,0.0,56.0)
for index, row in Coordinates.iterrows():
    #print row['Tilt'], row['Rotation']
    dip, strike =  row['Tilt'], (row['Rotation']-90.0)
    ax4.pole(strike, dip, 'cD', markersize=5, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)



ax5 = fig.add_subplot(335, projection='stereonet')
ax5.annotate('(e)', xy=(0, 0), xytext=(-3.0,0.1), fontsize=16)
ax5.plane(0.0, 90.0, 'k-', linewidth=1)
ax5.plane(90.0, 90.0, 'k-', linewidth=1)
ax5.set_azimuth_ticks([90,0], labels=['','']) 
SchemeName,Coordinates=TiltRotate("NoRotation-tilt60deg", 120.0, 1600.0, 0.0,60.0,56.0)
for index, row in Coordinates.iterrows():
    #print row['Tilt'], row['Rotation']
    dip, strike =  row['Tilt'], (row['Rotation']-90.0)
    ax5.pole(strike, dip, 'cD', markersize=3, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)


ax6 = fig.add_subplot(336, projection='stereonet')
ax6.annotate('(f)', xy=(0, 0), xytext=(-3.0,0.1), fontsize=16)
ax6.plane(0.0, 90.0, 'k-', linewidth=1)
ax6.plane(90.0, 90.0, 'k-', linewidth=1)
ax6.set_azimuth_ticks([90,0], labels=['','']) 
SchemeName,Coordinates=TiltRotate("Rotation-60detTilt", 120.0, 5000.0, 30.0,60.0,56.0)
for index, row in Coordinates.iterrows():
    #print row['Tilt'], row['Rotation']
    dip, strike =  row['Tilt'], (row['Rotation']-90.0)
    ax6.pole(strike, dip, 'cD', markersize=2, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)


    
# Hex Grids

ax7 = fig.add_subplot(337, projection='stereonet')
ax7.annotate('(g)', xy=(0, 0), xytext=(-3.0,0.1), fontsize=16)
ax7.plane(0.0, 90.0, 'k-', linewidth=1)
ax7.plane(90.0, 90.0, 'k-', linewidth=1)
ax7.set_azimuth_ticks([90,0], labels=['','']) 
SchemeName,Coordinates=HexGrid("HexGrid-90degTilt5degRes",90.0,5.0)
for index, row in Coordinates.iterrows():
    #print row['Tilt'], row['Rotation']
    dip, strike =  row['Tilt'], (row['Rotation']-90.0)
    ax7.pole(strike, dip, 'mh', markersize=2, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)


ax8 = fig.add_subplot(338, projection='stereonet')
ax8.annotate('(h)', xy=(0, 0), xytext=(-3.0,0.1), fontsize=16)
ax8.plane(0.0, 90.0, 'k-', linewidth=1)
ax8.plane(90.0, 90.0, 'k-', linewidth=1)
ax8.set_azimuth_ticks([90,0], labels=['','']) 
SchemeName,Coordinates=HexGrid("HexGrid-90degTilt22p5degRes",90.0,22.5)
for index, row in Coordinates.iterrows():
    #print row['Tilt'], row['Rotation']
    dip, strike =  row['Tilt'], (row['Rotation']-90.0)
    ax8.pole(strike, dip, 'mh', markersize=4, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)



ax9 = fig.add_subplot(339, projection='stereonet')
ax9.annotate('(i)', xy=(0, 0), xytext=(-3.0,0.1), fontsize=16)
ax9.plane(0.0, 90.0, 'k-', linewidth=1)
ax9.plane(90.0, 90.0, 'k-', linewidth=1)
ax9.set_azimuth_ticks([90,0], labels=['','']) 
SchemeName,Coordinates=HexGrid("HexGrid-60degTilt5degRes",60.0,5.0)
for index, row in Coordinates.iterrows():
    #print row['Tilt'], row['Rotation']
    dip, strike =  row['Tilt'], (row['Rotation']-90.0)
    ax9.pole(strike, dip, 'mh', markersize=2, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)
    


# Legend
import matplotlib.lines as mlines
blue_dot = mlines.Line2D([], [], color='white', marker='s', mfc='b',
                          markersize=15, label='',markeredgecolor='black', markeredgewidth=0.5)
red_dot = mlines.Line2D([], [], color='white', marker='s', mfc='r',
                          markersize=15, label='',markeredgecolor='black', markeredgewidth=0.5)
green_dot = mlines.Line2D([], [], color='white', marker='s', mfc='g',
                          markersize=15, label='',markeredgecolor='black', markeredgewidth=0.5)
yellow_dot = mlines.Line2D([], [], color='white', marker='s', mfc='y',
                          markersize=15, label='',markeredgecolor='black', markeredgewidth=0.5)

blue_plus = mlines.Line2D([], [], color='white', marker='+', mec='b', mew=2.0,
                          markersize=15, label='')
red_plus = mlines.Line2D([], [], color='white', marker='+', mec='r',mew=2.0,
                          markersize=15, label='')
green_plus = mlines.Line2D([], [], color='white', marker='+', mec='g',mew=2.0,
                          markersize=15, label='')

cyan_square = mlines.Line2D([], [], color='white', marker='D', mfc='c',
                          markersize=12, label='',markeredgecolor='black', markeredgewidth=0.5)

mag_hex = mlines.Line2D([], [], color='white', marker='h', mfc='m',
                          markersize=15, label='',markeredgecolor='black', markeredgewidth=0.5)

handles=[blue_dot, red_dot, green_dot,  blue_plus, red_plus, green_plus, yellow_dot, cyan_square, mag_hex]
labels=['ND Single','RD Single','TD Single', 'ND Ring','RD Ring','TD Ring','Morris Single', 'Tilt and Rotation', 'Hexagonal Grid']


#plt.figlegend((handles),(labels),'lower center', numpoints=1, ncol=3,fontsize=16,bbox_to_anchor=[0.46, -0.015])
plt.figlegend((handles),(labels),'lower center', numpoints=1, ncol=3,fontsize=16)


#swtich between saving figure and showing the figure

#plt.savefig("SamplingSchemes-Draft.eps", dpi=600,format="eps")
#plt.savefig("SamplingSchemes-Draft.pdf", dpi=600,format="pdf")
#plt.savefig("SamplingSchemes-Draft.png", dpi=600,format="png")
plt.show()


# ## Plotting of all the new Sampling Schemes

# In[ ]:


fig = plt.figure(figsize=(8,9), dpi=600)

#key
ax1 = fig.add_subplot(331, projection='stereonet')
ax1.set_azimuth_ticks([0,90], labels=['RD','-- TD'],fontsize=14) 
ax1.plane(0.0, 90.0, 'k-', linewidth=1)
ax1.plane(90.0, 90.0, 'k-', linewidth=1)
ax1.annotate('(a)', xy=(0, 0), xytext=(-3.0,0.1), fontsize=16)
ax1.annotate('ND', xy=(0, 0), xytext=(0.1,0.1), fontsize=14)

#Gaussian Quadrature
ax10 = fig.add_subplot(332, projection='stereonet')
ax10.annotate('(b)', xy=(0, 0), xytext=(-3.0,0.1), fontsize=16)
ax10.plane(0.0, 90.0, 'k-', linewidth=1)
ax10.plane(90.0, 90.0, 'k-', linewidth=1)
ax10.set_azimuth_ticks([90,0], labels=['','']) 

SchemeName,Coordinates=GaussQuad("Gaussian Quadrature")
dip, strike =Coordinates['Tilt'], (Coordinates['Rotation']-90.0)
l1=ax10.pole(strike, dip, marker ='^',mfc='orange', markersize=10, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)

ax11 = fig.add_subplot(333, projection='stereonet')
ax11.annotate('(c)', xy=(0, 0), xytext=(-3.0,0.1), fontsize=16)
ax11.plane(0.0, 90.0, 'k-', linewidth=1)
ax11.plane(90.0, 90.0, 'k-', linewidth=1)
ax11.set_azimuth_ticks([90,0], labels=['','']) 

SchemeName,Coordinates=Spiral("Spiral",90,10,5)
for index, row in Coordinates.iterrows():
    #print row['Tilt'], row['Rotation']
    dip, strike =  row['Tilt'], (row['Rotation']-90.0)
    ax11.pole(strike, dip, marker = 'o',mfc = 'purple', markersize=5, clip_on=False,markeredgecolor='black', markeredgewidth=0.5)
#legend
import matplotlib.lines as mlines
org_tri= mlines.Line2D([], [], color='white', marker='^', mfc='orange',
                          markersize=15, label='',markeredgecolor='black', markeredgewidth=0.5)
purp_circ= mlines.Line2D([], [], color='white', marker='o', mfc='purple',
                          markersize=15, label='',markeredgecolor='black', markeredgewidth=0.5)
handles = [org_tri, purp_circ]
labels = ['Gaussian Quadrature','Spiral']

plt.figlegend((handles),(labels),'lower center', numpoints=1, ncol=3,fontsize=16)

plt.show()

