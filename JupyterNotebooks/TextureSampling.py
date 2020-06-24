Define function to tilt and rotate:
def TiltRotate(name, time, datapoints, rpm,maxtilt,tiltcpm):
    #time=120  # in seconds
    #datapoints = 5000.0
    """A function used to represent the tilting angles and rotational angles
       of a sample for a given amount of datapoints, and outputs arrays of
       rotationposition and tiltposition that correspond to certain "time
       values"
       
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
    help(TiltRotate)

        #print coordsDF.sort_values('Tilt')    
    return name, coordinates


----------------------------------------------------------------------------------------------------------------------------------------


Define function to create a hexagonal grid:

def HexGrid(name, chi_max, angular_spacing):
 """A function used to create a hexagonal grid scheme, to be used for displaying the crystallographic texture of, in this case, various
 steels. The Hexagonal Grid Scheme has experimentally proven to be effective in the reduction of error associated with measuring Austenite
 Phase Fraction, compared to other potential schemes.   This function accepts the sample name, maximum "chi-angle" for measurements 
 (elaborated upon in the "Parameters" description), and the desired angle amount to increment each data point by. Outputs arrays of
 "Tilt" and "Rotation" angle values.
       
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
        """
    
#chi_max=90.0  #maximum tilt angle in degrees
#angular_spacing=7.0


    d_max=2.0*math.sin(math.radians(chi_max)/2.0) %

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





----------------------------------------------------------------------------------------------------------------------------------------


Define function to create a spiral grid scheme (INCOMPLETE):
import numpy as np
import pandas as pd
import math
# The angluar spacing is determined recursively based on the values of Phi and Chi
#on the run before, so I think these are the appropriate parameters...
def SpiralScheme(name, chi_max, phi_max):
    xaxis=[] #Chi Array of values to be outputted
    yaxis=[] #Phi Array of values to be outputted
    i=1 # Chi Counter, very important for the determination of Phi 
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
        if (phi_N<phi_max):
            chi_N=np.exp(b*phi_N)
        elif (phi_N==phi_max):
            chi_N=chi_max
            xaxis.append(chi_N)
            break
        
        xaxis.append(chi_N)
               
            # Need to perform the translations (arctan, scale down by 1/45 as the research paper, but would like to
            #understand why this needs to be done before implementing it...)
        
    i=i+1
    
    
    values = {'Tilt' : xaxis, 'Rotation' : yaxis}
    coordinates=pd.DataFrame(values)   
    DroppedVals = coordinates[coordinates['Tilt'] < 5.0].index #drops the chi values less than 5.0 as prescribed by the document
    coordinates.drop(DroppedVals , inplace=False)                      
    return name, coordinates
print(SpiralScheme("Austenite",90,6480)) # a test to see if the values are reasonable
