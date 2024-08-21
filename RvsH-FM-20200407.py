# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 12:00:47 2020

@author: 5ul
"""

"""
Created on Tue Sep 26 09:16:52 2017

@author: Matthew Brahlek_updated 9/26/2016
"""
import numpy as np
import matplotlib.pyplot as plt # plotting
from matplotlib import cm
import pandas as pd
from scipy.interpolate import interp1d
from scipy.misc import derivative
from scipy.interpolate import InterpolatedUnivariateSpline
#from scipy.misc import derivative
#plt.style.use('custom')

plt.close("all")


#################################
#################################
#######input information#########
#################################
#################################
filename = 'RvsHvsT.dat'#raw data from PPMS
direction = -1 #-1 takes +B->-B, +1 takes -B->+B
Bmax = 9.001 #max field used; important for interpolation fucntion
thickness = 1#thickness in centimeters
geo = 4.53236 #thickness in centimeters Van der Pauw: pi/ln(2) =4.53236; Hall bar = lengh/Area
Rxy_np =-1#sign of Rxy (increaseing with increase B or decreasing with increasing B)
#################################
#################################
#################################
#################################
#################################

data = np.genfromtxt(filename, # Specify filename
                     delimiter=',', # How is the data delimited?
                     skip_header=32, # Skip the first 32 rows
                     usecols=(3,4,19,20,21), # Only read these columns
                     comments='Measurement') # Ignore bad lines

# Determine where the magnetic field changes by more than 1000 Oe
change = np.abs(np.diff(data[:,0])) > 2

# Split the data where these changes occur
split_data = np.split(data, np.where(change)[0]+1)

#define dictionary of temperatures at which B is scanned, and plots the raw data
templist={}
for i in range(len(split_data)):
    templist[i]=round(split_data[i][0][0],1)
    dataplot = split_data[i]
    #plot raw Rxy,Rxx,Ryy to inspect the data
#    plt.figure(1)
#    plt.subplot(131)
#    plt.plot(dataplot[:,1]/10000, dataplot[:,4],'o')
#    plt.subplot(132)
#    plt.plot(dataplot[:,1]/10000, dataplot[:,3],'o')
#    plt.subplot(133)
#    plt.plot(dataplot[:,1]/10000, dataplot[:,2],'o')
#plt.show()
###############################################################  
###############################################################
plotbeg = 0
plotfin = len(split_data)

for j in range(3):   
    plt.figure(j)
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0,box.y0,box.width*0.8,box.height])
    for i in range(plotbeg,plotfin):
        name = '{:}K'.format(templist[i])
        ax.plot(split_data[i][:,1], split_data[i][:,2+j], label=name[:],marker='o',markersize=4)
            
    ax.legend(loc = 'center left',bbox_to_anchor=(1, 0.5))
    ax.set_xlabel('B(T)')
    ax.set_ylabel('R')
    plt.show()
##############################################################
############################################################## 

#since data goes from 0->+B,+B->-B, then -B->0; clips off the first and third legs
split_data_field ={}
upsweep = {}
downsweep = {}
changeBlist =[]
for i in range(len(split_data)):   
    datacenter = split_data[i]
    delta = np.diff(datacenter[:,1])*direction#finds difference between i and i+1
    datasign = np.sign(delta)#finds the sign for each element <0 = -1 and >0 =1
    changeB =np.where(((np.roll(datasign, 1) - datasign) != 0).astype(int)==1)#creates a list where these is a zero
    diff_changeB = np.diff(changeB)
    temp_changeB = [[]]
    temp_changeB[0].append(changeB[0][0])
    for k in range(len(diff_changeB[0])):
        if diff_changeB[0][k] > 10:  
            temp_changeB[0].append(changeB[0][k+1])
     
    temp_changeB[0].append(len(split_data[i]))             
    changeB = temp_changeB
        
    upsweep[i] = np.split(datacenter,[changeB[0][0],changeB[0][1]])[1]
    downsweep[i] = np.split(datacenter,[changeB[0][1],changeB[0][2]])[1]
##############################################################
##############################################################
    
#sym the data Rxx and Ryy are even[feven = (f(x)+f(-x))/2] Rxy is odd [fodd = (f(x)-f(-x))/2]
split_data_sym={}
final_data = {}
final_datacond = {}
halldata =np.empty([len(upsweep),5])#will store the resistance,high, medium, and low field slopes of Rxy
for i in range(len(upsweep)):
    T = upsweep[i][:,0]
    Bup = upsweep[i][:,1]/10000
    Bdown = downsweep[i][:,1]/10000
    Bup[0]=Bmax*direction#makes the endpount = to Bmax to aviod an error for the endpoints out of bounds
    Bdown[0]=Bmax*direction*-1#makes the endpount = to Bmax to aviod an error for the endpoints out of bounds
    Bup[len(Bup)-1]=-1*Bmax*direction#makes the endpount = to Bmax to aviod an error for the endpoints out of bounds
    Bdown[len(Bdown)-1]=Bmax*direction#makes the endpount = to Bmax to aviod an error for the endpoints out of bounds
    
    Rxx_up = upsweep[i][:,2]#Grabs resistance data on upsweep 
    Ryy_up = upsweep[i][:,3]
    Rxy_up = upsweep[i][:,4]
    
    #Will need to interpolate all of these so both up and down sweeps us the upsweep B field
    Rxx_down_pre = downsweep[i][:,2]#Grabs resistance data on downsweep
    Ryy_down_pre = downsweep[i][:,3]
    Rxy_down_pre = downsweep[i][:,4]
    
    Bupm = Bup*-1 #creates B-fields opposite sweep
    Bdownm = Bdown*-1
    
    Rxxint_up = interp1d(Bup, Rxx_up, kind='cubic')#sets up interpolation functions for upsweep
    Ryyint_up = interp1d(Bup, Ryy_up, kind='cubic')
    Rxyint_up = interp1d(Bup, Rxy_up, kind='cubic')
    
    Rxxint_down = interp1d(Bdown, Rxx_down_pre, kind='cubic')#sets up interpolation function for downsweep
    Ryyint_down = interp1d(Bdown, Ryy_down_pre, kind='cubic')
    Rxyint_down = interp1d(Bdown, Rxy_down_pre, kind='cubic')

##For symmetrization, interpolation is made so f(Rxx_upsweep) == f(-Rxx_downsweep), and f(Rxx_downsweep) == f(-Rxx_upsweep).  Always use Bup when plotting    
    Rxx_upm=Rxxint_up(Bupm)#interpolates -B values.
    Ryy_upm=Ryyint_up(Bupm)
    Rxy_upm=Rxyint_up(Bupm)
    
    Rxx_down=Rxxint_down(Bup) #interpolates -B values.  
    Ryy_down=Ryyint_down(Bup)
    Rxy_down=Rxyint_down(Bup)
    
    Rxx_downm=Rxxint_down(Bupm)#interpolates -B values.  
    Ryy_downm=Ryyint_down(Bupm)
    Rxy_downm=Rxyint_down(Bupm)
    
    Rxxsym_up = Rxx_downm/2+Rxx_up/2#even Rxx
    Rxxsym_down = Rxx_upm/2+Rxx_down/2#even Rxx
    Ryysym_up = Ryy_downm/2+Ryy_up/2#even Rxx
    Ryysym_down = Ryy_upm/2+Ryy_down/2#even Rxx
    Rsym_up = Rxxsym_up/2+Ryysym_up/2#average  
    Rsym_down = Rxxsym_down/2+Ryysym_down/2#average
    Rxysym_up = (-Rxy_downm/2 + Rxy_up/2)#odd Rxy
    Rxysym_down = (-Rxy_upm/2 + Rxy_down/2)#odd Rxy
    
##Combines up and down symmetrized scans.  B will be in the form of Bup + Bupm
    B = np.append(Bup, np.flip(Bup))
    Rxx = np.append(Rxx_up, np.flipud(Rxx_down))
    Rxxm = np.append(Rxx_upm, np.flipud(Rxx_downm))
    Ryy = np.append(Ryy_up, np.flipud(Ryy_down))
    Ryym = np.append(Ryy_upm, np.flipud(Ryy_downm))
    Rxy = np.append(Rxy_up, np.flipud(Rxy_down))
    Rxym = np.append(np.flipud(Rxy_downm), Rxy_upm)
    Rxxsym = np.append(Rxxsym_up, np.flipud(Rxxsym_down))
    Ryysym = np.append(Ryysym_up, np.flipud(Ryysym_down))
    Rxysym = np.append(Rxysym_up, np.flipud(Rxysym_down))
    plt.figure(4)
    plt.plot(Bup, Rxxsym_up, marker='.')
    plt.plot(Bup, Rxxsym_down, marker='.')
    plt.figure(5)
    plt.plot(Bup, Ryysym_up, marker='.')
    plt.plot(Bup, Ryysym_down, marker='.')
    plt.figure(6)
    plt.plot(Bup, Rxysym_up, marker='.')
    plt.plot(Bup, Rxysym_down, marker='.')
    
##############ALL THIS STUFF NEEDS REDONE################################################
    
    Rsym = Rxxsym/2+Ryysym/2#average
    Rsqr = Rsym*geo#resistance per square
    rho = Rsym*geo*thickness
    dGe2h=Rsqr/(Rsqr**2+Rxysym**2)/((1.6e-19)**2/6.626e-34)#conductance in units of e^2/h
    dGe2h_approx=1/Rsqr/((1.6e-19)**2/6.626e-34)#conductance in units of e^2/h
    #dGe2h=dGe2h_approx

    
    ##############derivative##################
    Rxysymint_up = interp1d(Bup, Rxysym_up, kind='cubic')#sets up interpolation function
    Bint_up=np.arange(-1*Bmax,Bmax,0.01)#creates a constant array
    Bint_up[0]=-Bmax#forces endpoints to match Bp's end points
    Bint_up[len(Bint_up)-1]=Bmax
    dRxydB_up=np.gradient(Rxysymint_up(Bint_up))*len(Bint_up)/2*1/Bmax#calculates derivative
    #note that the this is dRxy/dB = dB'/dB*dRxy/dB'; where B' = len(Rxy)/(2Bmax)*B+len(Rxy)/2
    dRxydBintoutput_up = interp1d(Bint_up, dRxydB_up, kind='cubic')#interpolates results to values of Bp
    dRxydBoutput_up= dRxydBintoutput_up(Bup)

    Rsqrsymint_up = interp1d(Bup, (Rxxsym_up/2+Ryysym_up/2)*geo, kind='cubic')#sets up interpolation function
    dRsqrdB_up=np.gradient(Rsqrsymint_up(Bint_up))*len(Bint_up)/2*1/Bmax#calculates derivative
    #note that the this is dRxy/dB = dB'/dB*dRxy/dB'; where B' = len(Rxy)/(2Bmax)*B+len(Rxy)/2
    dRsqrdBintoutput_up = interp1d(Bint_up, dRsqrdB_up, kind='cubic')#interpolates results to values of Bp
    dRsqrdBoutput_up= dRsqrdBintoutput_up(Bup)
    drhoBoutput_up=thickness*dRsqrdBoutput_up
    
    ############Derivative for down slope########
    Rxysymint_down = interp1d(Bup, Rxysym_down, kind='cubic')#sets up interpolation function
    Bint_down=np.arange(-1*Bmax,Bmax,0.01)#creates a constant array
    Bint_down[0]=-Bmax#forces endpoints to match Bp's end points
    Bint_down[len(Bint_down)-1]=Bmax
    dRxydB_down=np.gradient(Rxysymint_down(Bint_down))*len(Bint_down)/2*1/Bmax#calculates derivative
    #note that the this is dRxy/dB = dB'/dB*dRxy/dB'; where B' = len(Rxy)/(2Bmax)*B+len(Rxy)/2
    dRxydBintoutput_down = interp1d(Bint_down, dRxydB_down, kind='cubic')#interpolates results to values of Bp
    dRxydBoutput_down= dRxydBintoutput_down(Bup)

    Rsqrsymint_down = interp1d(Bup, (Rxxsym_down/2+Ryysym_down/2)*geo, kind='cubic')#sets up interpolation function
    dRsqrdB_down=np.gradient(Rsqrsymint_down(Bint_down))*len(Bint_down)/2*1/Bmax#calculates derivative
    #note that the this is dRxy/dB = dB'/dB*dRxy/dB'; where B' = len(Rxy)/(2Bmax)*B+len(Rxy)/2
    dRsqrdBintoutput_down = interp1d(Bint_down, dRsqrdB_down, kind='cubic')#interpolates results to values of Bp
    dRsqrdBoutput_down= dRsqrdBintoutput_down(Bup)
    drhoBoutput_down=thickness*dRsqrdBoutput_down
    
    ##############Combines Them###########################
    dRsqrdBoutput = np.append(dRsqrdBoutput_up, np.flipud(dRsqrdBoutput_down))
    
    dRxydBoutput = np.append(dRxydBoutput_up, np.flipud(dRxydBoutput_down))
    plt.figure(7)
    plt.plot(B, dRxydBoutput, marker='.')

    drhoBoutput = np.append(drhoBoutput_up, np.flipud(drhoBoutput_down))
    plt.figure(8)
    plt.plot(B, drhoBoutput, marker='.')    
    ################################
    
    ###########averages of slope of Rxy###########
    
    ###Just averages slopes not really right but worth trying
    dRxydB = dRxydB_up/2 + dRxydB_down/2    
    skip= 20
    span = 40
    quarter = round(len(dRxydB)/4)
    mid = round(len(dRxydB)/2)
    halldata[i,0]=templist[i]#store the temperature
    halldata[i,1]=Rsqr[int(round(len(Rsqr)/2,1))]#store the temperature
    
    mean = 0
    for x in range(skip,skip+span,1):#averages the high-field portion
        mean = mean +dRxydB[x]
    halldata[i,2]= mean/span#since 100-20=80
    
    mean = 0
    for x in range(int(quarter-span/2),int(quarter+span/2),1):#averages the mid-field portion
        mean = mean +dRxydB[x]
    halldata[i,3]= mean/span#since 100-20=80
    
    mean = 0    
    for x in range(int(mid-span/2),int(mid+span/2),1):#averages the 
        mean = mean +dRxydB[x]
    halldata[i,4]= mean/span#since 100-20=80        
    ##############################################
    
    ##################create plots to check results####################
    Tadj = np.append(T, np.flipud(T))    
    final_data[i]=np.stack((Tadj,B,Rxysym,rho,dRxydBoutput,dRsqrdBoutput,drhoBoutput,Rsqr,Rsym,Rxy,Rxx,Ryy),axis=1)
    final_datacond[i]=np.stack((B,dGe2h),axis=1)
    #plots the results to check that everything makes sense
    plt.figure(i+3)
    plt.subplot(131)
    plt.plot(B,Rxx)
    plt.plot(B,Rxxm)
    plt.plot(B,Rxxsym)
    plt.subplot(132)
    plt.plot(B,Ryy)
    plt.plot(B,Ryym)
    plt.plot(B,Ryysym)
    plt.subplot(133)
    plt.plot(B,Rxy)
    plt.plot(B,Rxym)
    plt.plot(B,Rxysym)
    plt.show()
#plot dG/(e^2/h)
plt.figure()
plt.plot(B,dGe2h)
plt.plot(B,dGe2h_approx)
plt.show()    
##############################################################
# Export each piece of the split data
head = 'T(K),B(T),Rxy(ohms),rho(ohm*cm),dRxy/dB,dRsqr/dB,drho/dB,Rave(ohm/sq),Rsym(ohm),Rxy_raw,Rxx_raw,Ryy_raw,'
for i in range(len(final_data)):
    # Put the temperature in the filename
    name = filename.strip('.dat')+'_output_{:}K.csv'.format(templist[i])
    # Save the file
    np.savetxt(name, final_data[i], delimiter=',', header=head)  

head = 'B(T),dGe2h'
for i in range(len(final_data)):
    # Put the temperature in the filename
    name = filename.strip('.dat')+'_outputconduc_{:}K.csv'.format(templist[i])
    # Save the file
    np.savetxt(name, final_datacond[i], delimiter=',')  


head = 'T,R(ohm/sqr),dRxydB(highfield),dRxydB(midfield),dRxydB(lowfield)'

# Save the file
name = filename.strip('.dat')+'_Hall data together.csv'
np.savetxt(name, halldata, delimiter=',', header=head)       
           

