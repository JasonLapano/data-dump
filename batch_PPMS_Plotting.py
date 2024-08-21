# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 14:12:10 2019

@author: 5ul
"""

# -*- coding: utf-8 -*-
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


##Filenames######
#Format## [name(0), filelocation(1), filename(2), field(3), Rxx column(4), Ryy column(5), Rxy column(6), prefactor for Ryy/Rxx if needed(7) (leave as 1 if Rxx and Ryy are both measured), slope (should be -1)(8), measurement direction(9)]

filenames = [['1031_1', r'C:\Users\5ul\Google Drive\Pd-Bi-Se\PPMS\20191031_1_14T_2K', r'\20191031_RvH_14T.dat', 14, 19, 20, 21, 1, 1, 1],
             ['1101_1', r'C:\Users\5ul\Google Drive\Pd-Bi-Se\PPMS\20191101_1_20191101_2\20191101_1', r'\RvH_14T.dat', 14, 19, 19, 20, 1, -1, 1],
             ['1101_2', r'C:\Users\5ul\Google Drive\Pd-Bi-Se\PPMS\20191101_1_20191101_2\20191101_2', r'\RvH_14T.dat', 14, 21, 21, 22, 1, -1, 1],
             ['1101_3', r'C:\Users\5ul\Google Drive\Pd-Bi-Se\PPMS\20191101_3_20191107_1_Rvst_RvsH_14T\20191101_3', r'\RvH_14T.dat', 14, 19, 19, 20, 2.2 , -1, 1],
             ['1107_1', r'C:\Users\5ul\Google Drive\Pd-Bi-Se\PPMS\20191101_3_20191107_1_Rvst_RvsH_14T\20191107_1', r'\RvH_14T.dat', 14, 21, 21, 22, .92 , -1, 1],
             ['1004_4', r'C:\Users\5ul\Google Drive\Pd-Bi-Se\PPMS\20191004_4_Pd3Bi2Se2_350C', r'\20191004_4_Pd3Bi2Se2_RvsH_b.dat', 9, 19, 21, 20, 1, -1, -1],
             ['1010_2', r'C:\Users\5ul\Google Drive\Pd-Bi-Se\PPMS\20191009_3_20191010_2', r'\RvsH_2K.dat', 9, 21, 21, 22, 1, -1, -1],
             ['1108_1', r'C:\Users\5ul\Google Drive\Pd-Bi-Se\PPMS\20191108_1_20191108_2\20191108_1', r'\RvsH_9T.dat', 9, 19, 19, 20, .746, -1, -1], 
             ['1108_2', r'C:\Users\5ul\Google Drive\Pd-Bi-Se\PPMS\20191108_1_20191108_2\20191108_2', r'\RvsH_9T.dat', 9, 21, 21, 22, 1.26, -1, -1]];
             
for n in range(len(filenames)):
    filename = filenames[n][1]+filenames[n][2]
    data = np.genfromtxt(filename, # Specify filename
                     delimiter=',', # How is the data delimited?
                     skip_header=32, # Skip the first 32 rows
                     usecols=(3,4, filenames[n][4],filenames[n][5], filenames[n][6]), # Only read these columns
                     comments='Measurement') # Ignore bad lines
    data[:,3]=data[:,2]*filenames[n][7]
    direction = filenames[n][9]
    Bmax = filenames[n][3]+.01 #max field used in measurement
    geo = 4.53236 
    thickness= 1 #thickness in centimeters Van der Pauw: pi/ln(2) =4.53236; Hall bar = lengh/Area
    Rxy_np =filenames[n][8]#sign of Rxy (increaseing with increase B or decreasing with increasing B)
    
    # Determine where the magnetic field changes by more than 1000 Oe
    change = np.abs(np.diff(data[:,0])) > 2

    # Split the data where these changes occur
    split_data = np.split(data, np.where(change)[0]+1)
    
    #define dictionary of temperatures at which B is scanned, and plots the raw data
    templist={}
    for q in range(len(split_data)):
        templist[q]=round(split_data[q][0][0],1)
        dataplot = split_data[q]
    #    plot raw Rxy,Rxx,Ryy to inspect the data
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
############################################################## #since data goes from 0->+B,+B->-B, then -B->0; clips off the first and third legs
    split_data_field ={}
    changeBlist =[]
    for i in range(len(split_data)):   
        datacenter = split_data[i]
        delta = np.diff(datacenter[:,1])*direction#finds difference between i and i+1
        datasign = np.sign(delta)#finds the sign for each element <0 = -1 and >0 =1
        changeB =np.where(((np.roll(datasign, 1) - datasign) != 0).astype(int)==1)#creates a list where these is a zero
        if len(changeB[0])==2:#this is just in case there are more than 2 data points in changeB 
            cutmin =changeB[0][0]
            cutmax = changeB[0][1]
        else:
            cutmin =changeB[0][1]#this assumes that only 4 are maximally present
            cutmax = changeB[0][2]
        split_data_field[i]=np.split(datacenter,[cutmin,cutmax])[1]
##############################################################
##############################################################
#sym the data Rxx and Ryy are even[feven = (f(x)+f(-x))/2] Rxy is odd [fodd = (f(x)-f(-x))/2]
    split_data_sym={}
    final_data = {}
    final_datacond = {}
    carrier_fit = {}
    halldata =np.empty([len(split_data_field),5])#will store the resistance,high, medium, and low field slopes of Rxy
    for i in range(len(split_data_field)):
        T = split_data_field[i][:,0]
        Bp = split_data_field[i][:,1]/10000
        Bp[0]=Bmax*direction#makes the endpount = to Bmax to aviod an error for the endpoints out of bounds
        Bp[len(Bp)-1]=-1*Bmax*direction#makes the endpount = to Bmax to aviod an error for the endpoints out of bounds
        Rxx = split_data_field[i][:,2]
        Ryy = split_data_field[i][:,3]
        Rxy = split_data_field[i][:,4]
        Bm = Bp*-1
        Rxxint = interp1d(Bp, Rxx, kind='cubic')#sets up interpolation function
        Ryyint = interp1d(Bp, Ryy, kind='cubic')
        Rxyint = interp1d(Bp, Rxy, kind='cubic')
        Rxxm=Rxxint(Bm)#interpolates -B values
        Ryym=Ryyint(Bm)
        Rxym=Rxyint(Bm)
        Rxxsym = Rxxm/2+Rxx/2#even Rxx
        Ryysym = Ryym/2+Ryy/2#even Ryy
        Rsym = Rxxsym/2+Ryysym/2#average
        Rsqr = Rsym*geo#resistance per square
        Rxysym = Rxy_np*(Rxym/2-Rxy/2)#odd Rxy
        rho = Rsym*geo*thickness
        dGe2h=Rsqr/(Rsqr**2+Rxysym**2)/((1.6e-19)**2/6.626e-34)#conductance in units of e^2/h
        dGe2h_approx=1/Rsqr/((1.6e-19)**2/6.626e-34)#conductance in units of e^2/h
        #dGe2h=dGe2h_approx
    
        
        ##############derivative##################
        Rxysymint = interp1d(Bp, Rxysym, kind='cubic')#sets up interpolation function
        Bint=np.arange(-1*Bmax,Bmax,0.01)#creates a constant array
        Bint[0]=-Bmax#forces endpoints to match Bp's end points
        Bint[len(Bint)-1]=Bmax
        dRxydB=np.gradient(Rxysymint(Bint))*len(Bint)/2*1/Bmax#calculates derivative
        #note that the this is dRxy/dB = dB'/dB*dRxy/dB'; where B' = len(Rxy)/(2Bmax)*B+len(Rxy)/2
        dRxydBintoutput = interp1d(Bint, dRxydB, kind='cubic')#interpolates results to values of Bp
        dRxydBoutput= dRxydBintoutput(Bp)
    
        Rsqrsymint = interp1d(Bp, Rsqr, kind='cubic')#sets up interpolation function
        dRsqrdB=np.gradient(Rsqrsymint(Bint))*len(Bint)/2*1/Bmax#calculates derivative
        #note that the this is dRxy/dB = dB'/dB*dRxy/dB'; where B' = len(Rxy)/(2Bmax)*B+len(Rxy)/2
        dRsqrdBintoutput = interp1d(Bint, dRsqrdB, kind='cubic')#interpolates results to values of Bp
        dRsqrdBoutput= dRsqrdBintoutput(Bp)
        drhoBoutput=thickness*dRsqrdBoutput
        ################################
        
        ###########averages of slope of Rxy###########
        
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
        final_data[i]=np.stack((T,Bp,Rxysym,rho,dRxydBoutput,dRsqrdBoutput,drhoBoutput,Rsqr,Rsym,Rxy,Rxx,Ryy),axis=1)
        final_datacond[i]=np.stack((Bp,dGe2h),axis=1)
        
        #create file for fitting in wolfram
        carrier_fit[i] = np.stack((Bp, Rxysym, dRxydBoutput, Rsqr), axis=1)
        
        #plots the results to check that everything makes sense
        plt.figure(i+3)
        plt.subplot(131)
        plt.plot(Bp,Rxx)
        plt.plot(Bp,Rxxm)
        plt.plot(Bp,Rxxsym)
        plt.subplot(132)
        plt.plot(Bp,Ryy)
        plt.plot(Bp,Ryym)
        plt.plot(Bp,Ryysym)
        plt.subplot(133)
        plt.plot(Bp,Rxy)
        plt.plot(Bp,Rxym)
        plt.plot(Bp,Rxysym)
        plt.show()
#    #plot dG/(e^2/h)
#    plt.figure()
#    plt.plot(Bp,dGe2h)
#    plt.plot(Bp,dGe2h_approx)
#    plt.show()    
        
        plt.figure(i+4)
        plt.subplot(121)
        plt.plot(Bp,Rsqr)
        plt.subplot(122)
        plt.plot(Bp, Rxysym+n*.1)
        plt.legend(filenames[n][0])
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
    
    for i in range(len(carrier_fit)):
        # Put the temperature in the filename
        name = filename.strip('.dat')+'_carrier_fit_{:}K.csv'.format(templist[i])
        # Save the file
        np.savetxt(name, carrier_fit[i], delimiter=',')  
    