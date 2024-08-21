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

#head = 'T(K),B(T),1/B, Rxy(ohms),Rxy-background(ohm), rho(ohm*cm), Magneto,dRxy/dB,dRsqr/dB,drho/dB,Rave(ohm/sq),Rsym(ohm),Rxy_raw,Rxx_raw,Ryy_raw, Delta-Rxx, Delta-Rxy, Rxx_ave, Rxy_ave'
# S:
#     0 = T
#     1 = B
#     2 = 1/B
#     3 = Rxy
#     4 = Rxy-AH
#     5 = rho
#     6 = MR
#     7 = dRxy/dB
#     8 = dRsquare/dB
#     9 = drho/dB
#     10 = Rave
#     11 = Rsym
#     12 = Rxy_raw
#     13 = Rxx_raw
#     14 = Ryy_raw
#     15 = delta_Rxx
#     16 = delta_Rxy
#     17 = Rxx_ave
#     18 = Rxy_ave


def RvsH_FM_import(filename, cols=[3,4,19,20,21], head=32, geo = 4.532, deli=',', Bmax=9):
    plt.close("all")


    direction = -1 #-1 takes +B->-B, +1 takes -B->+B
    # Bmax = 14 #max field used; important for interpolation fucntion
    thickness = 1#thickness in centimeters
    Rxy_np =-1#sign of Rxy (increaseing with increase B or decreasing with increasing B)
    saturation = 10 # % of data used to fit the linear background (i.e. for a 10T scan, 10 would give a fit over 9-10T data range)
    
    #################################
    #################################
    #################################
    #################################
    #################################
    
    scans = 1 #  Set to 1 if you are running a series of temperature scans
    
    #################################
    #################################
    #################################
    #################################
    #################################
    

    data = np.genfromtxt(filename, # Specify filename
                         delimiter=deli, # How is the data delimited, ',' for comma, '\t' for tab?
                         skip_header=head, # Skip the first 32 rows
                         usecols=cols, # Only read these columns
                         comments='Time') # Ignore bad lines
    
   
    # Determine where the magnetic field changes by more than 1000 Oe
    change = np.abs(np.diff(data[:,0])) > 1
    
    # Split the data where these changes occur
    split_data = np.split(data, np.where(change)[0]+1)
    
    #define dictionary of temperatures at which B is scanned, and plots the raw data
    templist={}
    for i in range(len(split_data)):
        templist[i]=round(split_data[i][0][0],1)
        dataplot = split_data[i]
    #    print(i)
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
        # plt.show()
    ##############################################################
    ############################################################## 
    
    #since data goes from 0->+B,+B->-B, then -B->0; clips off the first and third legs
    split_data_field ={}
    upsweep = {}
    downsweep = {}
    changeBlist =[]
    for i in range(len(split_data)):   
        datacentertest = split_data[i]
        #####Remove Duplicate Numbers in Array#####
        rounds = np.round(datacentertest[:, 1],-2)  
        nonduplicates = np.nonzero(np.diff(rounds))
        datacentertest = datacentertest[nonduplicates, :]
        datacenter = datacentertest[0]
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
         
        if 3 == np.size(changeB): 
            upsweep[i] = np.split(datacenter,[changeB[0][0],changeB[0][1]])[1]
            downsweep[i] = np.split(datacenter,[changeB[0][1],changeB[0][2]])[1]
            
        if 5 == np.size(changeB): 
            upsweep[i] = np.split(datacenter,[changeB[0][1],changeB[0][2]])[1]
            downsweep[i] = np.split(datacenter,[changeB[0][2],changeB[0][3]])[1]
                
    
    ##############################################################
    ##############################################################
        
    #sym the data Rxx and Ryy are even[feven = (f(x)+f(-x))/2] Rxy is odd [fodd = (f(x)-f(-x))/2]
    split_data_sym={}
    final_data = {}
    final_datacond = {}
    halldata =np.empty([len(upsweep),5])#will store the resistance,high, medium, and low field slopes of Rxy
    for i in range(len(upsweep)):
        print(i)
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
        B = np.append(Bup, np.flipud(Bup))
        Rxx = np.append(Rxx_up, np.flipud(Rxx_down))
        Rxxm = np.append(Rxx_upm, np.flipud(Rxx_downm))
        Ryy = np.append(Ryy_up, np.flipud(Ryy_down))
        Ryym = np.append(Ryy_upm, np.flipud(Ryy_downm))
        Rxy = np.append(Rxy_up, np.flipud(Rxy_down))
        Rxym = np.append(np.flipud(Rxy_downm), Rxy_upm)
        Rxxsym = np.append(Rxxsym_up, np.flipud(Rxxsym_down))
        Ryysym = np.append(Ryysym_up, np.flipud(Ryysym_down))
        Rxysym = np.append(Rxysym_up, np.flipud(Rxysym_down))
        # plt.figure(4)
        # plt.plot(Bup, Rxxsym_up, marker='.')
        # plt.plot(Bup, Rxxsym_down, marker='.')
        # plt.figure(5)
        # plt.plot(Bup, Ryysym_up, marker='.')
        # plt.plot(Bup, Ryysym_down, marker='.')
        # plt.figure(6)
        # plt.plot(Bup, Rxysym_up, marker='.')
        # plt.plot(Bup, Rxysym_down, marker='.')
        
    ##############ALL THIS STUFF NEEDS REDONE################################################
        
        Rsym = Rxxsym/2+Ryysym/2#average
        Rsqr = Rsym*geo#resistance per square
        rho = Rsym*geo*thickness
        R0 = Rsym[np.argmin(abs(Bup))]*geo
        rho_0 = Rsym[np.argmin(abs(Bup))]*geo*thickness
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
        # plt.figure(7)
        # plt.plot(B, dRxydBoutput, marker='.')
    
        drhoBoutput = np.append(drhoBoutput_up, np.flipud(drhoBoutput_down))
        # plt.figure(8)
        # plt.plot(B, drhoBoutput, marker='.')    
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
        
        #######Subtract Linear Background and calculate N2D###############
        
        linear = np.where(B>np.amax(B)*(1-saturation*.01)) #finds where data is saturated for linear fit
        linear = linear[0] #converts from tuple to float
        #substract linear background and plot
        z = np.polyfit(B[linear], Rxysym[linear], 1);
        Rxy_sub = Rxysym - B*z[0] 
        Rxy_subup = Rxysym_up - Bup*z[0] 
        Rxy_subdown = Rxysym_down-Bup*z[0]
        n2d = 1/(z[0] * 1.602*10**-19*10000)
        # plt.figure(9)
        # plt.plot(B, Rxy_sub, marker='.')    
        
         #######Plot dxx and dxy in 1/B###############
        # plt.figure(10)
        # plt.plot(1/B, drhoBoutput)
        # plt.figure(11)
        # plt.plot(1/B, dRxydBoutput)
        
        
        ############# %Magnetoresistance ###################
        mag_res = ((Rsym*geo-R0)/R0)*100
        inv_B = 1/B 
        
        ################ DIFFERENCE IN FIELD SWEEPS ################
        Rxx_diff = Rxxsym_up-Rxxsym_down
        Rxy_diff = Rxysym_up-Rxysym_down
        
 
        
        
        fRxx_diff = np.append(Rxx_diff, np.flipud(Rxx_diff))
        fRxy_diff = np.append(Rxy_diff, np.flipud(Rxy_diff))
        
        
        ################ Averaged up and Down #########################
        
        Rxx_ave = (Rxxsym_up+Rxxsym_down)/2
        Rxy_ave = (Rxysym_up+Rxysym_down)/2
        
        fRxx_ave = np.append(Rxx_ave, np.flipud(Rxx_ave))
        fRxy_ave = np.append(Rxy_ave, np.flipud(Rxy_ave))
        
        
        ##################create plots to check results####################
        Tadj = np.append(T, np.flipud(T))    
        final_data[i]=np.stack((Tadj,B,inv_B, Rxysym,Rxy_sub, rho, mag_res,dRxydBoutput,dRsqrdBoutput,drhoBoutput,Rsqr,Rsym,Rxy,Rxx,Ryy, fRxx_diff, fRxy_diff, fRxx_ave, fRxy_ave),axis=1)
        final_datacond[i]=np.stack((B,dGe2h),axis=1)
        #plots the results to check that everything makes sense
        # plt.figure(i+3)
        # plt.subplot(131)
        # plt.plot(B,Rxx)
        # plt.plot(B,Rxxm)
        # plt.plot(B,Rxxsym)
        # plt.subplot(132)
        # plt.plot(B,Ryy)
        # plt.plot(B,Ryym)
        # plt.plot(B,Ryysym)
        # plt.subplot(133)
        # plt.plot(B,Rxy)
        # plt.plot(B,Rxym)
        # plt.plot(B,Rxysym)
        # plt.show()
        
        #plot dG/(e^2/h)
        # plt.figure()
        # plt.plot(B,dGe2h)
        # plt.plot(B,dGe2h_approx)
        # plt.show()   
    
        plt.figure(2)
        plt.clf()
        plt.rcParams.update({'font.size': 6})
        ax1 = plt.subplot(431)
        plt.plot(Bup, Rxx_up)
        plt.plot(Bup, Rxx_down)
        ax1.set_title("Raw Rxx")
        ax1.set_xlabel("Field(B)")
        ax1.set_ylabel("Resistance(Ohm)")
        
        ax2 = plt.subplot(432)
        plt.plot(Bup, Ryy_up)
        plt.plot(Bup, Ryy_down)
        ax2.set_title("Raw Ryy")
        ax2.set_xlabel("Field(B)")
        ax2.set_ylabel("Resistance(Ohm)")
        
        ax3 = plt.subplot(433)
        plt.plot(Bup, Rxy_up)
        plt.plot(Bup, Rxy_down)
        ax3.set_title("Raw Hall")
        ax3.set_xlabel("Field(B)")
        ax3.set_ylabel("Resistance(Ohm)")
        
        ax4 = plt.subplot(434)
        plt.plot(Bup, Rsym_up*geo)
        plt.plot(Bup, Rsym_down*geo)
        ax4.set_title("R Symmeterized")
        ax4.set_xlabel("Field(B)")
        ax4.set_ylabel("Resistance(Ohm)")
        
    
        ax5 = plt.subplot(435)
        plt.plot(Bup, ((Rsym_up*geo-R0)/R0)*100)
        plt.plot(Bup, ((Rsym_down*geo-R0)/R0)*100)
        ax5.set_title("% Magnetoresistance")
        ax5.set_xlabel("Field(B)")
        ax5.set_ylabel("(R-R0)/R0 %")
        
        ax6 = plt.subplot(436)
        plt.plot(Bup[np.where(np.absolute(Bup)<9)], (Rsym_up[np.where(np.absolute(Bup)<9)]*geo-R0)/R0)
        plt.plot(Bup[np.where(np.absolute(Bup)<9)], (Rsym_down[np.where(np.absolute(Bup)<9)]*geo-R0)/R0)
        ax6.set_title("R Symmeterized")
        ax6.set_xlabel("Field(B)")
        ax6.set_ylabel("Resistance(Ohm)")
        ax6.set_xlim([-2, 2])
        plt.show
        
        
        
        
        ax7 = plt.subplot(437)
        plt.plot(Bup, (Rxysym_up))
        plt.plot(Bup, (Rxysym_down))
        ax7.set_title("Symmetrized Hall")
        ax7.set_xlabel("Field(B)")
        ax7.set_ylabel("R-hall")
        
        ax8 = plt.subplot(438)
        plt.plot(Bup[np.where(np.absolute(Bup)<9)], Rxy_subup[np.where(np.absolute(Bup)<9)])
        plt.plot(Bup[np.where(np.absolute(Bup)<9)], Rxy_subdown[np.where(np.absolute(Bup)<9)])
        ax8.set_title("Anamolous Hall")
        ax8.set_xlabel("Field(B)")
        ax8.set_ylabel("R-hall")
        
        
        ax9 = plt.subplot(439)
        s1 = 'R(0) = ' + str(np.round(R0, decimals=1)) + ' Ohm'
        plt.text(.05, .8, s1)
        s2 = 'High Field Slope = '+ str(np.round(z[0], decimals=1))
        plt.text(.05, .6, s2)
        s3 = 'N 2D = ' + "{:e}".format(n2d)
        plt.text(.05, .4, s3)
        s4 = 'Mobility = ' + str(abs(np.round((1/(rho_0/thickness))/n2d/(1.602*10**-19), decimals=1))) + ' cm2/Vs'
        plt.text(.05, .2, s4)
        ax9.get_xaxis().set_visible(False)
        ax9.get_yaxis().set_visible(False)
        
        ax10 = plt.subplot(4,3,10)
        Babs = np.absolute(B)
        plt.plot(1/Bup[np.where(((Bup)>.5) & ((Bup)<(Bmax-0.1)))], drhoBoutput_up[np.where(((Bup)>.5) & ((Bup)< (Bmax-0.1)))])
        plt.plot(1/Bup[np.where(((Bup)>.5) & ((Bup)<(Bmax-0.1)))], drhoBoutput_down[np.where(((Bup)>.5) & ((Bup)< (Bmax-0.1)))])
        ax10.set_xlim(left=0)
        ax10.set_title("dRxx/dB")
        ax10.set_ylabel("dRxx/Db")
        ax10.set_xlabel("1/B")
        
        ax11 = plt.subplot(4,3,11)
        plt.plot(Bup, Rxx_diff)
        ax11.set_title("Delta Magnetoresistance")
        ax11.set_ylabel("Delta Ohms")
        ax11.set_xlabel("B")
        
        ax11 = plt.subplot(4,3,12)
        plt.plot(Bup, Rxy_diff)
        ax11.set_title("Delta Hall Effect")
        ax11.set_ylabel("Delta Ohms")
        ax11.set_xlabel("B")
        
        
        plt.tight_layout(pad=.0, h_pad=.0, w_pad=None, rect=None)
        plt.rcParams.update({'font.size': 30})
        figname=filename.strip('.dat')+'_{:}K'.format(templist[i])
        plt.savefig(figname+'.png', dpi=400)
        plt.savefig(figname+'.svg')
        
        
    
    
    ################# COMPILED PLOT FOR ALL SCANS ####################################
        
        if scans == 1:       #only plots if you are doing a tempscan.  This is done stupidly and just adds the next plot over the last one. 
                             #This is because I didn't use a counter (i) to plot them, so it should be redone but works for nwo  
            
           
            
            plt.figure(3)    
            plt.rcParams.update({'font.size': 6})
            ax1 = plt.subplot(431)
            plt.plot(Bup, Rxx_up, color=((i/len(templist)), 0, (1-(i/len(templist)))), label=templist[i])
            plt.plot(Bup, Rxx_down, color=((i/len(templist)), 0, (1-(i/len(templist)))))
            ax1.set_title("Raw Rxx")
            ax1.set_xlabel("Field(B)")
            ax1.set_ylabel("Resistance(Ohm)")
            plt.legend()
            
            ax2 = plt.subplot(432)
            plt.plot(Bup, Ryy_up, color=((i/len(templist)), 0, (1-(i/len(templist)))), label='')
            plt.plot(Bup, Ryy_down, color=((i/len(templist)), 0, (1-(i/len(templist)))))
            ax2.set_title("Raw Ryy")
            ax2.set_xlabel("Field(B)")
            ax2.set_ylabel("Resistance(Ohm)")
            
            ax3 = plt.subplot(433)
            plt.plot(Bup, Rxy_up, color=((i/len(templist)), 0, (1-(i/len(templist)))))
            plt.plot(Bup, Rxy_down, color=((i/len(templist)), 0, (1-(i/len(templist)))))
            ax3.set_title("Raw Hall")
            ax3.set_xlabel("Field(B)")
            ax3.set_ylabel("Resistance(Ohm)")
            
            ax4 = plt.subplot(434)
            plt.plot(Bup, Rsym_up*geo, color=((i/len(templist)), 0, (1-(i/len(templist)))))
            plt.plot(Bup, Rsym_down*geo, color=((i/len(templist)), 0, (1-(i/len(templist)))))
            ax4.set_title("R Symmeterized")
            ax4.set_xlabel("Field(B)")
            ax4.set_ylabel("Resistance(Ohm)")
            
            ax5 = plt.subplot(435)
            plt.plot(Bup, ((Rsym_up*geo-R0)/R0)*100, color=((i/len(templist)), 0, (1-(i/len(templist)))))
            plt.plot(Bup, ((Rsym_down*geo-R0)/R0)*100, color=((i/len(templist)), 0, (1-(i/len(templist)))))
            ax5.set_title("% Magnetoresistance")
            ax5.set_xlabel("Field(B)")
            ax5.set_ylabel("(R-R0)/R0 %")
            
            ax6 = plt.subplot(436)
            plt.plot(Bup[np.where(np.absolute(Bup)<2)], (Rsym_up[np.where(np.absolute(Bup)<2)]*geo-R0)/R0*100, color=((i/len(templist)), 0, (1-(i/len(templist)))))
            plt.plot(Bup[np.where(np.absolute(Bup)<2)], (Rsym_down[np.where(np.absolute(Bup)<2)]*geo-R0)/R0*100, color=((i/len(templist)), 0, (1-(i/len(templist)))))
            ax6.set_title("% Magnetoresistance")
            ax6.set_xlabel("Field(B)")
            ax6.set_ylabel("(R-R0)/R0 %")
            ax6.set_xlim([-2, 2])
            plt.show     
            
            
            ax7 = plt.subplot(437)
            plt.plot(Bup, (Rxysym_up), color=((i/len(templist)), 0, (1-(i/len(templist)))))
            plt.plot(Bup, (Rxysym_down), color=((i/len(templist)), 0, (1-(i/len(templist)))))
            ax7.set_title("Symmetrized Hall")
            ax7.set_xlabel("Field(B)")
            ax7.set_ylabel("R-hall")
            
            ax8 = plt.subplot(438)
            plt.plot(Bup[np.where(np.absolute(Bup)<9)], Rxy_subup[np.where(np.absolute(Bup)<9)], color=((i/len(templist)), 0, (1-(i/len(templist)))))
            plt.plot(Bup[np.where(np.absolute(Bup)<9)], Rxy_subdown[np.where(np.absolute(Bup)<9)], color=((i/len(templist)), 0, (1-(i/len(templist)))))
            ax8.set_title("Anamolous Hall")
            ax8.set_xlabel("Field(B)")
            ax8.set_ylabel("R-hall")
            
            
            # ax9 = plt.subplot(439)
            # s1 = 'R(0) = ' + str(np.round(R0, decimals=1)) + ' Ohm'
            # plt.text(.05, .8, s1)
            # s2 = 'High Field Slope = '+ str(np.round(z[0], decimals=1))
            # plt.text(.05, .6, s2)
            # s3 = 'N 2D = ' + "{:e}".format(n2d)
            # plt.text(.05, .4, s3)
            # s4 = 'Mobility = ' + str(abs(np.round((1/(rho_0/thickness))/n2d/(1.602*10**-19), decimals=1))) + ' cm2/Vs'
            # plt.text(.05, .2, s4)
            # ax9.get_xaxis().set_visible(False)
            # ax9.get_yaxis().set_visible(False)
            
            ax10 = plt.subplot(4,3,10)
            Babs = np.absolute(B)
            plt.plot(1/Bup[np.where(((Bup)>0.5) & ((Bup)<(Bmax-0.1)))], drhoBoutput_up[np.where(((Bup)>0.5) & ((Bup)< (Bmax-0.1)))], color=((i/len(templist)), 0, (1-(i/len(templist)))))
            plt.plot(1/Bup[np.where(((Bup)>0.5) & ((Bup)<(Bmax-0.1)))], drhoBoutput_down[np.where(((Bup)>0.5) & ((Bup)< (Bmax-0.1)))], color=((i/len(templist)), 0, (1-(i/len(templist)))))
            ax10.set_xlim(left=0)
            ax10.set_title("dRxx/dB")
            ax10.set_ylabel("dRxx/Db")
            ax10.set_xlabel("1/B")
            
            ax11 = plt.subplot(4,3,11)
            plt.plot(Bup, Rxx_diff, color=((i/len(templist)), 0, (1-(i/len(templist)))))
            ax11.set_title("Delta Magnetoresistance")
            ax11.set_ylabel("Delta Ohms")
            ax11.set_xlabel("B")
            
            ax11 = plt.subplot(4,3,12)
            plt.plot(Bup, Rxy_diff, color=((i/len(templist)), 0, (1-(i/len(templist)))))
            ax11.set_title("Delta Hall Effect")
            ax11.set_ylabel("Delta Ohms")
            ax11.set_xlabel("B")
        
        #
    
    plt.figure(3)
    plt.tight_layout(pad=.0, h_pad=.0, w_pad=None, rect=None)
    plt.rcParams.update({'font.size': 30})
    figname=filename.strip('.dat') + '_compiled'
    plt.savefig(figname+'.png', dpi = 400)
    plt.savefig(figname+'.svg')
        
        
        
        
    
    ##############################################################
    # Export each piece of the split data
    head = 'T(K),B(T),1/B, Rxy(ohms),Rxy-background(ohm), rho(ohm*cm), Magneto,dRxy/dB,dRsqr/dB,drho/dB,Rave(ohm/sq),Rsym(ohm),Rxy_raw,Rxx_raw,Ryy_raw, Delta-Rxx, Delta-Rxy, Rxx_ave, Rxy_ave'
    for i in range(len(final_data)):
        # Put the temperature in the filename
        name = filename.strip('.dat')+'_{:}K_RvsH.csv'.format(templist[i])
        # Save the file
        np.savetxt(name, final_data[i], delimiter=',', header=head)  
    
    head = 'B(T),dGe2h'
    for i in range(len(final_data)):
        # Put the temperature in the filename
        name = filename.strip('.dat')+'_{:}K_cond.csv'.format(templist[i])
        # Save the file
        np.savetxt(name, final_datacond[i], delimiter=',')  
    
    
    head = 'T,R(ohm/sqr),dRxydB(highfield),dRxydB(midfield),dRxydB(lowfield)'
    
    # Save the file
    name = filename.strip('.dat')+'_Hall data together.csv'
    np.savetxt(name, halldata, delimiter=',', header=head)       
               
    
    return final_data








def RvsH_PM_import(filename, cols=[3,4,19,20,21], head=32, geo = 4.532, deli=',', Bmax=9):
    plt.close("all")


    direction = -1 #-1 takes +B->-B, +1 takes -B->+B
    # Bmax = 14 #max field used; important for interpolation fucntion
    thickness = 1#thickness in centimeters
    Rxy_np =-1#sign of Rxy (increaseing with increase B or decreasing with increasing B)
    saturation = 10 # % of data used to fit the linear background (i.e. for a 10T scan, 10 would give a fit over 9-10T data range)
    
    #################################
    #################################
    #################################
    #################################
    #################################
    
    scans = 1 #  Set to 1 if you are running a series of temperature scans
    
    #################################
    #################################
    #################################
    #################################
    #################################
    

    data = np.genfromtxt(filename, # Specify filename
                         delimiter=deli, # How is the data delimited, ',' for comma, '\t' for tab?
                         skip_header=head, # Skip the first 32 rows
                         usecols=cols, # Only read these columns
                         comments='Time') # Ignore bad lines
    
    #Correcting for Rxx and Ryy
    data[:,3] =  data[:,3]
    data[:, 1] = 10000*data[:, 1]
    
                         
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
    changeBlist =[]
    for i in range(len(split_data)):   
        datacentertest = split_data[i]
        #####Remove Duplicate Numbers in Array#####
        rounds = datacentertest[:, 1]  
        nonduplicates = np.nonzero(np.diff(rounds))
        datacentertest = datacentertest[nonduplicates, :]
        datacenter = datacentertest[0]
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
    print(i)
    split_data_sym={}
    final_data = {}
    final_datacond = {}
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
        rho_0 = Rsym[np.argmin(abs(Bp))]*geo*thickness
        R0 = Rsym[np.argmin(abs(Bp))]*geo
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
        
        #######Subtract Linear Background and calculate N2D###############
        
        linear = np.where(Bp>np.amax(Bp)*(1-saturation*.01)) #finds where data is saturated for linear fit
        linear = linear[0] #converts from tuple to float
        #substract linear background and plot
        z = np.polyfit(Bp[linear], Rxysym[linear], 1);
        Rxy_sub = Rxysym - Bp*z[0] 
    #    Rxy = Rxysym - Bp*z[0] 
        n2d = 1/(z[0] * 1.602*10**-19*10000)
        # plt.figure(9)
        # plt.plot(B, Rxy_sub, marker='.')    
        
        
        ############# %Magnetoresistance ###################
        mag_res = ((Rsym*geo-R0)/R0)*100
        inv_B = 1/Bp
        
        ##################create plots to check results####################    
        final_data[i]=np.stack((T,Bp,inv_B,Rxysym,Rxy_sub,mag_res,rho,dRxydBoutput,dRsqrdBoutput,drhoBoutput,Rsqr,Rsym,Rxy,Rxx,Ryy),axis=1)
        final_datacond[i]=np.stack((Bp,dGe2h),axis=1)
        #plots the results to check that everything makes sense
        # plt.figure(i+3)
        # plt.subplot(331)
        # plt.plot(Bp,Rxx)
        # plt.plot(Bp,Rxxm)
        # plt.plot(Bp,Rxxsym)
        # plt.subplot(132)
        # plt.plot(Bp,Ryy)
        # plt.plot(Bp,Ryym)
        # plt.plot(Bp,Ryysym)
        # plt.subplot(133)
        # plt.plot(Bp,Rxy)
        # plt.plot(Bp,Rxym)
        # plt.plot(Bp,Rxysym)
        # plt.show()
        
        plt.figure(2)
        plt.clf()
        plt.rcParams.update({'font.size': 6})
        ax1 = plt.subplot(431)
        plt.plot(Bp, Rxx)
        ax1.set_title("Raw Rxx")
        ax1.set_xlabel("Field(B)")
        ax1.set_ylabel("Resistance(Ohm)")
        
        ax2 = plt.subplot(432)
        plt.plot(Bp, Ryy)
        ax2.set_title("Raw Ryy")
        ax2.set_xlabel("Field(B)")
        ax2.set_ylabel("Resistance(Ohm)")
        
        ax3 = plt.subplot(433)
        plt.plot(Bp, Rxy)
        ax3.set_title("Raw Hall")
        ax3.set_xlabel("Field(B)")
        ax3.set_ylabel("Resistance(Ohm)")
        
        ax4 = plt.subplot(434)
        plt.plot(Bp, Rsym*geo)
        ax4.set_title("R Symmeterized")
        ax4.set_xlabel("Field(B)")
        ax4.set_ylabel("Resistance(Ohm)")
        
    
        ax5 = plt.subplot(435)
        plt.plot(Bp, ((Rsym*geo-R0)/R0)*100)
        ax5.set_title("% Magnetoresistance")
        ax5.set_xlabel("Field(B)")
        ax5.set_ylabel("(R-R0)/R0 %")
        
        ax6 = plt.subplot(436)
        plt.plot(Bp[np.where(np.absolute(Bp)<9)], (Rsym[np.where(np.absolute(Bp)<9)]*geo-R0)/R0)
        ax6.set_title("R Symmeterized")
        ax6.set_xlabel("Field(B)")
        ax6.set_ylabel("Resistance(Ohm)")
        ax6.set_xlim([-2, 2])
        plt.show
        
        
        
        
        ax7 = plt.subplot(437)
        plt.plot(Bp, (Rxysym))
        ax7.set_title("Symmetrized Hall")
        ax7.set_xlabel("Field(B)")
        ax7.set_ylabel("R-hall")
        
        ax8 = plt.subplot(438)
        plt.plot(Bp[np.where(np.absolute(Bp)<9)], Rxy_sub[np.where(np.absolute(Bp)<9)])
        ax8.set_title("Anamolous Hall")
        ax8.set_xlabel("Field(B)")
        ax8.set_ylabel("R-hall")
        
        
        ax9 = plt.subplot(439)
        s1 = 'R(0) = ' + str(np.round(R0, decimals=1)) + ' Ohm'
        plt.text(.05, .8, s1)
        s2 = 'High Field Slope = '+ str(np.round(z[0], decimals=1))
        plt.text(.05, .6, s2)
        s3 = 'N 2D = ' + "{:e}".format(n2d)
        plt.text(.05, .4, s3)
        s4 = 'Mobility = ' + str(abs(np.round((1/(rho_0/thickness))/n2d/(1.602*10**-19), decimals=1))) + ' cm2/Vs'
        plt.text(.05, .2, s4)
        ax9.get_xaxis().set_visible(False)
        ax9.get_yaxis().set_visible(False)
        
        ax10 = plt.subplot(4,3,10)
        Babs = np.absolute(Bp)
        plt.plot(1/Bp[np.where(((Bp)>1.5) & ((Bp)<(Bmax-0.2)))], drhoBoutput[np.where(((Bp)>1.5) & ((Bp)< (Bmax-0.2)))])
        ax10.set_xlim(left=0)
        ax10.set_title("dRxx/dB")
        ax10.set_ylabel("dRxx/Db")
        ax10.set_xlabel("1/B")
        
        plt.tight_layout(pad=.0, h_pad=.0, w_pad=None, rect=None)
        plt.rcParams.update({'font.size': 30})
        figname=filename.strip('.dat')+'_{:}K'.format(templist[i])
        plt.savefig(figname+'.png', dpi=400)
        plt.savefig(figname+'.svg')
        
        ################# COMPILED PLOT FOR ALL SCANS ####################################
        
        if scans == 1:       #only plots if you are doing a tempscan.  This is done stupidly and just adds the next plot over the last one. 
                             #This is because I didn't use a counter (i) to plot them, so it should be redone but works for nwo  
            
            plt.figure(3)
            plt.rcParams.update({'font.size': 6})
            ax1 = plt.subplot(431)
            plt.plot(Bp, Rxx)
            ax1.set_title("Raw Rxx")
            ax1.set_xlabel("Field(B)")
            ax1.set_ylabel("Resistance(Ohm)")
            
            ax2 = plt.subplot(432)
            plt.plot(Bp, Ryy)
            ax2.set_title("Raw Ryy")
            ax2.set_xlabel("Field(B)")
            ax2.set_ylabel("Resistance(Ohm)")
            
            ax3 = plt.subplot(433)
            plt.plot(Bp, Rxy)
            ax3.set_title("Raw Hall")
            ax3.set_xlabel("Field(B)")
            ax3.set_ylabel("Resistance(Ohm)")
            
            ax4 = plt.subplot(434)
            plt.plot(Bp, Rsym*geo)
            ax4.set_title("R Symmeterized")
            ax4.set_xlabel("Field(B)")
            ax4.set_ylabel("Resistance(Ohm)")
            
        
            ax5 = plt.subplot(435)
            plt.plot(Bp, ((Rsym*geo-R0)/R0)*100)
            ax5.set_title("% Magnetoresistance")
            ax5.set_xlabel("Field(B)")
            ax5.set_ylabel("(R-R0)/R0 %")
            
            ax6 = plt.subplot(436)
            plt.plot(Bp[np.where(np.absolute(Bp)<9)], (Rsym[np.where(np.absolute(Bp)<9)]*geo-R0)/R0)
            ax6.set_title("R Symmeterized")
            ax6.set_xlabel("Field(B)")
            ax6.set_ylabel("Resistance(Ohm)")
            ax6.set_xlim([-2, 2])
            plt.show
            
            
            
            
            ax7 = plt.subplot(437)
            plt.plot(Bp, (Rxysym))
            ax7.set_title("Symmetrized Hall")
            ax7.set_xlabel("Field(B)")
            ax7.set_ylabel("R-hall")
            
            ax8 = plt.subplot(438)
            plt.plot(Bp[np.where(np.absolute(Bp)<9)], Rxy_sub[np.where(np.absolute(Bp)<9)])
            ax8.set_title("Anamolous Hall")
            ax8.set_xlabel("Field(B)")
            ax8.set_ylabel("R-hall")
            
            
            ax9 = plt.subplot(439)
            s1 = 'R(0) = ' + str(np.round(R0, decimals=1)) + ' Ohm'
            plt.text(.05, .8, s1)
            s2 = 'High Field Slope = '+ str(np.round(z[0], decimals=1))
            plt.text(.05, .6, s2)
            s3 = 'N 2D = ' + "{:e}".format(n2d)
            plt.text(.05, .4, s3)
            s4 = 'Mobility = ' + str(abs(np.round((1/(rho_0/thickness))/n2d/(1.602*10**-19), decimals=1))) + ' cm2/Vs'
            plt.text(.05, .2, s4)
            ax9.get_xaxis().set_visible(False)
            ax9.get_yaxis().set_visible(False)
            
            ax10 = plt.subplot(4,3,10)
            Babs = np.absolute(Bp)
            plt.plot(1/Bp[np.where(((Bp)>1.5) & ((Bp)<(Bmax-0.2)))], drhoBoutput[np.where(((Bp)>1.5) & ((Bp)< (Bmax-0.2)))])
            ax10.set_xlim(left=0)
            ax10.set_title("dRxx/dB")
            ax10.set_ylabel("dRxx/Db")
            ax10.set_xlabel("1/B")
            
            
            # ax11 = plt.subplot(4,3,11)
            # plt.plot(Bp, Rxx_diff)
            # ax11.set_title("Delta Magnetoresistance")
            # ax11.set_ylabel("Delta Ohms")
            # ax11.set_xlabel("B")
            
            # ax11 = plt.subplot(4,3,12)
            # plt.plot(Bup, Rxy_diff)
            # ax11.set_title("Delta Hall Effect")
            # ax11.set_ylabel("Delta Ohms")
            # ax11.set_xlabel("B")
        
        
    plt.figure(3)
    plt.tight_layout(pad=.0, h_pad=.0, w_pad=None, rect=None)
    plt.rcParams.update({'font.size': 30})
    figname=filename.strip('.dat') + '_compiled'
    plt.savefig(figname+'.svg')
    plt.savefig(figname+'.png', dpi = 400)
        
        
        
        
        
        
        
        
        
        
        
        
    # #plot dG/(e^2/h)
    # plt.figure()
    # plt.plot(Bp,dGe2h)
    # plt.plot(Bp,dGe2h_approx)
    # plt.show()    
    ##############################################################
    # Export each piece of the split data
    head = 'T(K),B(T),1/B, Rxy(ohms),Rxy-background(ohm), rho(ohm*cm), Magneto,dRxy/dB,dRsqr/dB,drho/dB,Rave(ohm/sq),Rsym(ohm),Rxy_raw,Rxx_raw,Ryy_raw'
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
    
    return final_data