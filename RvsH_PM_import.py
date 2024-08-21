# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 10:30:27 2021

@author: 5ul
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

#head = 'T(K),B(T),1/B, Rxy(ohms),Rxy-background(ohm), rho(ohm*cm), Magneto,dRxy/dB,dRsqr/dB,drho/dB,Rave(ohm/sq),Rsym(ohm),Rxy_raw,Rxx_raw,Ryy_raw'
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