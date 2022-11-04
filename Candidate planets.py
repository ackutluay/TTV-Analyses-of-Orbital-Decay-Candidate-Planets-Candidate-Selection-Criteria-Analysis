# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 13:08:19 2022

@author: ackut
"""
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors
from astropy import constants as const

#exoplanet catalogue
nasaexo = pd.read_csv('nasa_exoplanet_archive_catalog_composite_20220113.csv',
                      comment="#", skipinitialspace=True)
# For Brown Dwarfs --> TEPCat
tepcat = pd.read_csv('tepcat_20220310.csv',
                     comment="#", skipinitialspace=True)

def t_a(m_p,m_s,a_0,r_s): #e-folding timescale for orbital decay, for one 
    M_p= m_p * const.M_jup.value #conversion of planetary mass to kg
    M_s= m_s * const.M_sun.value #conversion of star mass to kg
    a= (M_s / M_p)
    
    A= a_0 * const.au.value #conversion of a to meters
    R_s= r_s * const.R_sun.value #conversion of star radius to meters
    b= (A/R_s)**5
    return np.log10(a*b)

def t_s(m_p,m_s,a_0,r_s): #spin synchronization timescale of the spin of the orbit for one
    M_p= const.M_jup.value * m_p #conversion of planetary mass to kg
    M_s= const.M_sun.value * m_s #conversion of star mass to kg
    a= (M_s / M_p)**2 
    
    A= const.au.value * a_0 #conversion of a to meters
    R_s= const.R_sun.value * r_s #conversion of star radius to meters
    b= (A/R_s)**3       # ratio of semi-major axis to orbital radius
    return np.log10(a*b)

# Transiting planets are from NASA Exoplanet Archive
cond = (nasaexo['discoverymethod'] == 'Transit') & (
    nasaexo['pl_bmassj'] < 13.0)
tranet = nasaexo[cond]

# Brown dwarfs are from TEPCat
cond_BD = (tepcat['Type'] == 'BD')
BD = tepcat[cond_BD]
 
#CALCULATE THE TAU_S AND TAU_A FOR HATS-18, WASP-12 AND QATAR-4

#print(nasaexo['pl_name'][4665]) #for wasp-12b
#print(nasaexo['pl_name'][367]) #for hats-18b
#print(nasaexo['pl_name'][4483]) #for qatar-4b
#print orbital period values
print('ORBITAL PERIOD VALUES'),print('WASP-12b=',"{:.3f}".format(nasaexo['pl_orbper'][4665])),
print('HATS-18b=',"{:.3f}".format(nasaexo['pl_orbper'][367])),
print('WASP-12b=',"{:.3f}".format(nasaexo['pl_orbper'][4483]))

#now we are using the parameters pl_bmassj, st_rad, st_mass, pl_orbsmax
#these are the e-folding parameters of 4 planets
w_t= t_a(nasaexo['pl_bmassj'][4665],nasaexo['st_mass'][4665], #wasp-12 t_a
         nasaexo['pl_orbsmax'][4665],nasaexo['st_rad'][4665])   

h_t= t_a(nasaexo['pl_bmassj'][367],nasaexo['st_mass'][367], #hats-18 t_a
         nasaexo['pl_orbsmax'][367],nasaexo['st_rad'][367])

q_t= t_a(nasaexo['pl_bmassj'][4483],nasaexo['st_mass'][4483], #qatar-4 t_a
         nasaexo['pl_orbsmax'][4483],nasaexo['st_rad'][4483]) 

k_t=t_a(nasaexo['pl_bmassj'][1542],nasaexo['st_mass'][1542], #qatar-4 t_a
         nasaexo['pl_orbsmax'][1542],nasaexo['st_rad'][1542])

print('E-FOLDING VALUES'),print('WASP-12=',"{:.3f}".format(w_t)),print('HATS-18=',"{:.3f}".format(h_t)),
print('Qatar-4=',"{:.3f}".format(q_t)),print('KELT-1=',"{:.3f}".format(k_t))

opl=[nasaexo['pl_orbper'][4665],nasaexo['pl_orbper'][367],nasaexo['pl_orbper'][1542]]
list1= np.array([w_t,h_t,q_t])
#these are the spin_synchronization of 4 planets
w_s= t_s(nasaexo['pl_bmassj'][4665],nasaexo['st_mass'][4665], #wasp-12 t_s
         nasaexo['pl_orbsmax'][4665],nasaexo['st_rad'][4665])   

h_s= t_s(nasaexo['pl_bmassj'][367],nasaexo['st_mass'][367], #hats-18 t_s
         nasaexo['pl_orbsmax'][367],nasaexo['st_rad'][367])

q_s= t_s(nasaexo['pl_bmassj'][4483],nasaexo['st_mass'][4483], #qatar-4 t_s
         nasaexo['pl_orbsmax'][4483],nasaexo['st_rad'][4483]) 

k_s=t_s(nasaexo['pl_bmassj'][1542],nasaexo['st_mass'][1542], #qatar-4 t_s
         nasaexo['pl_orbsmax'][1542],nasaexo['st_rad'][1542])

list2= np.array([w_s,h_s,q_s])
print('SPIN SYNCHRONIZATION TIMESCALE VALUES'),print('WASP-12=',"{:.3f}".format(w_s)),print('HATS-18=',"{:.3f}".format(h_s)),
print('Qatar-4=',"{:.3f}".format(q_s)),print('KELT-1=',"{:.3f}".format(k_s))
#print(k_t,k_s) #just to check the values for KELT-1 from Baştürk et al(2022)


#NOW PLOT ORBITAL RADIUS (pl_orbper) vs t_a & t_s values
#first calculate for all catalogue
#first use the constants

jup_m= const.M_jup.value #jupiter mass
sun_m = const.M_sun.value #sun mass
sma= const.au.value #1 au in meters
s_rad= const.R_sun.value
#for tepcat
q2= (BD['R_A']*sun_m) / (BD['M_b']*jup_m)
j2= (BD['a(AU)']*sma) / (BD['R_A']*s_rad)

t2= np.log10(q2* (j2**5))
s2= np.log10((q2**2)*(j2**3))
 
#for Nasa exoplanetary archive
q3= (tranet['st_mass']*sun_m) / (tranet['pl_bmassj']*jup_m)
j3= (tranet['pl_orbsmax']*sma) / (tranet['st_rad']*s_rad)

t3= np.log10(q3* (j3**5))
s3= np.log10((q3**2)*(j3**3))

#now plot all of these things
plt.rc('font', size=12)
plt.rc('axes', labelsize= 10)
plt.rc('figure', titlesize= 16)
plt.rc('xtick', labelsize=12)    # fontsize of the tick labels
plt.rc('ytick', labelsize=12)    # fontsize of the tick labels
fig=plt.figure(figsize=(9,9))
plt.scatter(tranet['pl_orbper'],t3,marker='.',
            cmap='Greys',c=tranet['pl_bmassj'],
            norm=colors.PowerNorm(gamma=0.125),s=tranet['pl_radj']*120,)
plt.scatter(BD['Period(day)'],t2,color='brown',marker='*',s=BD['R_b']*120)
plt.ylim((3,15))
plt.xlim((0.5,50))
plt.xscale('log')
plt.xlabel('Orbital Period (days) ')
plt.ylabel('${\u03C4}_{a}$')
plt.text(1.1,5.4,'WASP-12',color='red',fontsize=12)
plt.text(0.9,5.7,'HATS-18',color='red',fontsize=12)
plt.text(1.9,6.2,'Qatar-4',color='red',fontsize=12)
plt.show()

fig=plt.figure(figsize=(9,9))
plt.scatter(tranet['pl_orbper'],s3,marker='.',
            cmap='Greys',c=tranet['pl_bmassj'],
            norm=colors.PowerNorm(gamma=0.125),s=tranet['pl_radj']*120,)
plt.scatter(BD['Period(day)'],s2,color='brown',marker='*')
plt.ylim((3,15))
plt.xlim((0.5,50))
plt.xscale('log')
plt.xlabel('Orbital Period (days)')
plt.ylabel('${\u03C4}_{\omega}$')
plt.text(1.1,7.3,'WASP-12',color='red',fontsize=12)
plt.text(0.9,7.0,'HATS-18',color='red',fontsize=12)
plt.text(1.9,6.7,'Qatar-4',color='red',fontsize=12)
plt.show()







