#!/usr/bin/env python
# -*- coding: utf-8 -*-
#This program includes a few functions, one of which will give you the diffusion coeficent for a given
# set of pitch angles (alpha), kinetic energy of the particles we're resonating with (E_kin),
# at a specific L-shell (L), magnetic latitude (lat), alpha star parameter (astar), and wave amplitude(dB)
#def getDs(alpha, E_kin, L, lat, astar, dB):
#the other function will give you the bounced averaged diffusion coeficient for a set of pitch angles and
#energies like above. This assumes a dipole field. This code is from the Summers 2003, 2005, 2007a, & 2007b papers):


# import some general useful stuff
import plasmaconst as pc
import numpy as np
from scipy.integrate import fixed_quad
import scipy.special as spec
from scipy.special.orthogonal import p_roots
from numpy import sum,isinf, real


#Here we are calling up the plasma constants that we will need
#Here we run the plasma constants to get out the SI and CGS values
const =pc.plasmaSI()
constcgs = pc.plasmaCGS()

#This just stores then the info about who wrote this code and what version it is in
 
__version__ = 0.01
__author__ = 'A.J. Halford, S. Morley, and J. Koller'
#modified/optimized by S. Morley, November 2010
#modified by Alexa Halford, January 2011, I've taken getDS from J.Kollers dicccoeff program 
#and put that in with all the obereved data to find the avereage diffusion due to EMIC waves
#during geomagnetic storms, their phases and quiet times. 
#Here we set the constants that will be used for EMIC waves resonating with electrons in the radiation belts
# set wave parameters
# calculate only for electrons
# change to m_p for protons
m_sig = const['m_e']
# set wave bandwidth
# if set to False, then find all resonances even outside waveband
select_bandwidth = True #for chorus waves
sig = 2.0 # 3.17 size of bandwitdh use 2 for chorus waves. 

# select the waves modes (DS2005, appendix A)
# R-mode waves (whistler mode chorus  and hiss)
s = 1 
#L-mode waves (EMIC) interacting 
#s = -1
#interacting with electrons
lam = -1 
#interacting with protons
#lam = const['epsilon']
# select:if True sum over all resonances if false then choose only backward propagating waves: y<0
sum_res = False #for chorus waves resonating with radiation belt electrons
          #True for EMIC waves resonating with radiation belt electrons



def getDs(alpha, E_kin, mlat, astar, dB, B, sbandwidth, w_center):
    """
    alpha = the pitch angles to get the Daa over for a given latitude
    E_kin = the energy of the particle
    lat = magnetic latitude of the observation
    astar_eq = the alpha star at the equator, I should think about how to have this pass through
    I want this to work for both on it's own and with the bounce average stuff.
    dB = the max wave amplitude
    B = the background magnetic field
    sbandwidth = the semi-bandwidth of the wave
    w_center is the max center frequency of the 
    calculate diffusion coefficient D_alpha,alpha given
    the range of alpha pitch angles in [rad],
    kinetic energy in [MeV]
    Lstar coordinate [1]
    latitude coordiantes
    the cold plasma parameter astar

    returned will be
    Daa, Dap, Dpp  where
    Dap=Dap/p and Dpp/p^2
    """
    if isinstance(alpha, np.ndarray):
        alpha = alpha
    else:
        alpha = np.array([alpha]) #These are just the pitch angles determined
                            #for the latitude passed
    # kinetic energy: unit transformation
    E_kin = E_kin*const['e']*1e6    # kinetic energy in [joule]
    E = E_kin/(m_sig*const['c']**2)  # dimensionless particle energy (Why divided by C**2?)
                      #don't know but this is from page 3 of Summers 2007.a in the
                      #paragraph after equation 8
    gamma = E+1  # Lorentz factor (DS2005, text after eq. 35, Summers 2007 after eq. 8)   
    beta = np.sqrt(E*(E+2))/(E+1) #also comes from DS2005 text after eq 35., Summers 2007
                                  # after equation 8
    mu = np.cos(alpha) #

    a = float(s*lam)/gamma #In DS2005 from Appendix A after equation A1 or equation 10 in
                           #Summers 2007
    #b = (1+const['epsilon'])/astar #this come from Appendix A after equation A1 or eq14
                                   # in Summers 2007

    # calc the equatorial B field assuming a dipole field as well as the wave amplitude to
    # background magnetic field ratio

    B0 = B*np.cos(mlat)**6/np.sqrt(4-3*np.cos(mlat)**2.)
    # ratio of wave amplitude over B field
    R = dB**2/(B**2) 

    # resonant wave domain
    Omega_sigl = const['e']*abs(B)/m_sig #this is the cyclotron freq. of the particle
                                         #we are looking at/resonating with 
    Omega_e = const['e']*abs(B)/const['m_e']#this is the electron cyclotron freq.
                                            # at a specific point in space 
    Omega_e0 = const['e']*abs(B0)/const['m_e']# this is the electron cyclotron freq. at the
		                                    #equator of the L-shell we are located on. 
    Omega_p = const['epsilon']*abs(Omega_e)  #Same for protrons
    Omega_p0 = const['epsilon']*abs(Omega_e0)
    nu = np.sqrt(np.pi)*spec.erf(sig) #this is in the text after equation 32
    xm = w_center*Omega_e0/Omega_sigl  # in text after eq 35 
		#wave gaussian center: omega_m = xm*|Omega_e|
    dx = (sbandwidth/sig)*Omega_e0/Omega_sigl  # gaussian semiwidth: domega = dx*|Omega_e|	

    b = (1+const['epsilon'])/astar #since astar has changed, now redeffining b.#this come from Appendix A after equation A1 or eq14
                                   # in Summers 2007

    # define a1 a2 a3 a4 this comes from appendix A in the DS2005 paper 
    a1 = (2*a+s*(-1+const['epsilon'])-(beta*mu)**2*s*(-1+ \
        const['epsilon']))/(1-(beta*mu)**2)
    a2 = (a**2+2*a*s*(-1+const['epsilon'])-const['epsilon']+ \
        (beta*mu)**2*(b+const['epsilon']))/(1-(beta*mu)**2)
    a3 = (a**2*s*(-1+const['epsilon'])-2*a*const['epsilon'])/(1-(beta*mu)**2)
    a4 = -a**2*const['epsilon']/(1-(beta*mu)**2)

    nPA = len(alpha)


    Daa = np.zeros(nPA)
    Dap = np.zeros(nPA)
    Dpp = np.zeros(nPA)

    for i in range(nPA):
        xall = np.roots([1, a1[i], a2[i], a3[i], a4[i]])
        
        # find real roots and ignore complex numbers
        # 1) real number only: isreal(xall)
        # 2) positive frequency only xall > 0
        # 3) y>0 forward propagating only, y<0 backward propagating waves
        # 4) within waveband
        y = (xall+a)/(beta*mu[i])
        #1?
        if (sum_res == False) & (select_bandwidth == True):
            # use this one for only backward prop, and within waveband
            idx = np.where( np.isreal(xall) & (xall>0) & (y<0) & (xall<xm+sig*dx) & (xall>xm-sig*dx) )[0]
        #2?
        elif (sum_res == False) & (select_bandwidth == False):
            # only backward prop, and but allowing all resonances even outside waveband
            idx = np.where( np.isreal(xall) & (xall>0) & (y<0)  )[0]
        #3?
        elif (sum_res == True) & (select_bandwidth == True):
            # use this one for including all forward+backward propagating waves
            idx = np.where( np.isreal(xall) & (xall>0) & (xall<xm+sig*dx) & (xall>xm-sig*dx) )[0]
        #4?
        else:
            # use this one for allowing all resonances even outside waveband
            idx = np.where( np.isreal(xall) & (xall>0) )[0]
            
        if (len(idx) == 0) & (mu[i] > 1e-16): # then no resonance found 
            Daa[i] = np.NaN
            Dap[i] = np.NaN
            Dpp[i] = np.NaN
        elif mu[i] < 1e-16: # then alpha = 90 deg; use limiting equation
            #this comes from the very last paragraph in Appendix A of Summers2005
            x0 = 1/gamma
            y0 = 1/gamma*np.sqrt(1+b*gamma**2/(gamma-1)/(1+const['epsilon']*gamma))
            #Dcore is the similar bit in equations 36 - 38 in Summers2005
            Dcore =  0.5*np.pi/nu*Omega_sigl**2/abs(Omega_e)/(E+1)**2 * \
                R/dx * np.exp( -(x0-xm)**2/dx**2 )
            #These are equations 36 - 38 from Summers2005
            if (sum_res == False):
                Daa[i] = Dcore
                Dap[i] = Dcore/beta*(x0/y0)
                Dpp[i] = Dcore/beta**2 * (x0/y0)**2
            else:  # sum over all resonances
                Daa[i] = 2.0*Dcore
                Dap[i] = np.NaN # the two resonances cancel each other
                Dpp[i] = 2.0*Dcore/beta**2 * (x0/y0)**2
        else:
            #Here we are calculating the critical values where a resonance line is tangential to a dispersion curve.
            #In other words, see Appendix C of Summers2005.
            # calc y with these roots;
            x = xall[idx]
            #print alpha[i]*180/pi, xall
            y = (xall[idx]+a)/beta/mu[i]
            # calc F with this x,y pair
            #C's found in Appendix C
            c1 = 2*s*(-1+const['epsilon'])
            c2 = 1-4*const['epsilon']+const['epsilon']**2
            c3 = -s*(-1+const['epsilon'])*(b+4*const['epsilon'])/2
            c4 = const['epsilon']*(b+const['epsilon'])
            #This comes from Appendix C equation right after equation C1
            g = x**4 + c1*x**3 + c2*x**2 + c3*x + c4
            #This comes from Appendix C equation C1
            F = y*(x-s)**2*(x+s*const['epsilon'])**2/(x*g)
            #Now we're back to the equations 36 - 28 in the Summers 2005 paper. 
            fac = 0.5*np.pi/nu*Omega_sigl**2/abs(Omega_e)/(E+1)**2
            
            Daa_temp = fac * R * (1-x*mu[i]/y/beta)**2*abs(F)/( dx * abs(beta*mu[i]-F)) \
                * np.exp( -(x-xm)**2/dx**2 )

            Dap_temp = -Daa_temp/beta*(x/y) * np.sin(alpha[i])
            Dpp_temp = Daa_temp/beta**2 * (x/y)**2 * np.sin(alpha[i])**2
            
            # sum over all resonances given in xall[idx]
            Daa[i] = np.real(np.sum(Daa_temp))
            Dap[i] = np.real(np.sum(Dap_temp))
            Dpp[i] = np.real(np.sum(Dpp_temp))

    return Daa, Dap, Dpp
   
   
# -----------------------------------------------
def getDs_ba(alpha_eq, E_kin, astar, dB, B0, sbandwidth, w_center, loop):
    """
    bounce averaged diffusion coefficients
    alpha0 ... equatoria PA in [rad]
    E_kin .... kinetic energy in [MeV]
    L      ... L value    
    astar .... cold plasma parameter
    dB    .....wave amplitude (assumed constant along field line
    B0 ........Equatorial field line
    sbandwidth The semi-bandwidth of the wave
    w_center the center frequency (often taken as the frequency where the peak amplitude is
    
    """

    # calculate bounce period, not sure this is correct. In the paper they don't divide by 4 .
    #tau = (1.3-0.56*np.sin(alpha_eq))
    # or more precisely
    tau = 1.38-0.32*( np.sin(alpha_eq) + np.sqrt(np.sin(alpha_eq)))

    # get mirror latitute and setup integration grid 
    # along field line
    lambda_m = [mirrlat(alpha_eq[i]) for i in range(len(alpha_eq))] # returns [rad]
    #lambda_m = min(lambda_m, 25/180*pi);  % equat. confinement of waves
    aveDaa = []
    avea_eq = []
    ave_area = []
    ave_area2 = []
    #Here we are just defining the array length for the bounce average Daa will be calculated
    Daa_ba = np.zeros(len(alpha_eq))
    Dap_ba = np.zeros(len(alpha_eq))
    Dpp_ba = np.zeros(len(alpha_eq))
    for i in range(len(alpha_eq)):
        #print 'this I think where the error is n = ', loop
        Daa_ba[i], dummy = fixed_quad_me(lambda lat: integ(lat, alpha_eq[i], E_kin, astar, dB, B0, sbandwidth, w_center, itype='Daa'), \
                    0, lambda_m[i], n = loop)
        Dap_ba[i], dummy = fixed_quad_me(lambda lat: integ(lat, alpha_eq[i], E_kin, astar, dB, B0, sbandwidth, w_center, itype='Dap'), \
                    0, lambda_m[i], n = loop)
        Dpp_ba[i], dummy = fixed_quad_me(lambda lat: integ(lat, alpha_eq[i], E_kin, astar, dB, B0, sbandwidth, w_center, itype='Dpp'), \
                    0, lambda_m[i], n = loop)
        

    v = np.sqrt(2.*E_kin/(m_sig))#Here we are just getting the velocity of the particle E_kin = 1/2 m v^2.
    E = E_kin/(m_sig*const['c']**2)  # dimensionless particle energy (Why divided by C**2?)
                      #don't know but this is from page 3 of Summers 2007.a in the
                      #paragraph after equation 8
    
    Daa_ba = Daa_ba/tau
    Dap_ba = Dap_ba/tau
    Dpp_ba = Dpp_ba/tau

    return Daa_ba, Dap_ba, Dpp_ba #The Dap and Dpp really are Dap/p and Dpp/p^2.  I don't need to find p
                                #, Dap_ba, Dpp_ba

# -----------------------------------------------
def integ(lat, alpha0, E_kin, astar_eq, dB, B0, sbandwidth, w_center, itype):
    alpha_eq =alpha0
    if isinstance(lat, np.ndarray):
        lat = lat
    else:
        lat = np.array([lat])
    #Check this equation to make sure that locPA is sin of the local pitch angle
    #locPA = np.array([np.sin(alpha_eq)*(4-3*np.cos(lat[i])*np.cos(lat[i]))**0.25 / np.cos(lat[i])**3 for i in range(len(lat))])

    sinlocPA = np.array([np.sin(alpha_eq)*(1. + 3.*np.sin(lat[i])**2.)**0.25/(np.cos(lat[i])**3.) for i in range(len(lat))]) 
    # get the diffusion coefficients
    #Here we are makeing sure that all sin(alpha) <= 1
    #If greater than 1 we're defining it as 1
    #then turning this back into an angle. 
    nlat = len(lat)
    sinlocPA[np.where(sinlocPA > 1.0)] = 1.0
    locPA = np.arcsin(sinlocPA)
    B = B0*np.sqrt(1.+ 3.*np.sin(lat)**2.)/np.cos(lat)**6
    Omega_sigl = const['e']*abs(B)/m_sig #this is the cyclotron freq. of the particle
                                         #we are looking at/resonating with 
    Omega_e = const['e']*abs(B)/const['m_e']#this is the electron cyclotron freq.
                                            # at a specific point in space 
    Omega_e0 = const['e']*abs(B0)/const['m_e']
    astar = astar_eq*Omega_sigl**2/Omega_e0**2
    #print 'nlat', nlat    
    Daa = np.zeros(nlat)
    Dap = np.zeros(nlat)
    Dpp = np.zeros(nlat)
    
    #print 'going over these latitudes', lat*180./np.pi
    for i in range(nlat):
        Daa[i], Dap[i], Dpp[i] = getDs(locPA[i], E_kin, lat[i], astar[i], dB, B[i], sbandwidth, w_center)
                
    Daa[np.where(np.isnan(Daa))] = np.nan#1e-16#np.nanmin(Daa)
    Dap[np.where(np.isnan(Dap))] = np.nan#1e-16#np.nanmin(Dap)
    Dpp[np.where(np.isnan(Dpp))] = np.nan#1e-16#np.nanmin(Dpp)#0#1e-16
    

    if itype == 'Daa':
        In = Daa* np.cos(locPA)*np.cos(lat)**7/np.cos(alpha_eq)**2
    elif itype == 'Dap':
        In = Dap* np.sin(locPA)*np.cos(lat)**7/np.cos(alpha_eq)/np.sin(alpha_eq)
    elif itype == 'Dpp':
        In = Dpp* np.sin(locPA)**2 * np.cos(lat)**7/np.cos(alpha_eq)/np.sin(alpha_eq)**2
    else:
        raise TypeError('Type not supported')
    
    return In
    
def test_ave(lat, alpha0, E_kin, astar_eq, dB, B0, sbandwidth, w_center, itype):
    alpha_eq =alpha0
    if isinstance(lat, np.ndarray):
        lat = lat
    else:
        lat = np.array([lat])
    #Here we are just finind the local pitch angles along the bounce path. 
    sinlocPA = np.array([np.sin(alpha_eq)*(1. + 3.*np.sin(lat[i])**2.)**0.25/(np.cos(lat[i])**3.) for i in range(len(lat))]) 
   
    #Here we are makeing sure that all sin(alpha) <= 1
    #If greater than 1 we're defining it as 1
    #then turning this back into an angle. 
    nlat = len(lat)
    sinlocPA[np.where(sinlocPA > 1.0)] = 1.0
    locPA = np.arcsin(sinlocPA)
    #Here we are getting the local B along the field line assuming a dipole field. Think about defining the local pitch angle
    #with respect to the local B-field instead of the other way above so that I can plug in any magnetic field model. 
    B = B0*np.sqrt(1.+ 3.*np.sin(lat)**2.)/np.cos(lat)**6
    Omega_sigl = const['e']*abs(B)/m_sig #this is the cyclotron freq. of the particle
                                         #we are looking at/resonating with 
    Omega_e = const['e']*abs(B)/const['m_e']#this is the electron cyclotron freq.
                                            # at a specific point in space 
    Omega_e0 = const['e']*abs(B0)/const['m_e']
    astar = astar_eq*Omega_sigl**2/Omega_e0**2

    #This is just creating the arrays that will be filled with the Daa's for this starting equatorial pitch angle.  
    Daa = np.zeros(nlat)
    Dap = np.zeros(nlat)
    Dpp = np.zeros(nlat)

    # get the diffusion coefficients
    for i in range(nlat):
        #print 'local pitch angle', locPA[i], 'E_kin', E_kin, 'lat', lat[i], 'astar', astar,'dB', dB, B0, sbandwidth, w_center
        Daa[i], Dap[i], Dpp[i] = getDs(locPA[i], E_kin, lat[i], astar[i], dB, B[i], sbandwidth, w_center)
           
    Daa[np.where(np.isnan(Daa))] = np.nan#1e-16#np.nanmin(Daa)
    Dap[np.where(np.isnan(Dap))] = np.nan#1e-16#np.nanmin(Dap)
    Dpp[np.where(np.isnan(Dpp))] = np.nan#1e-16#np.nanmin(Dpp)#0#1e-16


    #Here I'm going to get the "area underneath the curve".
    lat_step = lat[1] - lat[0]
    area = 0
    area2 = 0
    for i in range(nlat):
        if np.isfinite(Daa[i]):
            area = area + (Daa[i]* np.cos(locPA[i])*np.cos(lat[i])**7.*lat_step/np.cos(alpha_eq)**2.)
            area2 = area2 + (Daa[i]*np.sin(alpha_eq)*np.cos(lat[i])*(1.+3.*np.sin(lat[i])**2.)**0.25)/(np.cos(alpha_eq)*np.sin(locPA[i]))
    if itype == 'Daa':
        In = Daa* np.cos(locPA)*np.cos(lat)**7
    elif itype == 'Dap':
        In = Dap* np.sin(locPA)*np.cos(lat)**7/np.cos(alpha_eq)/np.sin(alpha_eq)
    elif itype == 'Dpp':
        In = Dpp* np.sin(locPA)**2 * np.cos(lat)**7/np.cos(alpha_eq)/np.sin(alpha_eq)**2
    else:
        raise TypeError('Type not supported')

    if len(Daa) > 1:

        aveDaa = np.nansum(In)
        avea_eq = alpha0
    
    return aveDaa, avea_eq, area, area2

# -----------------------------------------------
def mirrlat(alpha0):
    """
    calculate the mirror latutude given the equatorial
    PA alpha0 in [rad]
    returns: lam in [rad]


    can be found by solving a sixth order equation
    x^6 + (3*sin^4(a0))x - 4*sin^4(a0) = 0
    where x = cos(lambda_m)^2

    Okay - Use equations 19, 20, and 21 from Summers 2007 (paper 1). 
    """

    c1 = 3*(np.sin(alpha0))**4
    c0 = -4*(np.sin(alpha0))**4

    # get roots
    xall = np.roots([ 1, 0, 0, 0, 0, c1, c0])

    # find real and positive values
    idx = np.where( np.isreal(xall) & (xall>0) )[0]

    # get mirror latitude
    mlat = np.arccos(np.sqrt(xall[idx]))

    return np.real(mlat)

def fixed_quad_me(func, a, b, args=(), n=5):
    """
    My version of fixed quad which will then be able to sum over NaNs. 
    """
    [x,w] = p_roots(n)
    x = real(x)
    ainf, binf = map(isinf, (a,b))
    if ainf or binf:
        raise ValueError("Gaussian quadrature is only avalible for finite limits.")
    y = (b-a)*(x+1)/2.0 + a
    return (b-a)/2.0*np.nansum(w*func(y,*args),0), None



    
def getemin(B, w_center, den):
    """
    alpha = the pitch angles to get the Daa over for a given latitude
    E_kin = the energy of the particle
    lat = magnetic latitude of the observation
    astar_eq = the alpha star at the equator, I should think about how to have this pass through
    I want this to work for both on it's own and with the bounce average stuff.
    dB = the max wave amplitude
    B = the background magnetic field
    sbandwidth = the semi-bandwidth of the wave
    w_center is the max center frequency of the 
    calculate diffusion coefficient D_alpha,alpha given
    the range of alpha pitch angles in [rad],
    kinetic energy in [MeV]
    Lstar coordinate [1]
    latitude coordiantes
    the cold plasma parameter astar

    returned will be
    Daa, Dap, Dpp  where
    Dap=Dap/p and Dpp/p^2
    """
    
    q = const['e']
    c = const['c']
    me = const['m_e']
    

    omega_e =q*B/(const['m_e'])
    omega_p = den*q**2./(const['ep_0']*me)

    k2 = (w_center**2./c**2.) - ((omega_p**2/c**2.)/(1-(omega_e/w_center)))
    k = np.sqrt(k2)

    a = k2 + ((den**2.*omega_e**2.)/(c**2.))
    b = 2.*k*w_center
    d = w_center**2.-(den**2.*omega_e)

    vpos = (b+np.sqrt(b**2. -4.*a*d))/(2.*a)
    vneg = (b-np.sqrt(b**2. -4.*a*d))/(2.*a)
    Epos = 1.5*me*vpos**2.
    Eneg = 1.5*me*vneg**2.

    print( 'Epos = ', Epos*6.242*10**18)
    print( 'Eneg = ', Eneg*6.242*10**18)



    return np.nanmin([Epos*6.242*10**18,  Eneg*6.242*10**18])
