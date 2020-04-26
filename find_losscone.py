#Here we find the loss cone angle


import numpy as np
import plasmaconst as pc
const =pc.plasmaSI()
constcgs = pc.plasmaCGS()


#Here we are assuming a dipole
def find_losscone(L, H):
    """ This function finds the loss cone for a given L value and a given loss hight.
    To run import find_losscone as fl
    fl.find_losscone(L,H)
    This then outputs the loss cone angle in degrees """
    #L = L-vaue and H is loss height. 
    import numpy as np
    import plasmaconst as pc
    const =pc.plasmaSI()
    constcgs = pc.plasmaCGS()

    R = (H*1000.+const['RE'])/const['RE']

    cos2Mlat = R/L
    cosMlat = np.sqrt(R/L)
    Mlat = np.arccos(cosMlat)

    
    Bo = const['B_E']/(L**3.)
    Bm = Bo*((4.-3.*R/L)**(0.5))/((R/L)**3.)
    sin2theta = Bo/Bm
    sintheta = np.sqrt(sin2theta)
    theta = np.arcsin(sintheta)
    

    return theta*180./np.pi


def strong_diff(L, E):
    """ This function finds the strong difussion limit for a given l value and energy.
    to run type
    $import find_losscone as fl
    $fl.strong_diff(L,E)
    and it outpouts the relavent strong diffusion limit"""
    #E is in MeV
    import numpy as np
    import plasmaconst as pc
    const =pc.plasmaSI()
    constcgs = pc.plasmaCGS()

    frac1 = (0.205*constcgs['c'])/(L**4.*(E+1.)*constcgs['RE'])
    frac2a = (4.*L*E*(E+2.))/(4.*L -3.)
    frac2 = np.sqrt(frac2a)
    Dsd = frac1*frac2

    return Dsd

def dipole_drift(L, E, alpha):
    """ Here we find the dipole drift period for a given L, energy, and pitch angle.
    to run type
    fl.dipole_drift(L, E, alpha)
    and it returns the drift period in minutes"
    """
    #This is from eq 4 in Beharrell et al 2015. 
    import numpy as np
    import plasmaconst as pc
    const =pc.plasmaSI()
    constcgs = pc.plasmaCGS()
    frac = (3.*L*E)/(22.*10**(6.)) #L in RE and E in eV
    inside = (0.7 + 0.3*np.sin((alpha*np.pi/180.))) #alpha is the equatorial pitch angle
    #wd in degrees of magnetic longitude per second
    wd_rate = frac*inside + 4.17*10**(-3.) #4.17*10**(-3) is the contribution from Earth's rotation
    td = 360./wd_rate/60. # in minutes
    printe( td, ' minutes')
    return td
        
def dipole_bounce(L, E, alpha):

    #rest mass of electron
    m_0 = 511.*10**3 / const['c'] / const['c']
    #rest energy of electron in eV
    E_0 = 511.*10**3
    # Kinetic Energy in Rest Energy 
    T = E / E_0
    #relativistic correction factor
    gamma = T + 1
    # relativistic mass correction
    m = gamma * m_0
    
    tb = (L*const['RE']) * (3.7-1.6*np.sin(alpha*np.pi/180.) ) * (E / m)**(-0.5)

    print( tb, ' seconds')
    
    return tb
