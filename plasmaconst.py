#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This creats a dictionary of all the plasma constants used on a regular basis.
to run for example;
import plasmaconstant as pc

pc.plasmaSI() #This will return the values in SI units in dictionary pcsi
pc.plasmaCGS() #This will return the values in CGS units in dictionary pccgs

The values returned are 
'm_e' the electron mass
'm_p' the proton mass
'c' the speed of light
'e' the charge of an electron
'epsilon' the mass ratio of an electon to a proton
'B_E' the equatorial magnetic field at the surface of the Earth (1RE)
'kb' the boltzman constant
'mu_0' the permability of free space
'ep_0' the permittivity of free space
'G' the gravitational constant
'RE' the radius of the Earth
"""
import math as math
#here I am naming the function
def plasmaSI():
  """
  This creates a dictionary of all the plasma constants used on a regular basis.
  All constants are in SI units
  """
  # set physical constants here in SI units
  pcsi = {} #This is how he defines an empty dictionary.
  pcsi['m_e'] = 9.1094e-31     # electron mass in [kg]
  pcsi['m_p'] = 1.6726e-27     # proton mass in [kg]
  pcsi['c'] = 2.99792458e8         # speed of light in [m/s]
  pcsi['e'] = 1.6022e-19        # electron charge in [A.s]
  pcsi['epsilon'] = pcsi['m_e']/pcsi['m_p']  # mass ratio
  pcsi['B_E'] = 3.1255617e-5 # B field at surface = 1RE in [Tesla]
  pcsi['kB'] = 1.380658e-23 	#the boltzman constand in j/k
  pcsi['mu_0'] = math.pi*4.0e-7	# Permability of free space H/m
  pcsi['ep_0'] = 1./(pcsi['c']**2.*pcsi['mu_0']) #Permittivity of free space F/m
  pcsi['G'] = 6.6726e-11 	#Gravitational constant m^3/(s^2kg)
  pcsi['RE'] = 6378.1e3 	#the radius of the Earth in meters.
  return pcsi
 
 
  #Here I am naming the function
def plasmaCGS():
  """
  This creates a dictionary of all the plasma constants used on a regular basis.
  All constants are in CGS units
  """
  # set physical constants here in CGS units
  pccgs = {} #This is how he defines an empty dictionary.
  pccgs['m_e'] = 9.1094e-28     # electron mass in [g]
  pccgs['m_p'] = 1.6726e-24     # proton mass in [g]
  pccgs['c'] = 2.99792458e10       # speed of light in [cm/s]
  pccgs['e'] = 4.8032e-10        # electron charge in [esu]
  pccgs['epsilon'] = pccgs['m_e']/pccgs['m_p']  # mass ratio
  pccgs['B_E'] = 0.31255617 # B field at surface = 1RE in [Gauss]
  pccgs['kB'] = 1.380658e-16 	#the boltzman constand in e/k
  pccgs['mu_0'] = 1		#unitless
  pccgs['ep_0'] = 1 	#unitless
  pccgs['G'] = 6.6726e-8 	#Gravitational constant cm^3/(s^2g)
  pccgs['RE'] = 6378.1e5 	#the radius of the Earth in cm
  return pccgs
