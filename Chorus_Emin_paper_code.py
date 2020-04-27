###############################################################################################
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import find_losscone as LC
import matplotlib as mpl
import math
import matplotlib.dates as mdates
import whistler_daa_emin_diffcoef as dc
import plasmaconst as pc
const =pc.plasmaSI()
constcgs = pc.plasmaCGS()

#This will likely need to be updated to point where the
#cdf libaray is on your computer
import os
os.environ["CDF_LIB"] =  "/Applications/cdf/cdf37_1-dist/lib"


from matplotlib.pyplot import figure 
from matplotlib.colors import LogNorm
from spacepy import pycdf
from datetime import timedelta
from matplotlib.gridspec import GridSpec

###############################################################################################
#Here we are defining the start and stop times and other event times
starttime = dt.datetime(2014, 1, 9, 19, 55)
plttimestart = dt.datetime(2014, 1, 9, 20, 5)
endtime = dt.datetime(2014, 1, 9, 20, 35)
plttimeend = dt.datetime(2014, 1, 9, 20, 20)

compression = dt.datetime(2014, 1, 9, 20, 10, 15)
wave = dt.datetime(2014, 1, 9, 20, 11, 30)

compression_end = dt.datetime(2014, 1, 9, 20, 18, 00)
wave_end = dt.datetime(2014, 1, 9, 20, 14, 45)
##############################################################################################

# these are colorblind friendly on their own and many colormaps work well. But be sure to test on 
#https://www.color-blindness.com/coblis-color-blindness-simulator/

#First we are defining the primary colors used for these maps
cdi = '#093145'
cli = '#3c6478'
cda = '#107896'
cla = '#43abc9'
cdk = '#829356'
clk = '#b5c689'
cdd = '#bca136'
cld = '#efd469'
cdc = '#c2571a'
clc = '#f58b4c'
cdr = '#9a2617'
clr = '#cd594a'
clg = '#F3F4F6'
cdg = '#8B8E95'

#Now we are putting the ones together that we want to go from one to another
greycolors = [clg, cdg]                                          
greencolors = [clg, clk, cdk]
yellowcolors = [clg,cld, cdd]
redcolors = [clg, clr, cdr]
hotcolors = [cld, cdd, cdc, cdr]
colors = [cdi,  cdk, cld,  cdc,  cdr]
bluecolors = [clg, cla, cda, cdi]
trycolors = [cdi, cli, cdd, cld, clg]#[cdi, cli, cla, cld]#[cla, cda, cdi, cdr, cdc, clc]

#finally we make the actual maps
bluermap = mpl.colors.LinearSegmentedColormap.from_list("", bluecolors)
pltmap = mpl.colors.LinearSegmentedColormap.from_list("", hotcolors)
greenmap = mpl.colors.LinearSegmentedColormap.from_list("", greencolors)
yellowmap = mpl.colors.LinearSegmentedColormap.from_list("", yellowcolors)
redmap = mpl.colors.LinearSegmentedColormap.from_list("", redcolors)
trymap = mpl.colors.LinearSegmentedColormap.from_list("", trycolors)

#We chose trymap for this set of figures
colormap = trymap
#Here we are setting the text size for the figures. 
txtsize = 12
mpl.rcParams['font.size'] = txtsize
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#Setup the electron energies we want to test
#Here we are choosing the HOPE energy channels that we're going to plot.
#channels = np.array([6,11,17,20,31])
channels = np.array([7,13,24,33])

# setup pitch angles in [rad]
alpha = np.linspace(0.1, 90, 180)*np.pi/180.0# 180)*np.pi/180.0 #radians

#This section uses the observed values for the period during the chorus wave observations.  
dB = 0.025*10**(-9.)#Tesla's
mlat = 0. #radians
B0 = 167.*10**(-9.) #Teslas
den = 12.
astar = B0**2./(4*np.pi*den*const['m_e']*const['c']**2.)#

wave_freq = 2600.#Hz
sbandwidth = (2.*np.pi*500.)/(const['e']*abs(B0)/const['m_e'])#0.31
w_center =  (2.*np.pi*wave_freq)/(const['e']*abs(B0)/const['m_e'])#0.35


#This is just for plotting
W_Daamin = 10**(-18)
W_Daamax = 10**(0)


#here we are setting the line style/markers for plotting the results. 
style = '--'


#*********************************************************************************
#Here we are going to plot the Chorus wave observed at RBSPB
#**********************************************************************************

print('starting to plot the chorus wave Fig. 1 panel a')
#Here we are reading in the survey data
RBSPB_mag = 'data/rbsp-b_wna-survey_emfisis-L4_20140109_v1.2.2.cdf'
#Now we are setting the start and end times for the interval 
#we want to look at and what we will plot
starttime = dt.datetime(2014, 1, 9, 19, 55)
plttimestart = dt.datetime(2014, 1, 9, 20, 9)
endtime = dt.datetime(2014, 1, 9, 20, 35)
plttimeend = dt.datetime(2014, 1, 9, 20, 20)


#Here we are creating the plot of the wave power spectrogram. 
data = pycdf.CDF(RBSPB_mag) 
time = data['Epoch'][:]
Bsum = data['bsum'][:]
Bmag = Bsum
freq = data['WFR_frequencies'][0,:]
Buvw = data['Buvw'][:]
B = np.sqrt(Buvw[:,0]**2 + Buvw[:,1]**2 + Buvw[:,2]**2)

#Here are the times for the compression and all in mdate format for the plot labels
plttime = mdates.date2num(time)
starttime = mdates.date2num(plttimestart)
endtime = mdates.date2num(plttimeend)
comptime = mdates.date2num(compression)
wavetime = mdates.date2num(wave)
# 
# Now we are creating panel a) of figure 1
fig=plt.figure(figsize=(18.,11.69))
gs=GridSpec(6,2, hspace=0.3, wspace = 0.1)
ax=fig.add_subplot(gs[0:1,:])
xx, yy = np.meshgrid(plttime, freq)
plt.pcolormesh(xx,yy,Bmag.T, cmap = colormap,norm = LogNorm(vmin = 10**(-10), vmax = 10**(-4)))
ax.set_ylim(10**2,10**4)
ax.set_xlim(starttime, endtime)
cb = plt.colorbar()
plt.ylabel('Frequency Hz')
ax.set_yscale('log')
cb.set_ticks([10**(-10), 10**(-8), 10**(-6),10**(-4)])
cb.set_ticklabels([-10, -8, -6, -4 ])
cb.ax.set_ylabel(r'Log(nT$^{2}$/Hz)')
# 
# Finding the electron cyclotron frequency
ec = (const['e']*B*10**-9)/(2.*np.pi*const['m_e'])
# 
# Here we are overplotting the electron cycltron frequency and the half electron cyclotron frequency
ax.plot(plttime, ec, '--', c = clg, label = r'a)')
ax.plot(plttime, ec/2., '--', c = clg)
ax.plot((comptime,comptime), (10**1,10**4), linestyle = style, color = 'gold', linewidth = 2.)
ax.plot((wavetime,wavetime), (10**1,10**4), linestyle = style, color = 'white', linewidth = 2.)
plt.xlabel('9 January 2014')
date_format = mdates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(date_format)
leg = ax.legend(loc = 'upper right',  frameon=False, prop={'size': txtsize})
for text in leg.get_texts():
	text.set_color("w")
for item in leg.legendHandles:
	item.set_visible(False)

#now we are opening and plotting the HOPE data. 
HOPE_file = 'data/rbspb_rel03_ect-hope-PA-L3_20140109_v6.1.0.cdf'
data = pycdf.CDF(HOPE_file)
time = data['Epoch_Ele'][:]
index = np.where((time >= plttimestart ) & (time <= plttimeend))

#We are getting the flux and pitch angle data
HOPE_alpha =data['PITCH_ANGLE'][:]
energy = data['HOPE_ENERGY_Ele'][0,:]
electrons_data = data['FEDU'][:]
event = electrons_data[index, :, :]
electrons = event[0,:,:,:]

plttime = mdates.date2num(time)
#Making sure we are removing all the bad data
bad = np.where(electrons < 0.)
electrons[bad] = np.nan
#Now we plot it all
xx, yy = np.meshgrid(plttime[index], HOPE_alpha)# was used to try to center the bins-np.array([2.25, 6.75, 9, 9, 9, 4.5,-9, -9, -9, -6.75, -2.25]))

E_kin = energy[channels]/10**(6) # in MeV

panels_l = ['b)', 'd)', 'f)', 'h)']
panels_r = ['c)', 'e)', 'g)', 'i)']
daa_plt = np.zeros((len(E_kin), len(alpha)))
daab_plt = np.zeros((len(E_kin), len(alpha)))
# 
# Here we are now getting the Daa and <Daa> values and plotting them along with the observed data
for i in range(len(E_kin)):
    print( 'on loop', i, 'E_kin = ', np.int(E_kin[i]*10**6), 'eV')

#    Here we are getting the local diffusion coefficients. 
    daa, dap, dpp = dc.getDs(alpha, E_kin[i], mlat, astar, dB, B0, sbandwidth, w_center)
#    Now we get the bounced averaged stuff. The 180 is just so the number of loops in the integration. 
    daa_b, dap_b, dpp_b = dc.getDs_ba(alpha, E_kin[i], astar, dB, B0, sbandwidth, w_center, 180)
    print( '<daa> = ', daa_b)
    emin = dc.getemin(B0, w_center, den)
    print( 'E_min', emin)
    print( 'emin is ', emin, ' eV')
#    now we are plotting both
    daa_plt[i,:] = daa
    daab_plt[i,:] = daa_b
    
    dx=fig.add_subplot(gs[i+2,0])
    dx.plot(daa, alpha*180./np.pi, color = cli)
    dx.plot(daa, 180 - alpha*180./np.pi, color = cli)
    dx.plot(daa_b, alpha*180./np.pi, color = cdi, linestyle = style)
    dx.plot(daa_b, 180 - alpha*180./np.pi, color = cdi, linestyle = style)
    dx.plot((10**(-6),10**(0)), (90,90), c = 'k', linestyle = style, linewidth = 2., label = panels_l[i])
    dx.set_ylabel(r' $\alpha_{eq}$'+' \n'+np.str(np.int(energy[channels[i]]))+' eV')
    dx.set_ylim(0, 180)
    dx.set_xscale("log", nonposx='clip')
    dx.set_xlim(1./(10.*60.), 10**(-1) )
    leg = dx.legend(loc = 'upper right', frameon=False, prop={'size': txtsize})
    for text in leg.get_texts():
    	plt.setp(text, color = cdi)
    for item in leg.legendHandles:
    	item.set_visible(False)    
    if i != len(E_kin)-1:
        dx.xaxis.set_visible(False)
    print('finished plotting panel', panels_l[i])    

    ax=fig.add_subplot(gs[i+2, 1])
    plt.pcolormesh(xx,yy,(electrons[:,:,channels[i]].T)/(0.7*np.nanmax(electrons[:,:,channels[i]].T)), cmap = colormap,norm = LogNorm(vmin = 0.01, vmax = 1))
    
    plt.ylim(0,180)
    plt.plot((starttime,endtime), (90,90), c = 'w', linestyle = style, linewidth = 2., label = panels_r[i])
    plt.plot((comptime,comptime), (0,180), c = 'w', linestyle = style, linewidth = 2.)
    plt.plot((wavetime,wavetime), (0,180), c = 'w', linestyle = style, linewidth = 2.)
    plt.xlim(starttime, endtime)
    ax.yaxis.set_visible(False)
    if i != len(E_kin)-1:
        ax.xaxis.set_visible(False)
    plt.xlim(starttime, endtime)
    cb1 = plt.colorbar()
    cb1.ax.set_ylabel('norm. flux')
    leg = ax.legend(loc = 'upper right', frameon=False, prop={'size': txtsize})
    for text in leg.get_texts():
	    plt.setp(text, color = cdi)
    for item in leg.legendHandles:
	    item.set_visible(False)
    print('finished plotting panel', panels_r[i]) 

# Now we are finalizing the figures  
print('setting the time axis labels')      
ax.set_xlim(starttime, endtime)
ax.set_xlabel('9 January 2014')
date_format = mdates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(date_format)


print('setting the Daa axis labels')
dx.set_xlabel(r'D$_{aa}$ (solid) $<$D$_{aa}>$ (dashed)')
print('saving the figure to a png')
plt.savefig('figures/Fig1.png')
plt.close()
print('Finished Figure 1')


#Now we are making Figure 2 with the average pitch angle 
#distribution during different periods of the compression
#First period is before the compression hits
#Second period is between the compression and the start of the wave observations
#Third is during the period when the wave is observed
index = np.where((time >= plttimestart ) & (time <= compression))
precomp_event = electrons_data[index, :, :]
panel = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)', 'k)', 'l)']
fig22 = plt.figure(22, figsize=(8.5,22))
for i in range(len(E_kin)):
	print('now making the average pitch angle distribution Figure 2 panel',  panel[i])
	dx = plt.subplot(len(E_kin)*3+1, 1, i+1)
	temp_PAD = 0.*(event[0,0,:,channels[i]])
	dx.plot(alpha*180./np.pi, daab_plt[i,:]/np.nanmax(daab_plt[i,:]), linestyle = style,  c = cli)
	dx.plot(180.- alpha*180./np.pi, (daab_plt[i,:]/np.nanmax(daab_plt[i,:])), linestyle = style,  c = cli)
	for j in range(len(precomp_event[0,:,0,0])):
		temp_PAD = (precomp_event[0,j,:,channels[i]]/(np.nanmax(event[0,:,:,channels[i]]))) + temp_PAD
	dx.plot(HOPE_alpha, temp_PAD/len(precomp_event[0,:,0,0]), c = cdi, label = panel[i])
	dx.set_ylabel(np.str(np.int(E_kin[i]*10**6))+ ' eV')
	dx.set_ylim(0,1)
	if i != len(E_kin)-1:
		dx.xaxis.set_visible(False)
	if i == 0:
		dx.set_title('Prior to compression')
	leg = dx.legend(loc = 'upper right', frameon=False, prop={'size': txtsize})
	for item in leg.legendHandles:
		item.set_visible(False)  
	dx.set_xlim(0, 180)


index = np.where((time >= compression ) & (time <= wave))
comp_prewave_event = electrons_data[index, :, :]


for i in range(len(E_kin)):
	print('now making the average pitch angle distribution Figure 2 panel',  panel[i+4])
	dx = plt.subplot(len(E_kin)*3+1, 1, i+5)
	temp_PAD = 0.*(comp_prewave_event[0,0,:,channels[i]])
	dx.plot(alpha*180./np.pi, daab_plt[i,:]/np.nanmax(daab_plt[i,:]), linestyle = style, c = cli)
	dx.plot(180.- alpha*180./np.pi, (daab_plt[i,:]/np.nanmax(daab_plt[i,:])), linestyle = style, c = cli)
	for j in range(len(comp_prewave_event[0,:,0,0])):
		temp_PAD = (comp_prewave_event[0,j,:,channels[i]]/(np.nanmax(event[0,:,:,channels[i]]))) + temp_PAD
	dx.plot(HOPE_alpha, temp_PAD/len(comp_prewave_event[0,:,0,0]), c = cdk, label = panel[i+4])
	dx.set_ylabel(np.str(np.int(E_kin[i]*10**6))+ ' eV')
	dx.set_ylim(0,1)
	if i != len(E_kin)-1:
		dx.xaxis.set_visible(False)
	if i == 0:
		dx.set_title('Compression to start of wave')
	leg = dx.legend(loc = 'upper right', frameon=False, prop={'size': txtsize})
	for item in leg.legendHandles:
		item.set_visible(False) 
	dx.set_xlim(0, 180)


index = np.where((time >= wave ) & (time <= wave_end))
wave_event = electrons_data[index, :, :]


for i in range(len(E_kin)):
	print('now making the average pitch angle distribution Figure 2 panel',  panel[i+8])
	dx = plt.subplot(len(E_kin)*3+1, 1, (i+9))
	temp_PAD = 0.*(wave_event[0,0,:,channels[i]])
	dx.plot(alpha*180./np.pi, daab_plt[i,:]/np.nanmax(daab_plt[i,:]), linestyle = style, c = cli)
	dx.plot(180.- alpha*180./np.pi, (daab_plt[i,:]/np.nanmax(daab_plt[i,:])), linestyle = style, c = cli)
	for j in range(len(wave_event[0,:,0,0])):
		temp_PAD = (wave_event[0,j,:,channels[i]]/(np.nanmax(event[0,:,:,channels[i]]))) + temp_PAD
	dx.plot(HOPE_alpha, temp_PAD/len(wave_event[0,:,0,0]), c = cdr, label = panel[i+8])
	dx.set_ylabel(np.str(np.int(E_kin[i]*10**6))+ ' eV')
	dx.set_ylim(0,1)
	if i != len(E_kin)-1:
		dx.xaxis.set_visible(False)
	if i == 0:
		dx.set_title('During Chorus wave')
	leg = dx.legend(loc = 'upper right', frameon=False, prop={'size': txtsize})
	for item in leg.legendHandles:
		item.set_visible(False) 
	dx.set_xlim(0, 180)

plt.tight_layout()
plt.xlabel(r'pitch angle')
plt.savefig('figures/Fig2.png')
plt.close()
# 

#Here we are making a second version of Figure 2. 

index = np.where((time >= plttimestart ) & (time <= compression))
precomp_event = electrons_data[index, :, :]
#panel = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)', 'k)', 'l)']
panel = ['pre-compression', 'compression to wave', 'during wave']
panel2 = ['a)', 'b)', 'c)', 'd)']
fig22 = plt.figure(22, figsize=(8.5,11))
for i in range(len(E_kin)):
	print('now making the average pitch angle distribution Figure 2 panel',  panel[0])
	dx = plt.subplot(len(E_kin)+1, 1, i+1)
	temp_PAD = 0.*(event[0,0,:,channels[i]])
	dx.plot(alpha*180./np.pi, daab_plt[i,:]/np.nanmax(daab_plt[i,:]), linestyle = style,  c = cli, label = r'$<D_{\alpha\alpha}>$')
	dx.plot(180.- alpha*180./np.pi, (daab_plt[i,:]/np.nanmax(daab_plt[i,:])), linestyle = style,  c = cli)
	for j in range(len(precomp_event[0,:,0,0])):
		temp_PAD = (precomp_event[0,j,:,channels[i]]/(np.nanmax(event[0,:,:,channels[i]]))) + temp_PAD
	dx.plot(HOPE_alpha, temp_PAD/len(precomp_event[0,:,0,0]), c = cdi, label = panel[0])
	dx.set_ylabel(np.str(np.int(E_kin[i]*10**6))+ ' eV')
	dx.set_ylim(0,1)
	dx.annotate(panel2[i], xy=(170, 0.8))
	if i != len(E_kin)-1:
		dx.xaxis.set_visible(False)


index = np.where((time >= compression ) & (time <= wave))
comp_prewave_event = electrons_data[index, :, :]


for i in range(len(E_kin)):
	print('now making the average pitch angle distribution Figure 2 panel',  panel[1])
	dx = plt.subplot(len(E_kin)+1, 1, i+1)
	temp_PAD = 0.*(comp_prewave_event[0,0,:,channels[i]])
	dx.plot(alpha*180./np.pi, daab_plt[i,:]/np.nanmax(daab_plt[i,:]), linestyle = style, c = cli)
	dx.plot(180.- alpha*180./np.pi, (daab_plt[i,:]/np.nanmax(daab_plt[i,:])), linestyle = style, c = cli)
	for j in range(len(comp_prewave_event[0,:,0,0])):
		temp_PAD = (comp_prewave_event[0,j,:,channels[i]]/(np.nanmax(event[0,:,:,channels[i]]))) + temp_PAD
	dx.plot(HOPE_alpha, temp_PAD/len(comp_prewave_event[0,:,0,0]), c = cdk, label = panel[1])


index = np.where((time >= wave ) & (time <= wave_end))
wave_event = electrons_data[index, :, :]


for i in range(len(E_kin)):
	print('now making the average pitch angle distribution Figure 2 panel',  panel[2])
	dx = plt.subplot(len(E_kin)+1, 1, (i+1))
	temp_PAD = 0.*(wave_event[0,0,:,channels[i]])
	dx.plot(alpha*180./np.pi, daab_plt[i,:]/np.nanmax(daab_plt[i,:]), linestyle = style, c = cli)
	dx.plot(180.- alpha*180./np.pi, (daab_plt[i,:]/np.nanmax(daab_plt[i,:])), linestyle = style, c = cli)
	for j in range(len(wave_event[0,:,0,0])):
		temp_PAD = (wave_event[0,j,:,channels[i]]/(np.nanmax(event[0,:,:,channels[i]]))) + temp_PAD
	dx.plot(HOPE_alpha, temp_PAD/len(wave_event[0,:,0,0]), c = cdr, label = panel[2])
	leg = dx.legend(loc = 'upper left', frameon=False, prop={'size': txtsize})
	dx.set_xlim(0, 180)
	

plt.tight_layout()
plt.xlabel(r'pitch angle')
plt.savefig('figures/Fig2b.png')
plt.close()

