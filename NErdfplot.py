#!/usr/bin/python

# Import module for plot "matplotlib"
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import rc
import matplotlib
import matplotlib.ticker as mtick
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib import rcParams
from pylab import *
import pylab as pylab
import matplotlib.patches as patches
from matplotlib.pyplot import *

# Read data file
k1 = np.loadtxt("random_ISG-NE_RDF.xvg")

# Specifiy environmental parameter
rc('font',**{'family':'serif','serif':['Arial']})

# Create axes 
fig = plt.figure(figsize=(8.5,6)) #cm
fig.subplots_adjust(bottom=0.15)
fig.subplots_adjust(left=0.18)

# Main figure
ax1 = plt.subplot2grid((1,1), (0, 0))

#ax1.set_title("Plot title...")    
ax1.set_xlabel('$r$ (${\mu}m$)',fontname="Arial",fontweight="normal",fontsize="20")
ax1.set_ylabel('g($r$)',fontname="Arial",fontweight="normal",fontsize="20")
ax1.tick_params(direction='in', pad=6)
ax1.xaxis.set_label_coords(0.5, -0.1)
ax1.yaxis.set_label_coords(-0.1, 0.5)
ax1.set_xlim([0,0.2])
ax1.set_ylim([0,3])
xmajorLocator   = MultipleLocator(0.1)
xmajorFormatter = FormatStrFormatter('%.1f')
xminorLocator   = MultipleLocator(0.05)
ax1.xaxis.set_major_locator(xmajorLocator)
ax1.xaxis.set_major_formatter(xmajorFormatter)
ax1.xaxis.set_minor_locator(xminorLocator)
ymajorLocator   = MultipleLocator(0.5)
ymajorFormatter = FormatStrFormatter('%.1f')
yminorLocator   = MultipleLocator(0.25)
ax1.yaxis.set_major_locator(ymajorLocator)
ax1.yaxis.set_major_formatter(ymajorFormatter)
ax1.yaxis.set_minor_locator(yminorLocator)
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(20)
    tick.label.set_fontname("Arial")
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(20)
    tick.label.set_fontname("Arial")
for axis in ['bottom','left']:ax1.spines[axis].set_linewidth(2)
for axis in ['top','right']:ax1.spines[axis].set_linewidth(0)
for line in ax1.xaxis.get_ticklines():
    line.set_markersize(5)
    line.set_markeredgewidth(2)
for line in ax1.yaxis.get_ticklines():
    line.set_markersize(5)
    line.set_markeredgewidth(2)
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'
for line in ax1.yaxis.get_minorticklines():
    line.set_markersize(2.5)
    line.set_markeredgewidth(2)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
#ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
#ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
#xtick_locs=[0.00,0.03,0.05,0.10,0.20]
#xtick_lbls=['0.00','0.03',' 0.05','0.10','0.20']
#plt.xticks(xtick_locs,xtick_lbls)
#ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.3f'))
#ax1.set_yticks([0.015,0.025,0.035,0.045,0.055])
#ax1.set_axis_bgcolor('none')
#ax1.grid(True)
ax1.plot((0,0.6),(1,1),'k',linestyle="-",linewidth=1.5)
#ax1.plot((0,1000),(5.4,5.4),'grey',linestyle=":",linewidth=1.5)
#ax1.plot((0,1000),(5.6,5.6),'grey',linestyle=":",linewidth=1.5)
#ax1.plot((0,1000),(5.8,5.8),'grey',linestyle=":",linewidth=1.5)

#12:278

# Plot 
ax1.plot(k1[:,0]/100000,k1[:,1],linestyle='-',c='k',linewidth=1,alpha=1)
#ax1.plot((0,np.mean(k2[13:278,1]))(0.6,np.mean(k2[13:278,1])),linestyle='-',c='r',linewidth=2)
#ax1.errorbar(k2[:,0],k2[:,1], yerr=k2[:,2],marker='o',linestyle='none',markersize=6, fmt='-',capsize=4, elinewidth=2,linewidth=2,c='green')
#ax1.plot(k2[:,0],k2[:,1],marker='o',linestyle='none',c='green',markersize=8)

# build a rectangle in axes coords
left, width = .25, .5
bottom, height = .25, .5
right = left + width
top = bottom + height
ax1.text(0.04*(left+right), 0.96*(bottom+top), 'NE-ISG',horizontalalignment='left',verticalalignment='center',fontsize=20,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)
#ax1.text(0.02*(left+right), 0.89*(bottom+top), 'anodic layer',horizontalalignment='left',verticalalignment='center',fontsize=20,fontname="Arial",fontweight="normal",color='chocolate',transform=ax1.transAxes)
#ax1.text(0.02*(left+right), 0.83*(bottom+top), 'bilayer',horizontalalignment='left',verticalalignment='center',fontsize=20,fontname="Arial",fontweight="normal",color='k',transform=ax1.transAxes)

plt.show()

# Save figures, (PNG,JPG,EPS,SVG,PGF and PDF supported, PDF is recommanded or PGF)
fig.savefig("random_ISG-NE_RDF.eps",dpi=1200)
