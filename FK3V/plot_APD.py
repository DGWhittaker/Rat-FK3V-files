#This script can be used to plot the rat FK3V APDR curve
#Author: Dominic Whittaker

import numpy as np
import matplotlib.pyplot as p
from matplotlib import * 

p.rcdefaults()

p.rc('lines', markeredgewidth=2)
p.rc('lines', markersize=7)
p.rc('xtick.major', size=5)
p.rc('ytick.major', size=5) #changes size of y axis ticks on both sides
p.rc('xtick', direction='out')
p.rc('ytick', direction='out')

a = np.loadtxt( 'FK3V_APD_1D_LV.dat', unpack=True )
b = np.loadtxt( 'Benoist_exp_2012_APD.dat', unpack=True )

FK_DI = a[0]
FK_APD = a[2]

exp_DI = b[0]
exp_APD = b[1]

fig = p.figure(1, figsize=(3,2.25)) #(8,6) seems to be the default
fig.set_facecolor('white') #this changes the background colour from the default gray/blue to white

ax1 = fig.add_subplot(111)
ax1.plot( FK_DI, FK_APD, color='black', linewidth=1.5 )
ax1.scatter( exp_DI, exp_APD, s = 50, edgecolor='black',color='white', linewidth=1.5 )
ax1.set_xlabel( 'DI (ms)' )
ax1.set_ylabel( 'APD (ms)' )
ax1.axis([0, 160, 44, 53])
ax1.spines['right'].set_visible(False) #removes top and right borders
ax1.spines['top'].set_visible(False)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
p.xticks(np.arange(0,161,40))
p.yticks(np.arange(44,54+1,3))
p.legend(['Model', 'Exp.'], loc='best', frameon=False) #add frameon=False to remove the legend box
p.tight_layout()

p.show()

p.close()