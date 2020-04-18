# -*- Resonant Chamber-*-
# -*- coding: utf-8 -*-

# transmission around a 90-degree waveguide bend in 2d
from __future__ import division

import meep as mp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import h5py
from scipy import constants as S

####################################### MEEP SIMULATION #####################################################
#############################################################################################################
# - Constants
pi = S.pi
eps0 = S.epsilon_0
c = S.c
mu0 = S.mu_0

res = 5 # pixels/cm
sx = 38.  # size of cell in X direction
sy = 17.  # size of cell in Y direction
sz = 2.6  # size of cell in Y direction
dpml = 0.2
r = 5.
a = 0.1 #Meep unit (meters)
Wt = mp.Medium(epsilon=80, D_conductivity=0.2 * a/(c * 80 * eps0)) # Water dielectric properties (background) 
At = mp.Medium(epsilon=-1e20, D_conductivity=9.5 * a/(c * 80 * eps0)) # Antennas dielectric properties
fcen = 800e6*(a/c)  # pulse center frequency
df = 50e6*(a/c)     # pulse width (in frequency)
#Meep variable

#Domain
cell = mp.Vector3(sx,sy,sz)
pml_layers = [mp.PML(dpml)]
geometry = [mp.Block(mp.Vector3(sx,sy,mp.inf), center=mp.Vector3(),material=Wt), 
            mp.Block(mp.Vector3(sx,sy,0.2), center=mp.Vector3(0,0,0.6),material=mp.metal), #upper plate
            mp.Block(mp.Vector3(sx,sy,0.2), center=mp.Vector3(0,0,-0.6),material=mp.metal), #lower plate
#            mp.Block(mp.Vector3(0.2,0.2,0.8), center=mp.Vector3(-14.,5.,0),material=mp.metal), #Antenna 1
            mp.Block(mp.Vector3(0.2,0.2,0.8), center=mp.Vector3(-14.,-5.,0),material=mp.metal), #Antenna 2
            mp.Block(mp.Vector3(0.2,0.2,0.8), center=mp.Vector3(10.,5.,0),material=mp.metal), #Antenna 3
            mp.Block(mp.Vector3(0.2,0.2,0.8), center=mp.Vector3(10.,-5.,0),material=mp.metal), #Antenna 4
            mp.Cylinder(material=Wt, radius=r, height=sz, center=mp.Vector3())] #hole in the middle 

sources = [mp.Source(mp.GaussianSource(fcen,fwidth=df), component=mp.Ez, center=mp.Vector3(-14,5.,0.8))]
#sources = [mp.Source(mp.ContinuousSource(frequency=fcen, width=20), component=mp.Ez, center=mp.Vector3(-14,5.))]
sim = mp.Simulation(cell_size=cell, geometry=geometry, boundary_layers=pml_layers, sources=sources, resolution=res)#boundary_layers=pml_layers

#sim.run(until_after_sources=mp.stop_when_fields_decayed(1000,mp.Ez,mp.Vector3(-14,5.,),1e-3))
sim.run(mp.at_beginning(mp.output_epsilon),
        mp.in_volume(mp.Volume(center=mp.Vector3(0,0,0), size=mp.Vector3(sx,sy,0)), mp.to_appended("ez", mp.at_every(0.6, mp.output_efield_z))),
        mp.in_volume(mp.Volume(center=mp.Vector3(-14 ,5.,0), size=mp.Vector3(0,0,0)), mp.to_appended("EzTr", mp.at_every(0.6, mp.output_efield_z))),
        mp.in_volume(mp.Volume(center=mp.Vector3(-14.,-5.,0), size=mp.Vector3(0,0,0)), mp.to_appended("EzR2", mp.at_every(0.6, mp.output_efield_z))),
        mp.in_volume(mp.Volume(center=mp.Vector3(10.,5.,0), size=mp.Vector3(0,0,0)), mp.to_appended("EzR3", mp.at_every(0.6, mp.output_efield_z))),
        mp.in_volume(mp.Volume(center=mp.Vector3(10.,-5.,0), size=mp.Vector3(0,0,0)), mp.to_appended("EzR4", mp.at_every(0.6, mp.output_efield_z))),
        until=2000)
#############################################################################################################
#############################################################################################################

########### LOAD ARRAYS #############
problemname = "resonant_chamber"

#Tensor with electric field in time & space
filename = problemname+"-ez.h5"
f = h5py.File(filename, 'r')
a_group_key = list(f.keys())[0]
ezT = np.asarray(list(f[a_group_key]))

#Permittivity Tensor
filename = problemname+"-eps-000000.00.h5"
f = h5py.File(filename, 'r')
a_group_key = list(f.keys())[0]
eps = np.asarray(list(f[a_group_key]))

#Transmitter antenna field
filename = problemname+"-EzTr.h5"
f = h5py.File(filename, 'r')
a_group_key = list(f.keys())[0]
ezt = np.asarray(list(f[a_group_key]))

#Receiver2 antenna field
filename = problemname+"-EzR2.h5"
f = h5py.File(filename, 'r')
a_group_key = list(f.keys())[0]
ezr2 = np.asarray(list(f[a_group_key]))

#Receiver2 antenna field
filename = problemname+"-EzR3.h5"
f = h5py.File(filename, 'r')
a_group_key = list(f.keys())[0]
ezr3 = np.asarray(list(f[a_group_key]))

#Receiver4 antenna field
filename = problemname+"-EzR4.h5"
f = h5py.File(filename, 'r')
a_group_key = list(f.keys())[0]
ezr4 = np.asarray(list(f[a_group_key]))
#######################################


############## PLOTS ##################
#Plot Antennas data
res = 5.
courant = 0.5 #default in meep
a = 0.01
dt = a/S.c#(res/courant)*a/S.c
t = np.arange(0,len(ezt)*dt,dt)

fig, axs = plt.subplots(4, 1)
#axs.plot(t/1e-9, ezt*0.01, '-b', label = '|V$_{transmitted}$|')
#axs.plot(t/1e-9, ezr*0.01, '--r', label = '|V$_{received}$|')
#axs.set_xlabel('time (ns)')
#axs.grid(True)
axs[0].plot(t/1e-9, ezt*0.01, '-b', label = '|V$_{transmitted}$|')
axs[1].plot(t/1e-9, ezr2*0.01, '-r', label = '|V2$_{received}$|')
axs[2].plot(t/1e-9, ezr3*0.01, '-r', label = '|V3$_{received}$|')
axs[3].plot(t/1e-9, ezr4*0.01, '-r', label = '|V4$_{received}$|')
axs[1].set_xlabel('time (ns)')
axs[0].grid(True)
axs[1].grid(True)
axs[2].grid(True)
axs[3].grid(True)
plt.legend(loc='lower right', shadow='False')
plt.xlabel('time')
plt.ylabel('Ez')
plt.show()

#Plot Permittivity
i = 8 #z index
j = 42 #x index
plt.imshow(eps[:,:,i], cmap = 'seismic', animated=True)
plt.show()
plt.imshow(eps[:,j,:], cmap = 'seismic', animated=True)
plt.show()

#Permittivity 3D view
from mayavi import mlab
s = mlab.contour3d(eps, colormap="seismic")
mlab.show()

#Plot Animated field
fig = plt.figure()
i = 0
imax = len(ezT[0,0,:])
z = 2
im = plt.imshow(ezT[:,:,i], cmap = 'seismic', animated=True)
def animate(frame):
    im.set_data(ezT[:,:,frame])
    return im,
ani = animation.FuncAnimation(fig, animate, frames=imax, interval=50, repeat_delay=50, blit=True)
plt.show()
#######################################

################ FFT ################## 
#*****************************************
#t = t * 1e-9
vs = ezt
vr2 = ezr2
vr3 = ezr3
vr4 = ezr4

#Sampling Parameters
N = len(t)
tf = max(t)

deltat = tf / N #sampling interval
fs = 1 / deltat #sampling freq
x = np.linspace(0.0, fs, N)
fr = fs/len(t) #freq resolution

fft = np.fft.fft(vs)
fft2 = np.fft.fft(vr2)
fft3 = np.fft.fft(vr3)
fft4 = np.fft.fft(vr4)


fig, axs = plt.subplots()
axs.plot(x[0: N//2], 20*np.log(np.abs(fft2)[0:N // 2] * (1 / N) / np.abs(fft)[0:N // 2]), '-b', label = '|S$_{12}$|') # 
axs.grid(True)
axs.set_xlabel('freq (Hz)')
axs.grid(True)

axs.plot(x[0: N//2], 20*np.log(np.abs(fft3)[0:N // 2] * (1 / N) / np.abs(fft)[0:N // 2]), '-r', label = '|S$_{13}$|') # 
axs.grid(True)
axs.set_xlabel('freq (Hz)')
axs.grid(True)

axs.plot(x[0: N//2], 20*np.log(np.abs(fft4)[0:N // 2] * (1 / N) / np.abs(fft)[0:N // 2]), '-g', label = '|S$_{14}$|') # 
axs.grid(True)
axs.set_xlabel('freq (Hz)')
axs.grid(True)


plt.legend(loc='lower right', shadow='False')
plt.xlabel('Freq (GHz)')
plt.ylabel('dB')
plt.show()
#######################################


