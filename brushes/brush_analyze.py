import os
from sys import argv
import math
import numpy as np
import matplotlib
from pylab import *

script, polymer_file, ctr_file, salt_file = argv

#analyzes the polymers.trj file and compute mean brush height and dist function.  Coordinates must be wrapped.

def round_partial(value,res):
    return round( (value / res) * res)

tot_sample = 0
tot_atom = 0
xmin = 0
xmax = 15
ymin = 0
ymax = 15
zmin = 0
zmax = 62
grid_size = 1.0

#ANALYZE POLYMER TRAJECTORY FILE.  POLY_WRAP
poly_grid = np.zeros ( (int((xmax-xmin)/grid_size),int((ymax-ymin)/grid_size),int((zmax-zmin)/grid_size)) , dtype = float)

with open(polymer_file,'r') as f:
    for i,line in enumerate(f):
        templine = line.split()
        if len(templine) == 2 and templine[1] == "TIMESTEP":
            tot_sample += 1
        if i == 3:
            tot_atom = int(line)
        elif templine[0] != "ITEM:" and len(templine) == 6:
            xpos = int(round_partial(float(templine[3]),grid_size)/grid_size)
            ypos = int(round_partial(float(templine[4]),grid_size)/grid_size)
            zpos = int(round_partial(float(templine[5]),grid_size)/grid_size)
            poly_grid[xpos-1,ypos-1,zpos-1] += 1
        else:
            continue

#average over all atoms and snapshots
poly_grid = poly_grid/(float(tot_atom)*float(tot_sample))

#Compute z-distribution function of polymer density
z_dist = np.mean(poly_grid,axis=(0,1))
z_arr = np.arange(grid_size,zmax-zmin+grid_size,grid_size)

z_dist = np.insert(z_dist,0.0,0)
z_arr = np.insert(z_arr,0.0,0)

#normalize z_dist
z_dist = z_dist/(sum(z_dist))
#Compute mean brush height
z_mean = np.dot(z_dist,z_arr)/sum(z_dist)
print "The mean brush height is "+str(z_mean)+" sigma."

#compute xz face
xz_dist = np.mean(poly_grid,axis=0)
xz_dist = xz_dist/sum(xz_dist)
xz_dist = xz_dist.T

#compute yz face
yz_dist = np.mean(poly_grid,axis=1)
yz_dist = yz_dist/sum(yz_dist)
yz_dist = yz_dist.T

side_view = (xz_dist + yz_dist)/2.

#Plot z particle density distribution function
plt.figure(1)
plt.subplot(221)
plt.plot(z_dist,z_arr,color='black',linewidth=2)
plt.ylim(zmin,zmax/2)
#For probability axis limits, can use max(z_dist), but practically speaking 0.2 is good
plt.xlim(0,0.2)
plt.title('Polymer monomer density distribution',fontsize=12)
plt.ylabel('z-axis')
plt.xlabel('Probability')

#Plot side-view of polymer density
plt.subplot(222)
plt.imshow(side_view,origin='lower')
plt.title('Side face view of polymer density distribution.',fontsize=12)
plt.xlabel('<x/y> axis')
plt.ylabel('z axis')
plt.ylim(zmin,zmax/2.)

#Analyze Counterion Distribution 
ctr_grid = np.zeros ( (2,int((xmax-xmin)/grid_size),int((ymax-ymin)/grid_size),int((zmax-zmin)/grid_size)) , dtype = float)

with open(ctr_file,'r') as f:
    for i,line in enumerate(f):
        templine = line.split()
        if len(templine) == 2 and templine[1] == "TIMESTEP":
            tot_sample += 1
        if i == 3:
            tot_atom = int(line)
        elif templine[0] != "ITEM:" and len(templine) == 6:
            xpos = int(round_partial(float(templine[3]),grid_size)/grid_size)
            ypos = int(round_partial(float(templine[4]),grid_size)/grid_size)
            zpos = int(round_partial(float(templine[5]),grid_size)/grid_size)
            if templine[2] == "-1":
                ctr_grid[0,xpos-1,ypos-1,zpos-1] += 1
            elif templine[2] == "1":
                ctr_grid[1,xpos-1,ypos-1,zpos-1] += 1
            else:
                print "You have entered an invalid charge type for analysis."
        else:
            continue

#average over all atoms and snapshots
for i in range(2):
    if sum(ctr_grid[i]) > 0.0:
        ctr_grid[i] = ctr_grid[i]/(sum(ctr_grid[i])*float(tot_sample))

#Compute z-distribution function of ctrion density.  Implicity 0=-1,1=+1.
z_dist_arr = []
z_dist_arr.append(np.mean(ctr_grid[0],axis=(0,1)))
z_dist_arr.append(np.mean(ctr_grid[1],axis=(0,1)))
z_arr = np.arange(grid_size,zmax-zmin+grid_size,grid_size)

z_dist_arr[0] = np.insert(z_dist_arr[0],0.0,0)
z_dist_arr[1] = np.insert(z_dist_arr[1],0.0,0)

z_arr = np.insert(z_arr,0.0,0)

#normalize z_dist and compute mean ctr ion heights
z_mean_arr = []
for i in range(len(z_dist_arr)):
    if sum(z_dist_arr[i]) > 0.0:
        z_dist_arr[i] = z_dist_arr[i]/sum(z_dist_arr[i])
        z_mean_arr.append(np.dot(z_dist_arr[i],z_arr)/sum(z_dist_arr[i]))
    else:
        z_mean_arr.append(0.0)

color_list = ['blue','red']
legend_list = ['-1 anion','+1 cation']

#Plot z-distribution of counter ion density

plt.subplot(223)
for i in range(len(z_mean_arr)):
    if z_mean_arr[i] > 0.0:
        plt.plot(z_dist_arr[i],z_arr,color=color_list[i],linewidth=2,label=legend_list[i])
        plt.ylim(zmin,zmax/2.)
        plt.xlim(0,0.2)
    plt.title('Counter ion density distribution',fontsize=12)
    plt.ylabel('z-axis')
    plt.xlabel('Probability')
    plt.legend()

#Analyze Salt Distribution 
#salt_grid has an extra dimension to accomodate different salt types (-1,+1,+2,+3)
salt_grid = np.zeros ( (4,int((xmax-xmin)/grid_size),int((ymax-ymin)/grid_size),int((zmax-zmin)/grid_size)) , dtype = float)

with open(salt_file,'r') as f:
    for i,line in enumerate(f):
        templine = line.split()
        if len(templine) == 2 and templine[1] == "TIMESTEP":
            tot_sample += 1
        if i == 3:
            tot_atom = int(line)
        elif templine[0] != "ITEM:" and len(templine) == 6:
            xpos = int(round_partial(float(templine[3]),grid_size)/grid_size)
            ypos = int(round_partial(float(templine[4]),grid_size)/grid_size)
            zpos = int(round_partial(float(templine[5]),grid_size)/grid_size)
            if templine[2] == "-1":
                salt_grid[0,xpos-1,ypos-1,zpos-1] += 1
            elif templine[2] == "1":
                salt_grid[1,xpos-1,ypos-1,zpos-1] += 1
            elif templine[2] == "2":
                salt_grid[2,xpos-1,ypos-1,zpos-1] += 1
            elif templine[2] == "3":
                salt_grid[3,xpos-1,ypos-1,zpos-1] += 1
            else:
                print "Atom type with charge "+str(templine[1])+" not included in current version."
        else:
            continue


#average over all atoms and snapshots
for i in range(4):
    if sum(salt_grid[i]) > 0.0:
        salt_grid[i] = salt_grid[i]/(sum(salt_grid[i])*float(tot_sample))

#Compute z-distribution function of salt densities
#z_dist_arr corresponds implicitly to 0=-1,1=+1,2=+2,3=+3
z_dist_arr = []

z_dist_arr.append(np.mean(salt_grid[0],axis=(0,1)))
z_dist_arr.append(np.mean(salt_grid[1],axis=(0,1)))
z_dist_arr.append(np.mean(salt_grid[2],axis=(0,1)))
z_dist_arr.append(np.mean(salt_grid[3],axis=(0,1)))

z_arr = np.arange(grid_size,zmax-zmin+grid_size,grid_size)
z_arr = np.insert(z_arr,0.0,0)

z_dist_arr[0] = np.insert(z_dist_arr[0],0.0,0)
z_dist_arr[1] = np.insert(z_dist_arr[1],0.0,0)
z_dist_arr[2] = np.insert(z_dist_arr[2],0.0,0)
z_dist_arr[3] = np.insert(z_dist_arr[3],0.0,0)

#normalize z_dist and computer mean salt heights
z_mean_arr = []

for i in range(4):
    if sum(z_dist_arr[i]) > 0.0:
        z_dist_arr[i] = z_dist_arr[i]/(sum(z_dist_arr[i]))
        z_mean_arr.append(np.dot(z_dist_arr[i],z_arr)/sum(z_dist_arr[i]))
    else:
        z_mean_arr.append(0.0)

color_list = ['orange','purple','green','cyan']
legend_list = ['-1 anion','+1 cation','+2 cation','+3 cation']

plt.subplot(224)
for i in range(len(z_mean_arr)):
    if z_mean_arr[i] > 0.0:
        plt.plot(z_dist_arr[i],z_arr,color=color_list[i],linewidth=2,label=legend_list[i])
        plt.ylim(zmin,zmax/2.)
        plt.xlim(0,0.2)
    plt.title('Salt ion density distribution',fontsize=12)
    plt.ylabel('z-axis')
    plt.xlabel('Probability')
    plt.legend()

#plt.tight_layout()
plt.show()

#analyze the lammps.trj file and computes the number density of polymers, salt, and cations as a function of z position
