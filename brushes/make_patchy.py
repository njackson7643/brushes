import os
from functions import *
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from math import sqrt
import random
import scipy
from scipy.spatial import distance    

def find_neigh(dist_mat,part_list,thresh,seed):
    neigh_list = []
    for item,i in zip(dist_mat,range(len(dist_mat))):
        if item[0] < thresh:
            neigh_list.append(part_list[i])
    return neigh_list

def diff(first,second):
    return list(set(first) - set(second))
num_part = 0
part_list = []
tot_chg = -8

#Load in icosahedron coordinates into part_list
f = open('QET_data.xyz','r')
ico_lines = f.readlines()
for i,line in zip(range(len(ico_lines)),ico_lines):
    templine = line.split()
    if i == 0:
        num_part = int(templine[0])
    if len(templine) > 0 and i != 0:
        part_list.append((float(templine[1]),float(templine[2]),float(templine[3])))
    else: 
        continue

#Pick a particle at random from part_list and create a section of charged sites around it
rand1 = random.randint(0,num_part)
seed_arr = np.array([part_list[rand1]])
part_arr = np.array(part_list) 
dist_check = scipy.spatial.distance.cdist(part_arr,seed_arr)
patch1 = find_neigh(dist_check,part_list,3.0,rand1)

rest_list = diff(part_list,patch1)

num_patch_sites = len(patch1)

num_opp_sites = num_patch_sites + tot_chg

print str(len(patch1))+" negative sites identified in the patch"
print "Adding "+str(num_opp_sites)+" positive sites randomly to the particle"
