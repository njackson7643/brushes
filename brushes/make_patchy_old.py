import os
from functions import *
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from math import sqrt
import random
import scipy
from scipy.spatial import distance    

#def screen_dmat(dist_check,thresh):
#    if all(item >= thresh for item in dist_check):
#        return True
#    else:
#        return False

#def check_neigh(part_list,xt,yt,zt,thresh):
#    part_arr = np.array(part_list)
#    test_arr = np.array([[xt,yt,zt]])
#    dist_check = scipy.spatial.distance.cdist(part_arr,test_arr)
#    for item in dist_check:
#        if item[0] < thresh:
#            return False
#    else:
#        return True

#radius = 2.0 #in terms of LJ sigma
#part_size = 0.1 #radius of individual beads that make up the particle

#num_part = int((radius**2)/(part_size)**2)

#part_count = 0

#part_list = [[0.0,0.0,0.0]]
#theta_list = []
#phi_list = []

#while part_count < num_part:
#    rand1 = random.random()
#    rand2 = random.random()
#    x = 2.*radius*math.sqrt(rand1*(1-rand1))*math.cos(2.*math.pi*rand2)
#    y = 2.*radius*math.sqrt(rand1*(1-rand1))*math.sin(2.*math.pi*rand2)
#    z = radius*(1.-2.*rand1)
#    if check_neigh(part_list,x,y,z,2.*part_size) == True:
#        part_list.append([x,y,z])
#        part_count += 1
#    else:
#        continue


#fig = plt.figure()
#ax = fig.gca(projection='3d')
#for i in range(len(part_list)):
#    ax.scatter(part_list[i][0],part_list[i][1],part_list[i][2],linewidth=2.0)

#ax.set_xlim(-4.5,4.5)
#ax.set_ylim(-4.5,4.5)
#ax.set_zlim(-4.5,4.5)
#plt.show()

#make patches
#rand_int = random.randint(0,len(part_list))
#rand_site = poly_list[rand_int]
#part_arr = np.array(part_list)
#test_arr = np.array([rand_site])
#dist_check = scipy.spatial.distance.cdist(part_arr,test_arr)

