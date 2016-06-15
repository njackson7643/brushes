import os
from functions_patchy import *
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from math import sqrt
import random
import scipy
from scipy.spatial import distance    

#This script places a patchy particle model of a protein above a polymer brush in a LAMMPS file

def find_neigh(dist_mat,part_list,thresh,seed):
    neigh_list = []
    for item,i in zip(dist_mat,range(len(dist_mat))):
        if item[0] < thresh:
            neigh_list.append(part_list[i])
    return neigh_list

def diff(first,second):
    return list(set(first) - set(second))

#Parse the setup file (brush.setup)                                                                     
input_param = setup_parser("brush.setup")
#Write .data file                                                                                       
num_atom = 0
num_bond = 0
num_ang = 0
num_dih = 0
num_imp = 0
num_atom_typ = 0
num_bond_typ = 0
num_ang_typ = 0
num_dih_typ = 0
num_imp_typ = 0

#MAKE TOP AND BOTTOM SURFACES                                                                           
substr_len = int(input_param['substr_len'])
lat_spacing = float(input_param['lat_spacing'])
chain_len = int(input_param['chain_len'])
poly_bond_len = float(input_param['poly_bond_len'])
grid_disc = 1.0 #each grid spacing corresponds to 0.1 LJ sigma                                          
inv_grid_disc = 1.0/grid_disc
part_per_1D = int(substr_len) #Number of substrate sites per Lx                                         
abs_len = int(substr_len*lat_spacing) #Absolute length of Lx/Ly axes of simulation grid in LJ sigma     
Lx = int(abs_len/grid_disc) #Size of Lx in grid points                                                  
Ly = Lx #Size of Ly in grid points                                                                      
#Total number of particles to place in bottom and top substrates                                        
tot_sub_part = part_per_1D**2
#Position of top substrate in z-direction                                                               
top_bound = int(2*(chain_len+1)*poly_bond_len/grid_disc) #Position of top substrate in grid points      
Lz = top_bound*3 #Total simulation box heigh in grid points

#Sim grid divides the simulation box into cubic elements defaulted to grid_disc LJ sigma discretization
sim_grid = np.zeros( (Lx,Ly,top_bound) , dtype = str)
#Place substrate particles 'S' at substrate points in space                                             
#Create dictionaries for holding SUBSTRATE charge, mass, LJ, atom number information                    
S_dict = {}
S_count = 0
for i in range(0,Lx,int(lat_spacing/grid_disc)):
    for j in range(0,Ly,int(lat_spacing/grid_disc)):
        sim_grid[i,j,0] = "S"
        S_count += 1
#Dict entry is absolute atom number, dict value is the corresponding i,j,k in sim_grid                  
        S_dict[S_count] = str(i)+','+str(j)+',0'
for i in range(0,Lx,int(lat_spacing/grid_disc)):
    for j in range(0,Ly,int(lat_spacing/grid_disc)):
        sim_grid[i,j,top_bound-1] = "S"
        S_count += 1
        S_dict[S_count] = str(i)+','+str(j)+','+str(top_bound-1)

#MAKE GRAFTED POLYMER CHAINS.  ASSUMES UNIFORM GRAFTING DENSITY ON CUBIC LATTICE.                       
##1.Determine the substrate sites to be grafted to                                                      
##Assuming the bjerrum length is that of water (~0.7 nm) then we can say bjerrum length = sigma, and sigma ~ 0.7 nm                                                                                            
num_chain = int(input_param['num_chain'])
abs_graft_dens_2D = float(num_chain/(substr_len*lat_spacing*0.71)**2)
abs_graft_dens_1D = math.sqrt(abs_graft_dens_2D)
grid_graft_dens_2D = float(num_chain/(substr_len*lat_spacing/grid_disc)**2)
grid_graft_dens_1D = math.sqrt(grid_graft_dens_2D)
chain_sep_1D = int(1.0/grid_graft_dens_1D)

print "The 1D grafting density is: "+str(grid_graft_dens_1D)+"\n"
print "The 1D chain separation is: "+str(chain_sep_1D)+"\n"

print "There are "+str(num_chain)+" chains spread over "+str(substr_len**2)+" substrate sites."
print "This corresponds to 1 chain every "+str(chain_sep_1D)+" grid sites, and an absolute"
print "grafting density of "+str(round(abs_graft_dens_2D,3))+" sites/(nm)**2"

##2.Determine the repeat structure of the polymer from 'chain_typ'                                      
chain_typ = input_param['chain_typ']
chain_list = chain_typ.split(',')
seq_len = len(chain_list)

##3.Place polymer chains on the lattice                                                                 
branch_choice = input_param['branch']
poly_len = int(chain_len*poly_bond_len/grid_disc) #size of extended polymer chain in grid points        
bond_size = int(poly_bond_len/grid_disc) #size of polymer bond in grid points                           
offset = 0
seq_count = 0
N_count = 0
P_count = 0
Z_count = 0
#Don't mess with bond_scale                                                                             
bond_scale = 2
side_bond_size = bond_size/bond_scale

##3a.Place purely linear polymers with the chain sequence purely along the backbone                     
#Create dictionaries for holding P & N charge, mass, LJ, atom number information                        
P_dict = {}
N_dict = {}
Z_dict = {}

print "Building brush..."

if branch_choice == "no":
#offset = int(substr_len%(int(math.sqrt(num_chain))))                                                   
    for i in range(offset,Lx,chain_sep_1D):
        seq_count = 0
        for j in range(offset,Ly,chain_sep_1D):
            seq_count = 0
            for k in range(bond_size,poly_len+1,bond_size):
                sim_grid[i,j,k] = str(chain_list[seq_count%seq_len])
                seq_count += 1
                if sim_grid[i,j,k] == "P":
                    P_count += 1
                    P_dict[P_count] = str(i)+','+str(j)+','+str(k)
                elif sim_grid[i,j,k] == "N":
                    N_count += 1
                    N_dict[N_count] = str(i)+','+str(j)+','+str(k)
                elif sim_grid[i,j,k] == "Z":
                    Z_count += 1
                    Z_dict[Z_count] = str(i)+','+str(j)+','+str(k)

##3b.Place branched polymer with neutral backbone and the sequence on the sidechain                                                    
elif branch_choice == 'yes':
    branch_rep = input_param['branch_rep']
    branch_alt = input_param['branch_alt']
    alt_track = 0
    for i in range(offset,Lx,chain_sep_1D):
        seq_count = 0
        for j in range(offset,Ly,chain_sep_1D):
            seq_count = 0
            for k in range(bond_size,poly_len+1,bond_size):
                sim_grid[i,j,k] = "Z"
                Z_count += 1
                Z_dict[Z_count] = str(i)+','+str(j)+','+str(k)
                seq_count += 1
                #The alternating branching choice is totally fucked.  To make it sort of work, I could instead just project the alternating in +x then +y then +x.  Currently, +x then -x has some serious issues.                                                           
                if seq_count%int(branch_rep) == 0:
                    if branch_alt == 'yes':
                        if alt_track%2 == 0:
                            for ib in range(seq_len):
                                sim_grid[i+side_bond_size*(ib+1),j,k] = str(chain_list[ib])
                                if sim_grid[i+side_bond_size*(ib+1),j,k] == "P":
                                    P_count += 1
                                    P_dict[P_count] = str(i+side_bond_size*(ib+1))+','+str(j)+','+str(k)
                                elif sim_grid[i+side_bond_size*(ib+1),j,k] == "N":
                                    N_count += 1
                                    N_dict[N_count] = str(i+side_bond_size*(ib+1))+','+str(j)+','+str(k)
                                elif sim_grid[i+side_bond_size*(ib+1),j,k] == "Z":
                                    Z_count += 1
                                    Z_dict[Z_count] = str(i+side_bond_size*(ib+1))+','+str(j)+','+str(k)
                        elif alt_track%2 == 1:
                            for ib in range(seq_len):
                                sim_grid[i,j+side_bond_size*(ib+1),k] = str(chain_list[ib])
                                if sim_grid[i,j+side_bond_size*(ib+1),k] == "P":
                                    P_count += 1
                                    P_dict[P_count] = str(i)+','+str(j+side_bond_size*(ib+1))+','+str(k)
                                elif sim_grid[i,j+side_bond_size*(ib+1),k] == "N":
                                    N_count += 1
                                    N_dict[N_count] = str(i)+','+str(j+side_bond_size*(ib+1))+','+str(k)
                                elif sim_grid[i,j+side_bond_size*(ib+1),k] == "Z":
                                    Z_count += 1
                                    Z_dict[Z_count] = str(i)+','+str(j+side_bond_size*(ib+1))+','+str(k)
                        alt_track += 1
                    elif branch_alt == 'no':
                        for ib in range(seq_len):
                            sim_grid[i+bond_size*(ib+1),j,k] = str(chain_list[ib])
                            if sim_grid[i+side_bond_size*(ib+1),j,k] == "P":
                                P_count += 1
                                P_dict[P_count] = str(i+side_bond_size*(ib+1))+','+str(j)+','+str(k)
                            elif sim_grid[i+side_bond_size*(ib+1),j,k] == "N":
                                N_count += 1
                                N_dict[N_count] = str(i+side_bond_size*(ib+1))+','+str(j)+','+str(k)
                            elif sim_grid[i+side_bond_size*(ib+1),j,k] == "Z":
                                Z_count += 1
                                Z_dict[Z_count] = str(i+side_bond_size*(ib+1))+','+str(j)+','+str(k)
                else:
                    continue

#PLACE SALT RANDOMLY INTO THE EMPTY SPACES OF THE SIMULATION GRID                                                                      

salt_ani_chg_1 = float(input_param['salt_ani_chg_1'])
salt_cat_chg_1 = float(input_param['salt_cat_chg_1'])
salt_conc_1 = float(input_param['salt_conc_1'])
salt_pairs_1 = int(salt_conc_1*Lx*Ly*top_bound*(grid_disc**3))
cat_ani_ratio_1 = abs(salt_cat_chg_1/salt_ani_chg_1)
tot_salt_ani_1 = salt_pairs_1
tot_salt_cat_1 = salt_pairs_1

salt_ani_chg_2 = float(input_param['salt_ani_chg_2'])
salt_cat_chg_2 = float(input_param['salt_cat_chg_2'])
salt_conc_2 = float(input_param['salt_conc_2'])
salt_pairs_2 = int(salt_conc_2*Lx*Ly*top_bound*(grid_disc**3))
cat_ani_ratio_2 = abs(salt_cat_chg_2/salt_ani_chg_2)
tot_salt_ani_2 = salt_pairs_2
tot_salt_cat_2 = salt_pairs_2

if cat_ani_ratio_1 < 1.0:
    tot_salt_cat_1 = int(salt_pairs_1/cat_ani_ratio_1)
elif cat_ani_ratio_1 > 1.0:
    tot_salt_ani_1 = int(salt_pairs_1*cat_ani_ratio_1)

if cat_ani_ratio_2 < 1.0:
    tot_salt_cat_2 = int(salt_pairs_2/cat_ani_ratio_2)
elif cat_ani_ratio_2 > 1.0:
    tot_salt_ani_2 = int(salt_pairs_2*cat_ani_ratio_2)

print "The system has been specific to run at a salt concentration of "+str(salt_conc_1)+" LJ sigma^-3 with charges of "+str(salt_cat_chg_1)+" and "+str(salt_ani_chg_1)+". and a salt concentration of "+str(salt_conc_2)+" LJ sigma^-3 with charegs of "+str(salt_cat_chg_2)+" and "+str(salt_ani_chg_2)+"."


#CREATE PATCHY PARTICLE, SAVE INFO, AND WRITE

num_part = 0
part_list = []
tot_chg = -8

x_shift = Lx/2.
y_shift = Ly/2.
z_shift = top_bound*0.75

#Load in icosahedron coordinates into part_list
f = open('QET_data.xyz','r')
ico_lines = f.readlines()
for i,line in zip(range(len(ico_lines)),ico_lines):
    templine = line.split()
    if i == 0:
        num_part = int(templine[0])
    if len(templine) > 0 and i != 0:
        part_list.append((float(float(templine[1])+x_shift),float(float(templine[2])+y_shift),float(float(templine[3]))+z_shift))
    else: 
        continue

f.close()

#Pick a particle at random from part_list and create a section of charged sites around it
rand1 = random.randint(0,num_part)
seed_arr = np.array([part_list[rand1]])
part_arr = np.array(part_list) 
dist_check = scipy.spatial.distance.cdist(part_arr,seed_arr)
patch1 = find_neigh(dist_check,part_list,3.0,rand1)

rest_list = diff(part_list,patch1)

num_patch_sites = len(patch1)

num_opp_sites = num_patch_sites + tot_chg

#Create number of counterions needed for the patchy particle and add them to the counterion and salt list

patch_n_count = num_patch_sites
patch_p_count = num_opp_sites

print str(len(patch1))+" negative sites identified in the patch"
print "Adding "+str(num_opp_sites)+" positive sites randomly to the particle"


#Create random list of counter ions and salt ions                                                                                      
p_dict = {}
n_dict = {}
a_dict = {}
b_dict = {}
c_dict = {}
d_dict = {}
p_count = 0
n_count = 0
a_count = 0
b_count = 0
c_count = 0
d_count = 0
p_ctr_list = ['p'] * (N_count+patch_p_count)
n_ctr_list = ['n'] * (P_count+patch_n_count)
Sp_list_1 = ['a'] * tot_salt_cat_1
Sn_list_1 = ['b'] * tot_salt_ani_1
Sp_list_2 = ['c'] * tot_salt_cat_2
Sn_list_2 = ['d'] * tot_salt_ani_2
ctr_list = p_ctr_list + n_ctr_list + Sp_list_1 + Sn_list_1 + Sp_list_2 + Sn_list_2
for i in range(5):
    ctr_list = random.sample(ctr_list,len(ctr_list))

print 'Solvating counterions and salts...'
#PLACE COUNTERIONS RANDOMLY INTO THE EMPTY SPACES OF THE SIMULATION GRID                                                               
while len(ctr_list) > 0:
    xrand = random.randint(bond_size,Lx-bond_size)
    yrand = random.randint(bond_size,Ly-bond_size)
    zrand = random.randint(bond_size,int(0.40*top_bound))

    if sim_grid[xrand,yrand,zrand] == '':
        sim_grid[xrand,yrand,zrand] = ctr_list.pop()
        if sim_grid[xrand,yrand,zrand] == 'p':
            p_count += 1
            p_dict[p_count] = str(xrand)+','+str(yrand)+','+str(zrand)
        if sim_grid[xrand,yrand,zrand] == 'n':
            n_count += 1
            n_dict[n_count] = str(xrand)+','+str(yrand)+','+str(zrand)
        if sim_grid[xrand,yrand,zrand] == 'a':
            a_count += 1
            a_dict[a_count] = str(xrand)+','+str(yrand)+','+str(zrand)
        if sim_grid[xrand,yrand,zrand] == 'b':
            b_count += 1
            b_dict[b_count] = str(xrand)+','+str(yrand)+','+str(zrand)
        if sim_grid[xrand,yrand,zrand] == 'c':
            c_count += 1
            c_dict[c_count] = str(xrand)+','+str(yrand)+','+str(zrand)
        if sim_grid[xrand,yrand,zrand] == 'd':
            d_count += 1
            d_dict[d_count] = str(xrand)+','+str(yrand)+','+str(zrand)
    else:
        continue

e_dict = {}
f_dict = {}
g_dict = {}

num_brush_atom = S_count + Z_count + P_count + N_count + n_count + p_count + a_count + b_count + c_count + d_count


opp_site_list = []

i=0
while i < num_opp_sites:
    rest_list = random.sample(rest_list,len(rest_list))
    opp_site_list.append(rest_list.pop())
    i += 1

#wfile = open('patchy_part.xyz','w')
#wfile.write(str(num_part)+"\n\n")
#for i in range(len(rest_list)):
#    wfile.write("C \t "+str(rest_list[i][0])+"\t"+str(rest_list[i][1])+"\t"+str(rest_list[i][2])+"\n")
#for i in range(len(patch1)):
#    wfile.write("N \t "+str(patch1[i][0])+"\t"+str(patch1[i][1])+"\t"+str(patch1[i][2])+"\n")
#for i in range(len(opp_site_list)):
#    wfile.write("B \t "+str(opp_site_list[i][0])+"\t"+str(opp_site_list[i][1])+"\t"+str(opp_site_list[i][2])+"\n")

#wfile.close()

#NEED TO ADJUST THIS EVERY TIME
bond_thresh = 1.4

#Use a preliminary count of the number of bonds in the brush by multipling # of chains * chain molecular weight
bond_num_est = chain_len*num_chain
atom_num_est = num_brush_atom

fin_part_list = rest_list+patch1+opp_site_list
dist_check = scipy.spatial.distance.cdist(fin_part_list,fin_part_list)
patch_bond_dict = {}
patch_len_dict = {}
patch_bond_count = 0

for i in range(0,len(dist_check),1):
    for j in range(i+1,len(dist_check),1):
        if dist_check[i,j] < bond_thresh:
            if str(i+atom_num_est+1)+','+str(j+atom_num_est+1) not in patch_bond_dict and str(j+atom_num_est+1)+','+str(i+atom_num_est+1) not in patch_bond_dict:
                patch_bond_count += 1
                patch_bond_dict[bond_num_est+patch_bond_count] = str(i+atom_num_est+1)+','+str(j+atom_num_est+1)
                patch_len_dict[bond_num_est+patch_bond_count] = round(dist_check[i,j],4)
            else:
                continue
        else:
            continue
            
#Add patchy particle atoms to total atom and bond counts
num_bond = len(patch_bond_dict)
print "There are "+str(num_bond)+" associated with the patchy particle"

e_count = len(rest_list)
f_count = len(patch1)
g_count = len(opp_site_list)


for i in range(1,e_count+1,1):
    e_dict[i] = str(fin_part_list[i-1][0])+','+str(fin_part_list[i-1][1])+','+str(fin_part_list[i-1][2])
for i in range(1,f_count+1,1):
    f_dict[i] = str(fin_part_list[i+e_count-1][0])+','+str(fin_part_list[i+e_count-1][1])+','+str(fin_part_list[i+e_count][2])
for i in range(1,g_count+1,1):
    g_dict[i] = str(fin_part_list[i+e_count+f_count-1][0])+','+str(fin_part_list[i+e_count+f_count-1][1])+','+str(fin_part_list[i+e_count+f_count-1][2])

#print num_brush_atom
#print e_dict
#print f_dict
#print g_dict

#print e_dict

#Create the full atom dictionary                                                                                                       
tot_atom_dict = {}
for i in range(1,S_count+1,1):
    tot_atom_dict[i] = S_dict[i]
for i in range(1,Z_count+1,1):
    tot_atom_dict[i+S_count] = Z_dict[i]
for i in range(1,P_count+1,1):
    tot_atom_dict[i+S_count+Z_count] = P_dict[i]
for i in range(1,N_count+1,1):
    tot_atom_dict[i+S_count+Z_count+P_count] = N_dict[i]
for i in range(1,p_count+1,1):
    tot_atom_dict[i+S_count+Z_count+P_count+N_count] = p_dict[i]
for i in range(1,n_count+1,1):
    tot_atom_dict[i+S_count+Z_count+P_count+N_count+p_count] = n_dict[i]
for i in range(1,a_count+1,1):
    tot_atom_dict[i+S_count+Z_count+P_count+N_count+p_count+n_count] = a_dict[i]
for i in range(1,b_count+1,1):
    tot_atom_dict[i+S_count+Z_count+P_count+N_count+p_count+n_count+a_count] = b_dict[i]
for i in range(1,c_count+1,1):
    tot_atom_dict[i+S_count+Z_count+P_count+N_count+p_count+n_count+a_count+b_count] = c_dict[i]
for i in range(1,d_count+1,1):
    tot_atom_dict[i+S_count+Z_count+P_count+N_count+p_count+n_count+a_count+b_count+c_count] = d_dict[i]
for i in range(1,e_count+1,1):
    tot_atom_dict[i+S_count+Z_count+P_count+N_count+p_count+n_count+a_count+b_count+c_count+d_count] = e_dict[i]
for i in range(1,f_count+1,1):
    tot_atom_dict[i+S_count+Z_count+P_count+N_count+p_count+n_count+a_count+b_count+c_count+d_count+e_count] = f_dict[i]
for i in range(1,g_count+1,1):
    tot_atom_dict[i+S_count+Z_count+P_count+N_count+p_count+n_count+a_count+b_count+c_count+d_count+e_count+f_count] = g_dict[i]
#Invert this dictionary                                                                                                                
inv_tot_atom_dict = {v: k for k,v in tot_atom_dict.items()}

bond_dict = {}

#To correctly determine all of the bonds, need to add 'Z' to the chain_list                                                            

if branch_choice == 'yes':
    chain_list.append('Z')
    chain_list.append('S')

if branch_choice == 'no':
    chain_list.append('S')

#Create new dictionary that contains all neighbors to which it is bonded                                                               
bond_count = 1

bond_dict,bond_count = bond_find(chain_list,sim_grid,inv_tot_atom_dict,Lx,Ly,top_bound,bond_size,side_bond_size,branch_choice)

print 'Removing redundant bond definitions...'

#Remove redundant bond definitions                                                                                                     
rev_bond_dict = {}
new_count = 0
if branch_choice == 'yes':
    for i in range(1,len(bond_dict)+1,1):
        if i%1000 == 0:
            print float(i)/len(bond_dict)
        curr = bond_dict[i].split(',')
        val1 = curr[0]
        val2 = curr[1]
        if val1+','+val2 in rev_bond_dict.values() or val2+','+val1 in rev_bond_dict.values():
            continue
        else:
            new_count += 1
            rev_bond_dict[new_count] = val1+','+val2


print 'Building angle dict'

#create new dictionary that contains all angles in the system                                                                          
angle_count = 1
lin_angle_dict,per_angle_dict,angle_count = angle_find(chain_list,sim_grid,inv_tot_atom_dict,Lx,Ly,top_bound,bond_size,side_bond_size,\
branch_choice)

#print len(lin_angle_dict),len(per_angle_dict)

angle_dict = merge_two_dicts(lin_angle_dict,per_angle_dict)

print 'Removing redundant angle definitions...'

#Remove redundant angle definitions                                                                                                    

rev_angle_dict = {}
new_count = 0
if branch_choice == 'yes':
    for i in range(1,len(angle_dict)+1,1):
        curr = angle_dict[i].split(',')
        val1 = curr[0]
        val2 = curr[1]
        val3 = curr[2]
        if val1+','+val2+','+val3 in rev_angle_dict.values() or val3+','+val2+','+val1 in rev_angle_dict.values():
            continue
        else:
            new_count += 1
            rev_angle_dict[new_count] = val1+','+val2+','+val3

print 'The bond search algorithm found '+str(len(bond_dict))+' bonds in the brush'
print 'The angle search algorithm found '+str(len(angle_dict))+ ' angles in the brush'

num_ang += len(angle_dict)
num_bond += len(bond_dict)

print "There are "+str(num_bond)+" total bonds in the system"

print "Substrate atoms \t "+str(S_count)
print "poly atoms Z \t\t "+str(Z_count)
print "poly atoms P \t\t "+str(P_count)
print "poly atoms N \t\t "+str(N_count)
print "pos atoms cat \t "+str(p_count)
print "neg atoms cat \t "+str(n_count)
print "pos atoms salt_1 \t "+str(a_count)
print "neg atoms salt_1 \t "+str(b_count)
print "pos atoms salt_2 \t "+str(c_count)
print "neg atoms salt_2 \t "+str(d_count)
print "Total atoms \t "+str(S_count+Z_count+P_count+N_count+p_count+n_count+a_count+b_count+c_count+d_count+e_count+f_count+g_count)
num_atom = S_count+Z_count+P_count+N_count+p_count+n_count+a_count+b_count+c_count+d_count+e_count+f_count+g_count

#Actually write the .data file                                                                                                         
print "There are "+str(num_atom)+" atoms in this system."
print "There are "+str(num_bond)+" bonds in this system."
print "There are 2 bond types in this system."

#Determine number of atom types in the system                                                                                          
atom_typ_list = []
if S_count > 0:
    atom_typ_list.append('S')
if Z_count > 0:
    atom_typ_list.append('Z')
if P_count > 0:
    atom_typ_list.append('P')
if N_count > 0:
    atom_typ_list.append('N')
if p_count > 0:
    atom_typ_list.append('p')
if n_count > 0:
    atom_typ_list.append('n')
if a_count > 0:
    atom_typ_list.append('a')
if b_count > 0:
    atom_typ_list.append('b')
if c_count > 0:
    atom_typ_list.append('c')
if d_count > 0:
    atom_typ_list.append('d')
if e_count > 0:
    atom_typ_list.append('e')
if f_count > 0:
    atom_typ_list.append('f')
if g_count > 0:
    atom_typ_list.append('g')

num_atom_type = len(atom_typ_list)
num_bond_type = 2 #one type of bond in brush, one types of bond in icosahedron

if branch_choice == 'yes':
    num_ang_type = 2
if branch_choice == 'no':
    num_ang_type = 1
num_dih_type = 0
num_imp_type = 0

#create charge, LJ, and mass dictionaries                                                                                              
chg_dict = {}
LJ_dict = {}
m_dict = {}

for item in atom_typ_list:
    if item == 'S':
        chg_dict[item] = input_param['S_chg']
        LJ_dict[item] = input_param['S_LJ']
        m_dict[item] = input_param['S_mass']
    elif item == 'Z':
        chg_dict[item] = input_param['Z_chg']
        LJ_dict[item] = input_param['Z_LJ']
        m_dict[item] = input_param['Z_mass']
    elif item == 'P':
        chg_dict[item] = input_param['P_chg']
        LJ_dict[item] = input_param['P_LJ']
        m_dict[item] = input_param['P_mass']
    elif item == 'N':
        chg_dict[item] = input_param['N_chg']
        LJ_dict[item] = input_param['N_LJ']
        m_dict[item] = input_param['N_mass']
    elif item == 'p':
        chg_dict[item] = input_param['ctr_cat_chg']
        LJ_dict[item] = input_param['ctr_cat_LJ']
        m_dict[item] = input_param['ctr_cat_mass']
    elif item == 'n':
        chg_dict[item] = input_param['ctr_ani_chg']
        LJ_dict[item] = input_param['ctr_ani_LJ']
        m_dict[item] = input_param['ctr_ani_mass']
    elif item == 'a':
        chg_dict[item] = input_param['salt_cat_chg_1']
        LJ_dict[item] = input_param['salt_cat_LJ_1']
        m_dict[item] = input_param['salt_cat_mass_1']
    elif item == 'b':
        chg_dict[item] = input_param['salt_ani_chg_1']
        LJ_dict[item] = input_param['salt_ani_LJ_1']
        m_dict[item] = input_param['salt_ani_mass_1']
    elif item == 'c':
        chg_dict[item] = input_param['salt_cat_chg_2']
        LJ_dict[item] = input_param['salt_cat_LJ_2']
        m_dict[item] = input_param['salt_cat_mass_2']
    elif item == 'd':
        chg_dict[item] = input_param['salt_ani_chg_2']
        LJ_dict[item] = input_param['salt_ani_LJ_2']
        m_dict[item] = input_param['salt_ani_mass_2']
    elif item == 'e':
        chg_dict[item] = input_param['patch_Z_chg']
        LJ_dict[item] = input_param['patch_Z_LJ']
        m_dict[item] = input_param['patch_Z_mass']
    elif item == 'f':
        chg_dict[item] = input_param['patch_N_chg']
        LJ_dict[item] = input_param['patch_N_LJ']
        m_dict[item] = input_param['patch_N_mass']
    elif item == 'g':
        chg_dict[item] = input_param['patch_P_chg']
        LJ_dict[item] = input_param['patch_P_LJ']
        m_dict[item] = input_param['patch_P_mass']

#create angle dictionary                                                                                                               
angle_coeff_dict = {}
#first dict entry rest theta, second force constant                                                                                    
angle_coeff_dict['lin'] = input_param['poly_ang_k_lin']+','+input_param['poly_ang_theta_lin']
angle_coeff_dict['per'] = input_param['poly_ang_k_per']+','+input_param['poly_ang_theta_per']

#Determine all bonds in the system.                                                                                                    
##The way this is going to work is it searches the grid for things of a particular identity (Z,P,N) that are exactly one bond length away on the grid.                                                                                                                       

print "There are "+str(num_atom_type)+" atom types in this system."
print "There are 2 bond types in this system."
print "There are 0 angle types in this system."
print "There are 0 dihedral types in this system."
print "There are 0 improper types in this system."

#print len(bond_dict)
#print len(patch_bond_dict)

b_typ_1_len = len(bond_dict)
b_typ_2_len = len(patch_bond_dict)

#print b_typ_1_len
#print b_typ_2_len

#combine brush bonds and patchy particle bonds
bond_dict = merge_two_dicts(bond_dict,patch_bond_dict)


write_data('brush.data',num_atom,num_bond,num_ang,num_dih,num_imp,sim_grid,num_atom_type,num_bond_type,num_ang_type,num_dih_type,num_imp_type,Lx,Ly,top_bound,S_dict,Z_dict,P_dict,N_dict,p_dict,n_dict,a_dict,b_dict,c_dict,d_dict,S_count,Z_count,P_count,N_count,p_count,n_count,a_count,b_count,c_count,d_count,atom_typ_list,chg_dict,LJ_dict,m_dict,grid_disc,bond_dict,lin_angle_dict,per_angle_dict,angle_coeff_dict,num_part,fin_part_list,e_count,f_count,g_count,e_dict,f_dict,g_dict,b_typ_1_len,b_typ_2_len)

write_init('brush.init')

write_settings('brush.settings',LJ_dict,atom_typ_list,lin_angle_dict,per_angle_dict,angle_coeff_dict,input_param['poly_bond_len'],input_param['poly_bond_k'],input_param['FENE_max_len'],input_param['patch_bond_k'],input_param['patch_FENE_max'])

write_infile('brush.in',input_param['tstep'],input_param['equil_steps'],input_param['sample_steps'],input_param['temp'],input_param['substr_len'],atom_typ_list,input_param['dump_int'],input_param['dielectric'],input_param['thermo_step'],Lx,Ly,top_bound,grid_disc)


