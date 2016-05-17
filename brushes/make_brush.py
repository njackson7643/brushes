import os
from functions import *
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from numpy import *
from math import sqrt
import random
#This script generates the lammps input files necessary for running polymer brush trajectories.  Setup is modeled after the simulation of Dobrynin for polyelectrolyte brushes....

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
grid_disc = 0.5 #each grid spacing corresponds to 0.1 LJ sigma
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
num_chain = int(input_param['num_chain'])
abs_graft_dens_2D = float(num_chain/(substr_len*lat_spacing)**2)
abs_graft_dens_1D = math.sqrt(abs_graft_dens_2D)
grid_graft_dens_2D = float(num_chain/(substr_len*lat_spacing/grid_disc)**2)
grid_graft_dens_1D = math.sqrt(grid_graft_dens_2D)
chain_sep_1D = int(1.0/grid_graft_dens_1D)

print "There are "+str(num_chain)+" chains spread over "+str(substr_len**2)+" substrate sites."
print "This corresponds to 1 chain every "+str(chain_sep_1D)+" grid sites, and an absolute"
print "grafting density of "+str(round(abs_graft_dens_2D,3))+" sites/(LJ dist)**2" 

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
if branch_choice == "no":
#offset = int(substr_len%(int(math.sqrt(num_chain))))
    for i in range(offset,Lx,chain_sep_1D):
        seq_count = 0
        for j in range(offset,Ly,chain_sep_1D):
            seq_count = 0
            for k in range(bond_size,poly_len+1,bond_size):
                num_bond += 1
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
                num_bond += 1 #add a bond to the total for each along the linear section
                sim_grid[i,j,k] = "Z"
                Z_count += 1
                Z_dict[Z_count] = str(i)+','+str(j)+','+str(k)
                seq_count += 1
                #The alternating branching choice is totally fucked.  To make it sort of work, I could instead just project the alternating in +x then +y then +x.  Currently, +x then -x has some serious issues.
                if seq_count%int(branch_rep) == 0:
                    if branch_alt == 'yes':
                        if alt_track%2 == 0:
                            for ib in range(seq_len):
                                num_bond += 1
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
                                num_bond += 1
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
                            num_bond += 1
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

#print Z_dict
#print "\n\n"
#print P_dict
#print "\n\n"
#print N_dict

#PLACE SALT RANDOMLY INTO THE EMPTY SPACES OF THE SIMULATION GRID

salt_ani_chg = float(input_param['salt_ani_chg'])
salt_cat_chg = float(input_param['salt_cat_chg'])
salt_conc = float(input_param['salt_conc'])
salt_pairs = int(salt_conc*Lx*Ly*top_bound*(grid_disc**3))
cat_ani_ratio = abs(salt_cat_chg/salt_ani_chg)
tot_salt_ani = salt_pairs
tot_salt_cat = salt_pairs


if cat_ani_ratio < 1.0:
    tot_salt_cat = int(salt_pairs/cat_ani_ratio)
elif cat_ani_ratio > 1.0:
    tot_salt_ani = int(salt_pairs*cat_ani_ratio)

print "The system has been specific to run at a salt concentration of "+str(salt_conc)+"LJ sigma^-3."

#Create random list of counter ions and salt ions
p_dict = {}
n_dict = {}
a_dict = {}
b_dict = {}
p_count = 0
n_count = 0
a_count = 0
b_count = 0
p_ctr_list = ['p'] * N_count
n_ctr_list = ['n'] * P_count
Sp_list = ['a'] * tot_salt_cat
Sn_list = ['b'] * tot_salt_ani
ctr_list = p_ctr_list + n_ctr_list + Sp_list + Sn_list
for i in range(5):
    ctr_list = random.sample(ctr_list,len(ctr_list))


#PLACE COUNTERIONS RANDOMLY INTO THE EMPTY SPACES OF THE SIMULATION GRID
while len(ctr_list) > 0:
    xrand = random.randint(bond_size,Lx-bond_size)
    yrand = random.randint(bond_size,Ly-bond_size)
    zrand = random.randint(bond_size,top_bound-bond_size)
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
    else:
        continue

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

#Invert this dictionary
inv_tot_atom_dict = {v: k for k,v in tot_atom_dict.items()}

bond_dict = {}

#print chain_list


#print bond_size

curr_atom_type = ''
minus_i_nn_type = ''
plus_i_nn_type = ''
minus_j_nn_type = ''
plus_j_nn_type = ''
minus_k_nn_type = ''
plus_k_nn_type = ''

#To correctly determine all of the bonds, need to add 'Z' to the chain_list


if branch_choice == 'yes':
    chain_list.append('Z')
    chain_list.append('S')

if branch_choice == 'no':
    chain_list.append('S')


ZZ_count = 0
ZP_count = 0
ZN_count = 0
ZS_count = 0
PN_count = 0

#Create new dictionary that contains all neighbors to which it is bonded
bond_count = 1

bond_dict,bond_count = bond_find(chain_list,sim_grid,inv_tot_atom_dict,Lx,Ly,top_bound,bond_size,side_bond_size)

#Remove redundant bond definitions
rev_bond_dict = {}
new_count = 0
for i in range(1,len(bond_dict)+1,1):
    curr = bond_dict[i].split(',')
    val1 = curr[0]
    val2 = curr[1]
    if val1+','+val2 in rev_bond_dict.values() or val2+','+val1 in rev_bond_dict.values():
        continue
    else:
        new_count += 1
        rev_bond_dict[new_count] = val1+','+val2 

#create new dictionary that contains all angles in the system
angle_count = 1
lin_angle_dict,per_angle_dict,angle_count = angle_find(chain_list,sim_grid,inv_tot_atom_dict,Lx,Ly,top_bound,bond_size,side_bond_size)


angle_dict = merge_two_dicts(lin_angle_dict,per_angle_dict)


#Remove redundant angle definitions
rev_angle_dict = {}
new_count = 0
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

print 'The bond search algorithm found '+str(len(rev_bond_dict))+' bonds'
print 'The angle search algorithm found '+str(len(angle_dict))+ ' angles'

num_ang = len(angle_dict)

print "Substrate atoms \t "+str(S_count)
print "poly atoms Z \t\t "+str(Z_count)
print "poly atoms P \t\t "+str(P_count)
print "poly atoms N \t\t "+str(N_count)
print "pos atoms cat \t "+str(p_count)
print "neg atoms cat \t "+str(n_count)
print "pos atoms salt \t "+str(a_count)
print "neg atoms salt \t "+str(b_count)
print "Total atoms \t "+str(S_count+Z_count+P_count+N_count+p_count+n_count+a_count+b_count)

#PLOT THE ENTIRE SIMULATION BOX
#poly_size, ctr_size, and salt_size are all specification for the size of the plotted points below.

poly_size = 20
ctr_size = 5
salt_size = 30
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#Plot all monomers and particles
for i in range(Lx):
    for j in range(Ly):
        for k in range(top_bound):
            if sim_grid[i,j,k] == "S":
                ax.scatter(grid_disc*i,grid_disc*j,grid_disc*k,c='black',s=poly_size,marker='o')
            elif sim_grid[i,j,k] == "Z":
                ax.scatter(grid_disc*i,grid_disc*j,grid_disc*k,c='0.75',s=poly_size)
            elif sim_grid[i,j,k] == "P":
                ax.scatter(grid_disc*i,grid_disc*j,grid_disc*k,c='red',s=poly_size)
            elif sim_grid[i,j,k] == "N":
                ax.scatter(grid_disc*i,grid_disc*j,grid_disc*k,c='blue',s=poly_size)
            elif sim_grid[i,j,k] == "p":
                ax.scatter(grid_disc*i,grid_disc*j,grid_disc*k,c='red',alpha=0.3,s=ctr_size)
            elif sim_grid[i,j,k] == "n":
                ax.scatter(grid_disc*i,grid_disc*j,grid_disc*k,c='blue',alpha=0.3,s=ctr_size)       
            elif sim_grid[i,j,k] == "a":
                ax.scatter(grid_disc*i,grid_disc*j,grid_disc*k,c='purple',alpha=0.8,s=salt_size)
            elif sim_grid[i,j,k] == "b":
                ax.scatter(grid_disc*i,grid_disc*j,grid_disc*k,c='green',alpha=0.8,s=salt_size)        
            else:
                continue

#Plot all bonds
line_thick = 2.
for i in range(1,len(rev_bond_dict)+1,1):
    curr = rev_bond_dict[i].split(',')
    xyz_1 = tot_atom_dict[int(curr[0])].split(',')
    xyz_2 = tot_atom_dict[int(curr[1])].split(',')
    x1 = grid_disc*float(xyz_1[0])
    x2 = grid_disc*float(xyz_2[0])
    y1 = grid_disc*float(xyz_1[1])
    y2 = grid_disc*float(xyz_2[1])
    z1 = grid_disc*float(xyz_1[2])
    z2 = grid_disc*float(xyz_2[2])
    ax.plot([x1,x2],[y1,y2],[z1,z2],linewidth=line_thick,c='black')

ax.view_init(elev=0.,azim=45)
ax.set_xlim3d(0,Lx*grid_disc/2)
ax.set_ylim3d(0,Lx*grid_disc/2)
ax.set_zlim3d(0,top_bound*grid_disc)
plt.show()

num_atom = S_count+Z_count+P_count+N_count+p_count+n_count+a_count+b_count

#Actually write the .data file
print "There are "+str(num_atom)+" atoms in this system."
print "There are "+str(num_bond)+" bonds in this system."
print "There is only one bond type in this system."

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

num_atom_type = len(atom_typ_list)
num_bond_type = 1

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
        chg_dict[item] = input_param['salt_cat_chg']
        LJ_dict[item] = input_param['salt_cat_LJ']
        m_dict[item] = input_param['salt_cat_mass']
    elif item == 'b':
        chg_dict[item] = input_param['salt_ani_chg']
        LJ_dict[item] = input_param['salt_ani_LJ']
        m_dict[item] = input_param['salt_ani_mass']

#create angle dictionary
angle_coeff_dict = {}
#first dict entry rest theta, second force constant
angle_coeff_dict['lin'] = input_param['poly_ang_k_lin']+','+input_param['poly_ang_theta_lin']
angle_coeff_dict['per'] = input_param['poly_ang_k_per']+','+input_param['poly_ang_theta_per']

#Determine all bonds in the system.
##The way this is going to work is it searches the grid for things of a particular identity (Z,P,N) that are exactly one bond length away on the grid.
            

print "There are "+str(num_atom_type)+" atom types in this system."
print "There are 1 bond types in this system."
print "There are 0 angle types in this system."
print "There are 0 dihedral types in this system."
print "There are 0 improper types in this system."

write_data('brush.data',num_atom,num_bond,num_ang,num_dih,num_imp,sim_grid,num_atom_type,num_bond_type,num_ang_type,num_dih_type,num_imp_type,Lx,Ly,top_bound,S_dict,Z_dict,P_dict,N_dict,p_dict,n_dict,a_dict,b_dict,S_count,Z_count,P_count,N_count,p_count,n_count,a_count,b_count,atom_typ_list,chg_dict,LJ_dict,m_dict,grid_disc,rev_bond_dict,lin_angle_dict,per_angle_dict,angle_coeff_dict)

write_init('brush.init')

write_settings('brush.settings',LJ_dict,atom_typ_list,lin_angle_dict,per_angle_dict,angle_coeff_dict,input_param['poly_bond_len'],input_param['poly_bond_k'],input_param['FENE_max_len'])

write_infile('brush.in',input_param['tstep'],input_param['equil_steps'],input_param['sample_steps'],input_param['temp'],input_param['substr_len'],atom_typ_list,input_param['dump_int'],input_param['dielectric'])
#Write .init file

#Write .settings file

#Write .in file

#Define bottom substrate group and freeze particles.
#Define top substrate group and freeze particles.
#Define polymer group
#Define counterion group



