import os
import sys
import math
import random

def setup_parser(filename):
    if os.path.isfile(filename) == False:
        print "Could not locate "+str(filename)+", exiting program..."
        sys.exit()
    keys = ['num_chain','chain_len','substr_len','substr_patt','chain_typ','equil_steps','sample_steps','tstep','temp','lat_spacing','poly_LJ','ctr_cat_chg','ctr_ani_chg','ctr_cat_LJ','ctr_ani_LJ','salt_cat_chg_1','salt_ani_chg_1','salt_cat_LJ_1','salt_ani_LJ_1','poly_bond_len','branch','branch_rep','branch_alt','salt_conc_1','P_mass','N_mass','Z_mass','P_LJ','N_LJ','Z_LJ','ctr_ani_mass','ctr_cat_mass','salt_ani_mass_1','salt_cat_mass_1','S_mass','S_chg','S_LJ','Z_chg','P_chg','N_chg','poly_bond_k','poly_ang_theta_lin','poly_ang_theta_per','poly_ang_k_lin','poly_ang_k_per','dump_int','FENE_max_len','dielectric','salt_conc_2','salt_cat_chg_2','salt_ani_chg_2','salt_cat_LJ_2','salt_ani_LJ_2','salt_cat_mass_2','salt_ani_mass_2','thermo_step','init_substr_sep','final_substr_sep']
    indices = list(range(len(keys)))
    hash = {k:i for k,i in zip(keys,indices)}
    input_param = range(len(keys))
    setup_f = open(filename,'r')
    for line in setup_f:
        templine = line.split()
        if line[0] == '#':
            continue
        elif len(templine) == 3:
            curr_key = templine[0]
            curr_val = templine[2]
            if curr_key in keys:
                input_param[hash[curr_key]] = curr_val
            else:
                print "Error: "+str(curr_key)+" is not a valid keyword"
        else:
            continue
    #return a hashed table with the input parameters and the read-in values
    return {k:v for k,v in zip(keys,input_param)}

def write_data(filename,num_atom,num_bond,num_ang,num_dih,num_imp,sim_grid,num_atom_type,num_bond_type,num_ang_type,num_dih_type,num_imp_type,Lx,Ly,top_bound,S_dict,Z_dict,P_dict,N_dict,p_dict,n_dict,a_dict,b_dict,c_dict,d_dict,S_count,Z_count,P_count,N_count,p_count,n_count,a_count,b_count,c_count,d_count,atom_type_list,chg_dict,LJ_dict,m_dict,grid_disc,rev_bond_dict,lin_angle_dict,per_angle_dict,angle_coeff_dict):
    wfile = open(filename,'w')
    wfile.write("LAMMPS Description\n\n")
    wfile.write("\t "+str(num_atom)+" atoms\n")
    wfile.write("\t "+str(num_bond)+" bonds\n")
    wfile.write("\t "+str(num_ang)+" angles\n")
    wfile.write("\t "+str(num_dih)+" dihedrals\n")
    wfile.write("\t "+str(num_imp)+" impropers\n\n")
    wfile.write("\t "+str(num_atom_type)+" atom types\n")
    wfile.write("\t "+str(num_bond_type)+" bond types\n")
    wfile.write("\t "+str(num_ang_type)+" angle types\n")
    wfile.write("\t "+str(num_dih_type)+" dihedral types\n")
    wfile.write("\t "+str(num_imp_type)+" improper types\n\n")
    wfile.write(str(0.0)+"  "+str(Lx*grid_disc)+" xlo xhi"+"\n")
    wfile.write(str(0.0)+"  "+str(Ly*grid_disc)+" ylo yhi"+"\n")
    wfile.write(str(0.0)+"  "+str(top_bound*grid_disc)+" zlo zhi \n"+"\n")
    wfile.write('Masses\n\n')
    for i in range(1,num_atom_type+1,1):
        wfile.write(str(i)+"\t "+str(m_dict[atom_type_list[i-1]]))
        wfile.write("\n")
    wfile.write("\nAtoms\n")
    wfile.write("\n")
    typ_count = 0
    mol_count = 1
#Write bottom substrate atoms
    if S_count > 0:
        for i in range(1,(S_count/2)+1,1):
            curr = S_dict[i].split(',')
            wfile.write(str(i)+" "+str(mol_count)+" "+str(typ_count+1)+"\t"+str(chg_dict[atom_type_list[typ_count]])+"\t"+str(float(curr[0])*grid_disc)+"\t"+str(float(curr[1])*grid_disc)+"\t"+str(float(curr[2])*grid_disc)+"\n")
#Write top substrate atoms
        for i in range((S_count/2)+1,S_count+1,1):
            curr = S_dict[i].split(',')
            wfile.write(str(i)+" "+str(mol_count)+" "+str(typ_count+1)+"\t"+str(chg_dict[atom_type_list[typ_count]])+"\t"+str(float(curr[0])*grid_disc)+"\t"+str(float(curr[1])*grid_disc)+"\t"+str(float(curr[2])*grid_disc)+"\n")
#Write Z atoms
    if Z_count > 0:    
        mol_count += 1
        typ_count += 1
        for i in range(1,Z_count+1,1):
            curr = Z_dict[i].split(',')
            wfile.write(str(i+S_count)+" "+str(mol_count)+" "+str(typ_count+1)+"\t"+str(chg_dict[atom_type_list[typ_count]])+"\t"+str(float(curr[0])*grid_disc)+"\t"+str(float(curr[1])*grid_disc)+"\t"+str(float(curr[2])*grid_disc)+"\n")
#Write P atoms
    if P_count > 0:
        mol_count += 1
        typ_count += 1
        for i in range(1,P_count+1,1):
            curr = P_dict[i].split(',')
            wfile.write(str(i+S_count+Z_count)+" "+str(mol_count)+" "+str(typ_count+1)+"\t"+str(chg_dict[atom_type_list[typ_count]])+"\t"+str(float(curr[0])*grid_disc)+"\t"+str(float(curr[1])*grid_disc)+"\t"+str(float(curr[2])*grid_disc)+"\n")
#Write N atoms
    if N_count > 0:
        mol_count += 1
        typ_count += 1
        for i in range(1,N_count+1,1):
            curr = N_dict[i].split(',')
            wfile.write(str(i+S_count+Z_count+P_count)+" "+str(mol_count)+" "+str(typ_count+1)+"\t"+str(chg_dict[atom_type_list[typ_count]])+"\t"+str(float(curr[0])*grid_disc)+"\t"+str(float(curr[1])*grid_disc)+"\t"+str(float(curr[2])*grid_disc)+"\n")
#Write p atoms
    if p_count > 0:
        mol_count += 1
        typ_count += 1
        for i in range(1,p_count+1,1):
            curr = p_dict[i].split(',')
            wfile.write(str(i+S_count+Z_count+P_count+N_count)+" "+str(mol_count)+" "+str(typ_count+1)+"\t"+str(chg_dict[atom_type_list[typ_count]])+"\t"+str(float(curr[0])*grid_disc)+"\t"+str(float(curr[1])*grid_disc)+"\t"+str(float(curr[2])*grid_disc)+"\n")
#Write n atoms
    if n_count > 0:
        mol_count += 1
        typ_count += 1
        for i in range(1,n_count+1,1):
            curr = n_dict[i].split(',')
            wfile.write(str(i+S_count+Z_count+P_count+N_count+p_count)+" "+str(mol_count)+" "+str(typ_count+1)+"\t"+str(chg_dict[atom_type_list[typ_count]])+"\t"+str(float(curr[0])*grid_disc)+"\t"+str(float(curr[1])*grid_disc)+"\t"+str(float(curr[2])*grid_disc)+"\n")
#Write a atoms
    if a_count > 0:
        mol_count += 1
        typ_count += 1
        for i in range(1,a_count+1,1):
            curr = a_dict[i].split(',')
            wfile.write(str(i+S_count+Z_count+P_count+N_count+p_count+n_count)+" "+str(mol_count)+" "+str(typ_count+1)+"\t"+str(chg_dict[atom_type_list[typ_count]])+"\t"+str(float(curr[0])*grid_disc)+"\t"+str(float(curr[1])*grid_disc)+"\t"+str(float(curr[2])*grid_disc)+"\n")
#Write b atoms
    if b_count > 0:
        mol_count += 1
        typ_count += 1
        for i in range(1,b_count+1,1):
            curr = b_dict[i].split(',')
            wfile.write(str(i+S_count+Z_count+P_count+N_count+p_count+n_count+a_count)+" "+str(mol_count)+" "+str(typ_count+1)+"\t"+str(chg_dict[atom_type_list[typ_count]])+"\t"+str(float(curr[0])*grid_disc)+"\t"+str(float(curr[1])*grid_disc)+"\t"+str(float(curr[2])*grid_disc)+"\n")
#Write c atoms
    if c_count > 0:
        mol_count += 1
        typ_count += 1
        for i in range(1,c_count+1,1):
            curr = c_dict[i].split(',')
            wfile.write(str(i+S_count+Z_count+P_count+N_count+p_count+n_count+a_count+b_count)+" "+str(mol_count)+" "+str(typ_count+1)+"\t"+str(chg_dict[atom_type_list[typ_count]])+"\t"+str(float(curr[0])*grid_disc)+"\t"+str(float(curr[1])*grid_disc)+"\t"+str(float(curr[2])*grid_disc)+"\n")
#Write d atoms
    if d_count > 0:
        mol_count += 1
        typ_count += 1
        for i in range(1,d_count+1,1):
            curr = d_dict[i].split(',')
            wfile.write(str(i+S_count+Z_count+P_count+N_count+p_count+n_count+a_count+b_count+c_count)+" "+str(mol_count)+" "+str(typ_count+1)+"\t"+str(chg_dict[atom_type_list[typ_count]])+"\t"+str(float(curr[0])*grid_disc)+"\t"+str(float(curr[1])*grid_disc)+"\t"+str(float(curr[2])*grid_disc)+"\n")
#Write BONDS section
    wfile.write("\n")
    wfile.write("Bonds\n\n")
#    print rev_bond_dict
    for i in range(1,len(rev_bond_dict)+1,1):
        curr = rev_bond_dict[i].split(',')
        wfile.write(str(i)+"\t"+str(1)+"\t"+str(curr[0])+"\t"+str(curr[1])+"\n")
    wfile.write("\n")
    wfile.write("Angles\n\n")
    for i in range(1,len(lin_angle_dict)+len(per_angle_dict)+1,1):
        if i in lin_angle_dict:
            curr = lin_angle_dict[i].split(',')
            wfile.write(str(i)+" "+str(1)+" "+str(curr[0])+" "+str(curr[1])+" "+str(curr[2])+"\n")
        elif i in per_angle_dict:
            curr = per_angle_dict[i].split(',')
            wfile.write(str(i)+" "+str(2)+" "+str(curr[0])+" "+str(curr[1])+" "+str(curr[2])+"\n")
        
def write_init(filename):
    wfile = open(filename,'w')
    wfile.write("units \t\t\t lj \n")
    wfile.write("atom_style \t\t full \n")
#Set short range LJ cutoff to 2.5 and use shift.  Electrostatics cutoff set to > 2x LJ cutoff
    wfile.write("pair_style \t\t lj/cut/coul/long 2.5 5.5 \n")
    wfile.write("pair_modify \t\t shift yes mix arithmetic \n")
    wfile.write("bond_style \t\t fene \n")
    wfile.write("angle_style \t\t cosine/delta \n")
    wfile.write("boundary \t\t p p f \n")
    wfile.write("neighbor \t\t 0.5 bin \n")
    wfile.write("neigh_modify \t\t every 1 delay 3 check yes one 10000 \n")
    wfile.write("kspace_style \t\t pppm   0.005 \n")
    wfile.write("kspace_modify \t\t slab   3.0 \n")
    wfile.write("special_bonds \t\t fene \n\n")

def write_settings(filename,LJ_dict,atom_type_list,lin_angle_dict,per_angle_dict,angle_coeff_dict,poly_bond_len,poly_bond_k,FENE_bond_len):
    wfile = open(filename,'w')
    wfile.write("#Non-bonded interactions (pair-wise) \n")
    for i in range(len(atom_type_list)):
        for j in range(i,len(atom_type_list)):
            LJ_eps_i = float(LJ_dict[atom_type_list[i]].split(',')[0])
            LJ_sig_i = float(LJ_dict[atom_type_list[i]].split(',')[1])
            LJ_eps_j = float(LJ_dict[atom_type_list[j]].split(',')[0])
            LJ_sig_j = float(LJ_dict[atom_type_list[j]].split(',')[1])
            LJ_eps_ij = round(math.sqrt(LJ_eps_i*LJ_eps_j),6)
            LJ_sig_ij = round((LJ_sig_i+LJ_sig_j)/2.,6)
            wfile.write("pair_coeff "+str(i+1)+" "+str(j+1)+" "+str(LJ_eps_ij)+" "+str(LJ_sig_ij)+"\n")
    wfile.write("\n#FENE Stretching Interactions \n")
    wfile.write("#bond_coeff bondtype k R0 epsilon sigma\n")
    wfile.write("bond_coeff 1 "+str(poly_bond_k)+" "+str(FENE_bond_len)+" 1.0 1.0\n\n")
    wfile.write("\n#Harmonic Angle Interaction \n")
    wfile.write("angle_coeff 1 "+angle_coeff_dict['lin'].split(',')[0]+" "+angle_coeff_dict['lin'].split(',')[1]+"\n")
    if len(lin_angle_dict)+len(per_angle_dict) != len(lin_angle_dict):
        wfile.write("angle_coeff 2 "+angle_coeff_dict['per'].split(',')[0]+" "+angle_coeff_dict['per'].split(',')[1]+"\n")

def write_infile(filename,tstep,equil_steps,sample_steps,temp,substr_len,atom_type_list,dump_int,dielectric,thermo_step):
    vel_seed1 = random.randint(1,99999999)
    vel_seed2 = random.randint(1,99999999)
    wfile = open(filename,'w')
    wfile.write("#------------- Init Section ---------------\n\n")
    wfile.write("include "+'"'+filename[:-3]+'.init"\n')
    wfile.write("#------------- Particle Definitions ------------ \n\n")
    wfile.write("read_data "+'"'+filename[:-3]+'.data"\n')
    wfile.write("#------------- Settings Section ----------- \n\n")
    wfile.write("include "+'"'+filename[:-3]+'.settings"\n')
    wfile.write("#------------- Run Section ----------- \n\n")
    wfile.write("thermo "+str(thermo_step)+"\n")
    wfile.write("run_style verlet\n")
    wfile.write("timestep "+str(tstep)+"\n\n")
    wfile.write("restart 50000 brush.restart\n\n")
    wfile.write("#SIMULATION BOX FIXES\n\n")
    wfile.write("group substrates type 1 \n")
    wfile.write("group bot_substr id <= "+str(int(substr_len)**2)+'\n')
    wfile.write("group top_substr subtract substrates bot_substr\n")
    if "Z" in atom_type_list and "P" in atom_type_list and "N" in atom_type_list:
        wfile.write("group polymers type 2:4\n")
        wfile.write("group ctr_ions type 5:6\n")
    elif "Z" in atom_type_list and "P" in atom_type_list:
        wfile.write("group polymers type 2:3\n")
        wfile.write("group ctr_ions type 4\n")
    elif "Z" in atom_type_list and "N" in atom_type_list:
        wfile.write("group polymers type 2:3\n")
        wfile.write("group ctr_ions type 4\n")
    elif "P" in atom_type_list and "N" in atom_type_list:
        wfile.write("group polymers type 2:3\n")
        wfile.write("group ctr_ions type 4:5\n")
    elif "Z" in atom_type_list:
        wfile.write("group polymers type 2\n")
    elif "P" in atom_type_list:
        wfile.write("group polymers type 2\n")
        wfile.write("group ctr_ions type 3\n")
    elif "N" in atom_type_list:
        wfile.write("group polymers type 2\n")
        wfile.write("group ctr_ions type 3\n")
    if "a" in atom_type_list or "b" in atom_type_list:
        wfile.write("group salt subtract all substrates polymers ctr_ions\n")
    wfile.write("group dump_group subtract all top_substr\n")
    wfile.write("fix 1 substrates setforce 0.0 0.0 0.0\n")
    wfile.write("group not_substr subtract all substrates\n")
    wfile.write("fix wall1 not_substr wall/lj126 zlo EDGE 0.1 1.0 2.5 \n")
    wfile.write("fix wall2 not_substr wall/lj126 zhi EDGE 0.1 1.0 2.5 \n\n")
    wfile.write("compute real_temp not_substr temp\n")
    wfile.write("thermo_style custom step dt c_real_temp press vol etotal ke pe ebond eangle evdwl ecoul elong\n\n")
    wfile.write("#Minimize the simulation box. \n")
    wfile.write("fix poly_hold polymers setforce 0.0 0.0 0.0\n")
    wfile.write("minimize 1.0e-6 1.0e-6 2000 2000\n\n")
    wfile.write("unfix poly_hold\n\n")
    wfile.write("#Initial Safe Equilibration to remove bad contacts\n")
    wfile.write("velocity not_substr create "+str(temp)+" "+str(vel_seed1)+"\n")
    wfile.write("fix temper not_substr nve/limit 0.1\n")
    wfile.write("fix temper2 not_substr langevin "+str(temp)+" "+str(temp)+" 100.0 986537\n")
    wfile.write("fix rescale0 not_substr temp/rescale 2 1.0 1.0 0.2 1.0\n")
    wfile.write("dump 1 dump_group custom "+str(dump_int)+" equil.trj id type x y z\n")
    wfile.write("run "+str(equil_steps)+"\n")
    wfile.write("unfix rescale0\n")
    wfile.write("unfix temper2\n")
    wfile.write("unfix temper\n")
    wfile.write("undump 1\n\n")
    wfile.write("#Run NVT Sampling\n")
    wfile.write("fix 11 not_substr nve\n")
    wfile.write("fix 3 not_substr langevin "+str(temp)+" "+str(temp)+ " 100.0 "+str(vel_seed2)+"\n")
    wfile.write("dump 2 polymers custom "+str(dump_int)+" polymers.trj id type q xu yu zu\n")
    wfile.write("dump 55 polymers custom "+str(dump_int)+" poly_wrap.trj id type q x y z\n")
    if "P" in atom_type_list or "N" in atom_type_list:
        wfile.write("dump 3 ctr_ions custom "+str(dump_int)+" ctrions.trj id type q x y z\n")
    if "a" in atom_type_list or "b" in atom_type_list or "c" in atom_type_list or "d" in atom_type_list:
        wfile.write("dump 4 salt custom "+str(dump_int)+" salt.trj id type q x y z\n")
    wfile.write("run "+str(sample_steps)+"\n")
    wfile.write("unfix 3\n")
    wfile.write("unfix 11\n")
    wfile.write("undump 2\n\n")
    wfile.write("undump 55\n\n")
    if "P" in atom_type_list or "N" in atom_type_list:
        wfile.write("undump 3\n")
    if "a" in atom_type_list or "b" in atom_type_list:
        wfile.write("undump 4\n")
    

#    wfile.write("velocity not_substr all create 1.0 094376\n")
    
 
def bond_find(chain_list,sim_grid,inv_tot_atom_dict,Lx,Ly,top_bound,bond_size,side_bond_size,branch_choice):
    if branch_choice == 'yes':
        bond_count = 1
        bond_dict = {}
        for i in range(Lx):
            for j in range(Ly):
                for k in range(top_bound/2):
                    if sim_grid[i,j,k] != '' and sim_grid[i,j,k] != 'p' and sim_grid[i,j,k] != 'n' and sim_grid[i,j,k] != 'a' and sim_grid[i,j,k] != 'b' and sim_grid[i,j,k] != 'c' and sim_grid[i,j,k] != 'd':
                        curr_num = inv_tot_atom_dict[str(i)+','+str(j)+','+str(k)]
                        curr_type = sim_grid[i,j,k]
                        i_m_num = ''
                        i_m_type = ''
                        i_p_num = ''
                        i_p_type = ''
                        j_m_num = ''
                        j_m_type = ''
                        j_p_num = ''
                        j_p_type = ''
                        k_m_num = ''
                        k_m_type = ''
                        k_p_num = ''
                        k_p_type = ''
                        if i < side_bond_size:
                            i_m = Lx-side_bond_size-i 
                            i_p = i+side_bond_size
                        elif i > Lx-side_bond_size-1:
                            i_m = i-side_bond_size
                            i_p = i+side_bond_size-Lx
                        else:
                            i_m = i-side_bond_size
                            i_p = i+side_bond_size
                        if str(i_m)+','+str(j)+','+str(k) in inv_tot_atom_dict:
                            i_m_num = inv_tot_atom_dict[str(i_m)+','+str(j)+','+str(k)]
                            i_m_type = sim_grid[i_m,j,k]
                        if str(i_p)+','+str(j)+','+str(k) in inv_tot_atom_dict:
                            i_p_num = inv_tot_atom_dict[str(i_p)+','+str(j)+','+str(k)]
                            i_p_type = sim_grid[i_p,j,k]
                        if j < bond_size:
                            j_m = Ly-side_bond_size-j 
                            j_p = j+side_bond_size
                        elif j > Ly-side_bond_size-1:
                            j_m = j-side_bond_size
                            j_p = j+side_bond_size-Ly
                        else:
                            j_m = j-side_bond_size
                            j_p = j+side_bond_size
                        if str(i)+','+str(j_m)+','+str(k) in inv_tot_atom_dict:
                            j_m_num = inv_tot_atom_dict[str(i)+','+str(j_m)+','+str(k)]
                            j_m_type = sim_grid[i,j_m,k]
                        if str(i)+','+str(j_p)+','+str(k) in inv_tot_atom_dict:
                            j_p_num = inv_tot_atom_dict[str(i)+','+str(j_p)+','+str(k)]
                            j_p_type = sim_grid[i,j_p,k]
                        if k < bond_size:
                            k_m = top_bound-bond_size-k 
                            k_p = k+bond_size
                        elif k > top_bound-bond_size-1:
                            k_m = k-bond_size
                            k_p = k+bond_size-top_bound
                        else:
                            k_m = k-bond_size
                            k_p = k+bond_size
                        if str(i)+','+str(j)+','+str(k_m) in inv_tot_atom_dict:
                            k_m_num = inv_tot_atom_dict[str(i)+','+str(j)+','+str(k_m)]
                            k_m_type = sim_grid[i,j,k_m]
                        if str(i)+','+str(j)+','+str(k_p) in inv_tot_atom_dict:
                            k_p_num = inv_tot_atom_dict[str(i)+','+str(j)+','+str(k_p)]
                            k_p_type = sim_grid[i,j,k_p]
                    #Add bonds to dictionary
                        if curr_type in chain_list and i_m_type in chain_list:
                            bond_dict[bond_count] = str(curr_num)+','+str(i_m_num)
                            bond_count += 1
                        if curr_type in chain_list and i_p_type in chain_list:
                            bond_dict[bond_count] = str(curr_num)+','+str(i_p_num)
                            bond_count += 1
                        if curr_type in chain_list and j_m_type in chain_list:
                            bond_dict[bond_count] = str(curr_num)+','+str(j_m_num)
                            bond_count += 1
                        if curr_type in chain_list and j_p_type in chain_list:
                            bond_dict[bond_count] = str(curr_num)+','+str(j_p_num)
                            bond_count += 1
                        if curr_type in chain_list and k_m_type in chain_list:
                            bond_dict[bond_count] = str(curr_num)+','+str(k_m_num)
                            bond_count += 1
                        if curr_type in chain_list and k_p_type in chain_list:
                            bond_dict[bond_count] = str(curr_num)+','+str(k_p_num)
                            bond_count += 1                        
                    else:
                        continue
    elif branch_choice == 'no':
        bond_count = 1
        bond_dict = {}
        for i in range(Lx):
            for j in range(Ly):
                for k in range(top_bound/2):
                    if sim_grid[i,j,k] != '' and sim_grid[i,j,k] != 'p' and sim_grid[i,j,k] != 'n' and sim_grid[i,j,k] != 'a' and sim_grid[i,j,k] != 'b' and sim_grid[i,j,k] != 'c' and sim_grid[i,j,k] != 'd':
                        curr_num = inv_tot_atom_dict[str(i)+','+str(j)+','+str(k)]
                        curr_type = sim_grid[i,j,k]
                        k_p_num = ''
                        k_p_type = ''
                        if k < bond_size:
                            k_p = k+bond_size
                        elif k > top_bound-bond_size-1:
                            k_p = i+bond_size-top_bound
                        else:
                            k_p = k+bond_size
                        if str(i)+','+str(j)+','+str(k_p) in inv_tot_atom_dict:
                            k_p_num = inv_tot_atom_dict[str(i)+','+str(j)+','+str(k_p)]
                            k_p_type = sim_grid[i,j,k_p]
                    #Add bonds to dictionary
                        if curr_type in chain_list and k_p_type in chain_list:
                            bond_dict[bond_count] = str(curr_num)+','+str(k_p_num)
                            bond_count += 1                        
                    else:
                        continue

    return bond_dict,bond_count

def angle_find(chain_list,sim_grid,inv_tot_atom_dict,Lx,Ly,top_bound,bond_size,side_bond_size,branch_choice):
    if branch_choice == 'yes':
        angle_count = 1
        lin_angle_dict = {}
        per_angle_dict = {}
        for i in range(Lx):
            for j in range(Ly):
            #Starting at 3 in the z direction ensures that the angle potentials don't connect to the substrate
            #Start k at 1 if you want it to connect through the substrate
                for k in range(3,top_bound/2):
                    if sim_grid[i,j,k] != '' and sim_grid[i,j,k] != 'p' and sim_grid[i,j,k] != 'n' and sim_grid[i,j,k] != 'a' and sim_grid[i,j,k] != 'b' and sim_grid[i,j,k] != 'c' and sim_grid[i,j,k] != 'd':
                        curr_num = inv_tot_atom_dict[str(i)+','+str(j)+','+str(k)]
                        curr_type = sim_grid[i,j,k]
                        i_m_num = ''
                        i_m_type = ''
                        i_p_num = ''
                        i_p_type = ''
                        j_m_num = ''
                        j_m_type = ''
                        j_p_num = ''
                        j_p_type = ''
                        k_m_num = ''
                        k_m_type = ''
                        k_p_num = ''
                        k_p_type = ''
                        if i < side_bond_size:
                            i_m = Lx-side_bond_size-i 
                            i_p = i+side_bond_size
                        elif i > Lx-side_bond_size-1:
                            i_m = i-side_bond_size
                            i_p = i+side_bond_size-Lx
                        else:
                            i_m = i-side_bond_size
                            i_p = i+side_bond_size
                        if str(i_m)+','+str(j)+','+str(k) in inv_tot_atom_dict:
                            i_m_num = inv_tot_atom_dict[str(i_m)+','+str(j)+','+str(k)]
                            i_m_type = sim_grid[i_m,j,k]
                        if str(i_p)+','+str(j)+','+str(k) in inv_tot_atom_dict:
                            i_p_num = inv_tot_atom_dict[str(i_p)+','+str(j)+','+str(k)]
                            i_p_type = sim_grid[i_p,j,k]
                        if j < side_bond_size:
                            j_m = Ly-side_bond_size-j 
                            j_p = j+side_bond_size
                        elif j > Ly-side_bond_size-1:
                            j_m = j-side_bond_size
                            j_p = j+side_bond_size-Ly
                        else:
                            j_m = j-side_bond_size
                            j_p = j+side_bond_size
                        if str(i)+','+str(j_m)+','+str(k) in inv_tot_atom_dict:
                            j_m_num = inv_tot_atom_dict[str(i)+','+str(j_m)+','+str(k)]
                            j_m_type = sim_grid[i,j_m,k]
                        if str(i)+','+str(j_p)+','+str(k) in inv_tot_atom_dict:
                            j_p_num = inv_tot_atom_dict[str(i)+','+str(j_p)+','+str(k)]
                            j_p_type = sim_grid[i,j_p,k]
                        if k < bond_size:
                            k_m = top_bound-bond_size-k 
                            k_p = k+bond_size
                        elif k > top_bound-bond_size-1:
                            k_m = k-bond_size
                            k_p = k+bond_size-top_bound
                        else:
                            k_m = k-bond_size
                            k_p = k+bond_size
                        if str(i)+','+str(j)+','+str(k_m) in inv_tot_atom_dict:
                            k_m_num = inv_tot_atom_dict[str(i)+','+str(j)+','+str(k_m)]
                            k_m_type = sim_grid[i,j,k_m]
                        if str(i)+','+str(j)+','+str(k_p) in inv_tot_atom_dict:
                            k_p_num = inv_tot_atom_dict[str(i)+','+str(j)+','+str(k_p)]
                            k_p_type = sim_grid[i,j,k_p]
                    #Add angles to dictionary
                    #Form all angle potentials along linear chain segments
                        if curr_type in chain_list and i_m_type in chain_list and i_p_type in chain_list:
                            lin_angle_dict[angle_count] = str(i_m_num)+','+str(curr_num)+','+str(i_p_num)
                            angle_count += 1
                        if curr_type in chain_list and j_m_type in chain_list and j_p_type in chain_list:
                            lin_angle_dict[angle_count] = str(j_m_num)+','+str(curr_num)+','+str(j_p_num)
                        angle_count += 1
                        if curr_type in chain_list and k_m_type in chain_list and k_p_type in chain_list:
                            lin_angle_dict[angle_count] = str(k_m_num)+','+str(curr_num)+','+str(k_p_num)
                            angle_count += 1
                    #Form all angle potential at branch points
                        if curr_type in chain_list and i_m_type in chain_list and k_m_type in chain_list:
                            per_angle_dict[angle_count] = str(i_m_num)+','+str(curr_num)+','+str(k_m_num)
                            angle_count += 1
                        if curr_type in chain_list and i_p_type in chain_list and k_m_type in chain_list:
                            per_angle_dict[angle_count] = str(i_p_num)+','+str(curr_num)+','+str(k_m_num)
                            angle_count += 1
                        if curr_type in chain_list and j_m_type in chain_list and k_m_type in chain_list:
                            per_angle_dict[angle_count] = str(j_m_num)+','+str(curr_num)+','+str(k_m_num)
                            angle_count += 1
                        if curr_type in chain_list and j_p_type in chain_list and k_m_type in chain_list:
                            per_angle_dict[angle_count] = str(j_p_num)+','+str(curr_num)+','+str(k_m_num)
                            angle_count += 1
                        if curr_type in chain_list and i_m_type in chain_list and k_p_type in chain_list:
                            per_angle_dict[angle_count] = str(i_m_num)+','+str(curr_num)+','+str(k_p_num)
                            angle_count += 1
                        if curr_type in chain_list and i_p_type in chain_list and k_p_type in chain_list:
                            per_angle_dict[angle_count] = str(i_p_num)+','+str(curr_num)+','+str(k_p_num)
                            angle_count += 1
                        if curr_type in chain_list and j_m_type in chain_list and k_p_type in chain_list:
                            per_angle_dict[angle_count] = str(j_m_num)+','+str(curr_num)+','+str(k_p_num)
                            angle_count += 1
                        if curr_type in chain_list and j_p_type in chain_list and k_p_type in chain_list:
                            per_angle_dict[angle_count] = str(j_p_num)+','+str(curr_num)+','+str(k_p_num)
                            angle_count += 1
                    else:
                        continue
    elif branch_choice == 'no':
        angle_count = 1
        lin_angle_dict = {}
        per_angle_dict = {}
        for i in range(Lx):
            for j in range(Ly):
            #Starting at 3 in the z direction ensures that the angle potentials don't connect to the substrate
            #Start k at 1 if you want it to connect through the substrate
                for k in range(2,top_bound/2):
                    if sim_grid[i,j,k] != '' and sim_grid[i,j,k] != 'p' and sim_grid[i,j,k] != 'n' and sim_grid[i,j,k] != 'a' and sim_grid[i,j,k] != 'b' and sim_grid[i,j,k] != 'c' and sim_grid[i,j,k] != 'd':
                        curr_num = inv_tot_atom_dict[str(i)+','+str(j)+','+str(k)]
                        curr_type = sim_grid[i,j,k]
                        k_m_num = ''
                        k_m_type = ''
                        k_p_num = ''
                        k_p_type = ''
                        if k < bond_size:
                            k_m = top_bound-bond_size-k 
                            k_p = k+bond_size
                        elif k > top_bound-bond_size-1:
                            k_m = k-bond_size
                            k_p = k+bond_size-top_bound
                        else:
                            k_m = k-bond_size
                            k_p = k+bond_size
                        if str(i)+','+str(j)+','+str(k_m) in inv_tot_atom_dict:
                            k_m_num = inv_tot_atom_dict[str(i)+','+str(j)+','+str(k_m)]
                            k_m_type = sim_grid[i,j,k_m]
                        if str(i)+','+str(j)+','+str(k_p) in inv_tot_atom_dict:
                            k_p_num = inv_tot_atom_dict[str(i)+','+str(j)+','+str(k_p)]
                            k_p_type = sim_grid[i,j,k_p]
                    #Add angles to dictionary
                    #Form all angle potentials along linear chain segments
                        if curr_type in chain_list and k_m_type in chain_list and k_p_type in chain_list:
                            lin_angle_dict[angle_count] = str(k_m_num)+','+str(curr_num)+','+str(k_p_num)
                            angle_count += 1
                    #Form all angle potential at branch points
                    else:
                        continue
    return lin_angle_dict,per_angle_dict,angle_count
    
def merge_two_dicts(x,y):
    z = x.copy()
    z.update(y)
    return z
