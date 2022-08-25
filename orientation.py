# %%
import MDAnalysis as mda 
import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sn 
import math
import pandas as pd 
from MDAnalysis.analysis import contacts
from MDAnalysis.analysis import rdf
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})
# for Palatino and other serif fonts use:
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
})
# It's also possible to use the reduced notation by directly setting font.family:
plt.rcParams.update({
  "text.usetex": True,
  "font.family": "Helvetica"
})
import scipy
from scipy.optimize import curve_fit
import math
from MDAnalysis.analysis import polymer
import os

# %%
u_pro = mda.Universe("data.psf","dump_atom.prod",format = 'LAMMPSDUMP')

# %%
#Function for the Minimum Image Convention
def MIC(d):
    global box
    if d > box/2.:
        d = d-box
    elif d < -box/2.:
        d = d + box
    return d

# %%
def cluster_search(j):
    global dist,pos,cluster_status,cluster_index,cluster_nr,d_cut,monomer
    for i in range(0,len(pos)):
        if (i != j) & (abs(dist[i,j]) < d_cut) & (cluster_status[i] == 0):
            cluster_status[i] = 1
            cluster_index.append(cluster_nr)
            cluster_search(i)
            

# %%
box = 46.45 #check data.packed

#tot_pos = u_total.select_atoms('name B')
#eq_pos = u_eq.select_atoms('name B')
pro_pos = u_pro.select_atoms('name B')
#pos1 = tot_pos.positions
eq_Avg_chain_size = []
pro_Avg_chain_size = []
eq_number_of_chains = []
pro_number_of_chains = []
eq_time = []
pro_time = []
eq_max_cluster = []
pro_max_cluster = []
eq_frame_dt = 1000000  #check input.lmp
pro_frame_dt = 100000 #check input.lmp
dt = 0.001 #check input.lmp
eq_monomer_number = []
pro_monomer_number = []
d_cut = 1.1
dist = np.zeros_like(np.arange(len(pro_pos)*len(pro_pos)).reshape(len(pro_pos),len(pro_pos)),dtype = float)    

# %%
def cluster_search_new(j):
    global dist,pos,cluster_status,cluster_index,cluster_nr,d_cut,monomer,p1,p2,posA1,posA2,xyz_file
    for i in range(0,len(pos)):
        if (i != j) & (abs(dist[i,j]) < d_cut) & (cluster_status[i] == 0):
            cluster_status[i] = 1
            p1.append(posA1[i])
            p2.append(posA2[i])
            write_pos(posA1[i],xyz_file)
            cluster_index.append(cluster_nr)
            cluster_search_new(i)           

# %%
def unit_vector(vec):
    mag = np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
    return np.array(vec)/mag

# %%
def magnitude(vec):
    return(np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2))

# %%
def av_angle(pos_A1,pos_A2):
    cos_t = 0
    brack = 0
    Np = 0
    for i in range(1,len(pos_A1)):
        u1 = unit_vector(np.array([(pos_A1[i-1][0]-pos_A2[i-1][0]),(pos_A1[i-1][1]-pos_A2[i-1][1]),(pos_A1[i-1][2]-pos_A2[i-1][2])]))
        u2 = unit_vector(np.array([(pos_A1[i][0]-pos_A2[i][0]),(pos_A1[i][1]-pos_A2[i][1]),(pos_A1[i][2]-pos_A2[i][2])]))
        cos_t = np.dot(u1,u2)
        brack = brack + (3*(cos_t)**2 - 1)
        Np += 1
   
    return brack/(2*Np)

# %%
def write_pos(vec,f_handle):
    print("H ",vec[0]," ",vec[1]," ",vec[2],file = f_handle)

# %%
B = u_pro.select_atoms('name B')
A1 = u_pro.atoms[list(range(0,3000,3))]
A2 = u_pro.atoms[list(range(2,3000,3))]
pos = pro_pos.positions
for ts in u_pro.trajectory[-300:]:
    Q_m = []
    posB = B.positions
    posA1 = A1.positions
    posA2 = A2.positions
    for i in range(0,len(pos)):
        for j in range(0,len(pos)):
            if i != j:
                dx = MIC(posB[i,0] - posB[j,0])
                dy = MIC(posB[i,1] - posB[j,1])
                dz = MIC(posB[i,2] - posB[j,2])
                dist[i,j] = math.sqrt(dx*dx + dy*dy + dz*dz)
            if i == j:
                dist[i,j] = 0
    cluster_nr = 0
    cluster_status = np.zeros_like(np.arange(len(pos)))
    cluster_index=[]
    p1 = []
    p2 = []
    Q = []
    file_counter = 0
    xyz_file = open("analysis_files/chain_pos/chain%d.xyz"%file_counter,'w')
    for i in range(0,len(posB)):
        if cluster_status[i] == 0:
            cluster_status[i] = 1
            p1.append(posA1[i])
            p2.append(posA2[i])
            write_pos(posA1[i],xyz_file)
            cluster_search_new(i)
            if len(p1)>1:
                Q.append(av_angle(p1,p2))
            cluster_index.append(cluster_nr)
            cluster_nr = cluster_nr + 1
            xyz_file.close()
            file_counter += 1
            xyz_file = open("analysis_files/chain_pos/chain%d.xyz"%file_counter,'w')
            p1 = []
            p2 = []
    Q_m.append(np.mean(Q))
xyz_file.close()

# %%
qfile = open('analysis_files/orientation_factor.dat','w')
print("orientation factor:",file = qfile)
print(np.mean(Q_m)," +- ",np.std(Q_m),file = qfile)
qfile.close()
print(np.mean(Q_m),np.std(Q_m))


