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
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})
# for Palatino and other serif fonts use:
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
    "font.serif": ["Palatino"],
})
# It's also possible to use the reduced notation by directly setting font.family:
plt.rcParams.update({
  "text.usetex": False,
  "font.family": "Helvetica"
})
import scipy
from scipy.optimize import curve_fit
import math
from MDAnalysis.analysis import polymer
import os

# %%
try:
    os.mkdir("analysis_files")
except OSError as error:
    print(error)
try:
    os.mkdir("analysis_files/chain_pos")
except OSError as error:
    print(error)
try:
    os.mkdir("figures")
except OSError as error:    
    print(error)
try:
    os.mkdir("post_process_files")
except OSError as error:
    print(error)

# %%
u_eq = mda.Universe("data.psf","dump_atom.eq1","dump_atom.eq2",
                    format = 'LAMMPSDUMP')

u_pro = mda.Universe("data.psf","dump_atom.prod",format = 'LAMMPSDUMP')

u_total = mda.Universe("data.psf","dump_atom.eq1","dump_atom.eq2",
                       "dump_atom.prod",format = 'LAMMPSDUMP')

# %%
print("equilibration no of frames  :  ",len(u_eq.trajectory))
print("production no of frames  :  ",len(u_pro.trajectory))

# %%
def block(vector):
    n = len(vector)
    n_prime = math.ceil(int(n/2))
    x_prime = np.zeros(n_prime)
    for i in range(0,n_prime):
        x_prime[i] = 0.5*(vector[2*i] + vector[2*i+1])
    return x_prime,n_prime

def BlockAv(temp_vec):
    n_prime = math.ceil(int(len(temp_vec)/2))
    count = []
    var_test = []
    sig_test = []
    Av = np.mean(temp_vec)
    i=0
    while n_prime>2:
        t_vec,n_prime = block(temp_vec)
        var_test.append(np.var(t_vec)/(n_prime-1))
        sig_test.append(np.std(t_vec)/np.sqrt(n_prime-1))
        temp_vec = t_vec
        i = i+1
        count.append(i)
    plt.plot(count,sig_test,'o')
    plt.ylabel('std-deviation')
    plt.xticks(np.arange(1,count[-1]+1,2))
    plt.xlabel('block operation No.')
    return Av,sig_test[-2]
    #plt.savefig('block.png',dpi = 600)

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
            

# %% [markdown]
# # Key to the algorithm:
# 
# ## **cluster_nr** : 
# ### Number of chains. The final value of the **cluster_nr** for each frame is the number of individual clusters. An individual particle is also called a cluster for simplicity. For example at $t=0$, the number of clusters is the total number of particles itself.
# 
# ## **cluster_status[i]** : 
# ### It is just a number to indicate whether the cluster status of the $i^{th}$ partcile has been checked or not. 0indicated not checked and 1 indicated checked. At the end of the cluster search, all the particles should have a cluster_status of 1.
# 
# ## **cluster_index** : 
# ### Its a vector that stores the cluster identity of the particle. For example, if 
#     cluster_index = [0,0,0,0,1,1,1,1,1,2,2,3,4,4,4] 
# ### then it means that the 4 particles belong to the first cluster, 5particle belong to the 2 cluster, 2 particles belong to the third cluster, 1 particle belongs to the fourth cluster,3 particles belong to the fifth cluster and so on... The length of this vector is equal to the number of particles and the max(cluster_index)+1 gives the total number of clusters for that time frame.
# 
# ## **particle_in_cluster** : 
# ### This vector gives the number of particles that are present in each cluster. For example
#     particle_in_cluster = [16,8,3,2,1,7,1] 
# ### means there are 7 clusters in increasing order of their size and their respective sizes are given by the each element in the vector. A size =1 means that it is an individual particle.
# 
# ## **number_cluster** : 
# ### This vector gives us an idea about the distribution of the cluster sizes. For example
#     number_cluster = [10,3,5,0,6,1] 
# ### means that there are 10 clusters of size 1, 3clusters of size 2, 5 clusters of size 3, 0 clusters of size 4, 6 clusters of size 5 and 1 cluster of size 6 and so on...In detail, the $i^{th}$ element of the cluster gives the number of clusters of size/length $i$
# 
# ## number_of_chains :
# ### This is just to store the number of clusters/chains at each time frame. This is just to see hot the total number of clusters evolve with time.
# 
# ## monomer_number :
# ### This stores the number of free monomers at each time frame. This is to see the kinetics of the consumption of monomers and fit it with a rate eqaution and get a rate constant.
# 
# 
# 

# %% [markdown]
# # Some data to get from lammps input file:

# %%
box = 46.45 #check data.packed

tot_pos = u_total.select_atoms('name B')
eq_pos = u_eq.select_atoms('name B')
pro_pos = u_pro.select_atoms('name B')
pos1 = tot_pos.positions
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
dist = np.zeros_like(np.arange(len(pos1)*len(pos1)).reshape(len(pos1),len(pos1)),dtype = float)    

# %%
#Equilibrium Analysis:
for ts in u_eq.trajectory[0:-1:10]:
    eq_time.append(ts.time*eq_frame_dt*dt)
    pos = eq_pos.positions # this stores the positions of the particles in [Nx3] matrix/array
    for i in range(0,len(pos)):
        for j in range(0,len(pos)):
            if i != j:
                dx = MIC(pos[i,0] - pos[j,0])
                dy = MIC(pos[i,1] - pos[j,1])
                dz = MIC(pos[i,2] - pos[j,2])
                dist[i,j] = math.sqrt(dx*dx + dy*dy + dz*dz)
            if i == j:
                dist[i,j] = 0
    cluster_nr = 0
    
    cluster_status = np.zeros_like(np.arange(len(pos)))
    cluster_index=[]
    for i in range(0,len(pos)):
        if cluster_status[i] == 0:
            cluster_status[i] = 1
            cluster_search(i)
            cluster_index.append(cluster_nr)
            cluster_nr = cluster_nr + 1 
    
    particle_in_cluster = np.zeros_like(np.arange(cluster_nr))

    for i in range(0,len(pos)):
        particle_in_cluster[cluster_index[i]] = particle_in_cluster[cluster_index[i]] + 1
    start_size = 1
    end_size = np.max(particle_in_cluster)
    number_cluster = []
    for i in range(start_size,end_size+1):
        number_cluster.append(particle_in_cluster.tolist().count(i))
    eq_number_of_chains.append(np.sum(number_cluster[2:]))
    eq_monomer_number.append(number_cluster[0])
    eq_max_cluster.append(np.max(particle_in_cluster)) 
    num = []
    for k in range(0,len(number_cluster)):
        num.append(k+1)
    eq_Avg_chain_size.append(np.average(num,weights=number_cluster))

# %%
type(dist)
Min = []
for i in range(0,len(dist)):
    for j in range(0,len(dist)):
        if (abs(dist[i,j]) <= 1.1) & (dist[i,j] > 0.):
            Min.append(abs(dist[i,j]))
            #print(i,j) # will print pairs of particles
plt.plot(Min)
print(len(Min))
print(np.mean(Min),np.std(Min))

# %%
#Production Analysis
for ts in u_pro.trajectory[0:-1:10]:
    pro_time.append(ts.time*pro_frame_dt*dt)
    pos = pro_pos.positions # this stores the positions of the particles in [Nx3] matrix/array
    for i in range(0,len(pos)):
        for j in range(0,len(pos)):
            if i != j:
                dx = MIC(pos[i,0] - pos[j,0])
                dy = MIC(pos[i,1] - pos[j,1])
                dz = MIC(pos[i,2] - pos[j,2])
                dist[i,j] = math.sqrt(dx*dx + dy*dy + dz*dz)
            if i == j:
                dist[i,j] = 0
    cluster_nr = 0
    
    cluster_status = np.zeros_like(np.arange(len(pos)))
    cluster_index=[]
    for i in range(0,len(pos)):
        if cluster_status[i] == 0:
            cluster_status[i] = 1
            cluster_search(i)
            cluster_index.append(cluster_nr)
            cluster_nr = cluster_nr + 1
        
    particle_in_cluster = np.zeros_like(np.arange(cluster_nr))

    for i in range(0,len(pos)):
        particle_in_cluster[cluster_index[i]] = particle_in_cluster[cluster_index[i]] + 1
    start_size = 1
    end_size = np.max(particle_in_cluster)
    number_cluster = []
    for i in range(start_size,end_size+1):
        number_cluster.append(particle_in_cluster.tolist().count(i))
    pro_number_of_chains.append(np.sum(number_cluster[2:]))
    pro_monomer_number.append(number_cluster[0])
    pro_max_cluster.append(np.max(particle_in_cluster))  
    num = []
    for k in range(0,len(number_cluster)):
        num.append(k+1)
    pro_Avg_chain_size.append(np.average(num,weights=number_cluster))

# %%
pro_time = [eq_time[-1]+x for x in pro_time]
total_time = np.concatenate((eq_time,pro_time))
total_monomer_number = np.concatenate((eq_monomer_number,pro_monomer_number))

# %%
av_l,s_l = BlockAv(pro_number_of_chains)
av,s = BlockAv(pro_Avg_chain_size)
max_av,max_s = BlockAv(pro_max_cluster)

# %%
#writing the information to a file
print("average number of clusters (production) : ",av_l, " +- ",s_l)
print(number_cluster)
print(np.sum(particle_in_cluster)) #- to check whether all particles are accounte for
print("weighted average chain length (production): ", av , " +- ", s)
print("The production average maximum cluster size :",max_av," +- ",max_s)

# %%
#writing xhain size distribution to a file:
fileh = open("analysis_files/chain_size_distribution.dat","w")
fileh.write("Chain Size Distribution:\n")
fileh.write("Num of particles    Number of clusters\n")
num = []
for i in range(0,len(number_cluster)):
    fileh.write("%d              %d\n"%(i+1,number_cluster[i]))
    num.append(i+1)
fileh.close()

# %%
fig = plt.figure()
ax = fig.add_axes([0.15,0.15,0.80,0.80])
n, bins, patches = plt.hist(particle_in_cluster,color = 'green',edgecolor = 'black',bins = 28,width = 1)
plt.xlabel("Chain Size (no. of particles)", fontsize = 16)
plt.ylabel("No. Of Chains",fontsize = 16)
#plt.title(r"$ \epsilon = 40,\rho = 0.01 $",fontsize = 16)
plt.xticks(fontsize =14)
plt.yticks(fontsize = 14)
plt.savefig("figures/chain_size_barplot.png",dpi = 600)

# %%
# writing chain statistics to a file:
stat_f = open('analysis_files/chain_stats.dat','w')
print("average number of chains (production)",file = stat_f)
print(av_l, " +- ",s_l,"\n",file = stat_f)
print("weighted average chain length (production):",file = stat_f)
print(av , " +- ", s,"\n",file = stat_f)
print("production average maximum chain size:",file = stat_f)
print(max_av," +- ",max_s,"\n",file = stat_f)
stat_f.close()

# %% [markdown]
# # Monomer concentration wrt time
# The monomer concentration is fitted to a PBE correspondin to A Brownian Coagualation model devised by Soluchowski (Taken from *Smoke, Dust and Haze by Sheldon Freidlander*). The monomer fraction is given by:
# $$ m = \frac{1}{(1+kt)^2}$$
# 

# %%
len(eq_monomer_number)

# %%
def objective_fun(x,k):
    m0 = 1.0
    return m0/(1+k*x)**2

X = eq_time
Y = np.array(eq_monomer_number)/1000
pars, cov = curve_fit(f=objective_fun, xdata=X, ydata=Y, p0=[0])
print(pars)
y = objective_fun(eq_time,pars)

# %%
fig = plt.figure()
ax = fig.add_axes([0.15,0.15,0.80,0.80])
ax.plot(eq_time,Y,'.',markersize = 15,alpha = 0.2)
ax.plot(eq_time,y,'--',linewidth = 2)
#plt.text(0.7e6,200,"$\frac{m_0}{(1+kt)^2}$",fontsize = 14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xlabel('Time',fontsize = 20)
plt.ylabel("Monomer fraction",fontsize = 20)
plt.text(0.6e6, 0.9, '$k=%e$'%(pars), fontsize=14)
plt.savefig("figures/monomer_rate.png",dpi = 600)
r_f = open("analysis_files/rate.dat",'w')
print("monomer depletion rate:",file = r_f)
print(pars[0],file = r_f)
r_f.close()

# %%
#calculating and plotting radial distribution function (g(r))
B_select = u_pro.select_atoms('name B')
irdf = rdf.InterRDF(B_select,B_select,exclusion_block=(1,1),nbins=200)
irdf.run()
fig = plt.figure()
ax = fig.add_axes([0.15,0.15,0.78,0.8])
ax.plot(irdf.bins, irdf.rdf,linewidth = 2)
plt.xlabel('$r$ ($\sigma$)',fontsize = 20)
plt.ylabel('$g(r)$',fontsize = 20)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.savefig("figures/rdf.png",dpi = 600)

# %%
# Calculating coordination number from radial distribution function
rdf_r = irdf.bins;
rdf_g = irdf.rdf;
plt.plot(rdf_r[:16],rdf_g[:16])
CN = 0
number_density = 1000/(box**3) #check equilibration.dat
mass_density = 0.01 #check equilibration.dat
for i in range(0,16):
    CN = CN + (4*3.14*number_density*rdf_g[i]*(rdf_r[i]**2))*(rdf_r[i+1]-rdf_r[i])
print(CN)  

# %%
# calculating coordination number from average number of nearest neighbours
def contacts_within_cutoff(u, group_a, group_b, radius=1.1):
    timeseries = []
    for ts in u.trajectory[-10:]:
        # calculate distances between group_a and group_b
        d = contacts.distance_array(group_a.positions, group_b.positions)
        # determine which distances <= radius
        n_contacts = contacts.contact_matrix(d, radius).sum()
        timeseries.append([ts.frame, n_contacts])
    return np.array(timeseries)

# %%
ca = contacts_within_cutoff(u_pro,u_pro.select_atoms('name B'),u_pro.select_atoms('name B'), radius=1.1)
neigh = (ca[:,1]-1000)/1000
av_neigh,sig_neigh = BlockAv(neigh)
print("coordination number from average neighbour:")
print(av_neigh, "+-", sig_neigh)
print("coordination number from g(r):")
print(CN)
f = open('analysis_files/coordination.dat','w')
print("coordination number from average neighbour:",file = f)
print(av_neigh, "+-", sig_neigh,"\n",file = f)
print("coordination number from g(r):",file = f)
print(CN,file =f)
f.close()

# %% [markdown]
# # Vector Plot

# %%
A = u_pro.select_atoms('name A')
B = u_pro.select_atoms('name B')
A1 = u_pro.atoms[list(range(0,3000,3))]
for ts in u_pro.trajectory[-1]:
    posB = B.positions
    posA = A1.positions
ax = plt.figure(figsize = (10,10)).add_subplot(projection='3d')
x = posB[:,0]
y = posB[:,1]
z = posB[:,2]
u = posB[:,0]-posA[:,0]
v = posB[:,1]-posA[:,1]
w = posB[:,2]-posA[:,2]

ax.quiver(x,y,z,u,v,w,normalize=True,color = "red",length = 1.5)
ax.set_xlabel("x",fontsize = 20)
ax.set_ylabel("y",fontsize=20)
ax.set_zlabel("z",fontsize = 20)
plt.savefig("figures/quiver.png",dpi = 600)


# %% [markdown]
# # Chain Orientation Factor:

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
for ts in u_pro.trajectory[-100:]:
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

# %%
print(np.mean(Q_m))

# %%
def countX(lst, x):
    if lst.count(x):
        return lst.count(x)
    else:
        return 1

# %%
frame_no = []
line_counter = 0
iter = 0
part_id = [] #stores the particle id for each frame
currcluster_id = [] #stores the cluster ID for the current frame
oldcluster_id=[] #stores the the cluster ID for the previous frame
curr_coord = [] #stores the coordination number for the current frame 
old_coord = [] #stores the coordination number for the previous frame
size = [] #stores the sizes of each cluster ID
avg_len = [] #average length of the chain
max_len = [] #maximum length of the chain
migration = [] #stores the number of times there is a change between the surrent and the previous cluster ID
curr_particle_cluster_size = [] #stores the chain size that the particle belongs to
old_particle_cluster_size = []
dimer = [] #stores the number of dimers in each frame
trimer = [] #stores the number of trimers in each frame
hexamer = [] #stores the number of 6-mers in each frame
elongation = [] #stores elongation instance for each frame
fragmentation = [] #stores fragmentation instance for each frame
depoly = [] #stores depolymerization instance for each frame
f_elon = 0
f_frag = 0
f_depo = 0
iter = 0
file1 = open("cluster.dat","r")
file2 = open("coord.dat","r")
c_file = open("post_process_files/cluster_dynamic_stat.dat","w")
print("time_frame   elongation  fragmentation   depolymerization",file = c_file)
p_file = open("post_process_files/particle_cluster_size.dat","w")
print("Size of the cluster each particle ID belongs to",file = p_file)
m_file = open("analysis_files/monomer_conc.dat","w")
print("time_frame   Monomer_No",file = m_file)
for line,line2 in zip(file1,file2):
    line_counter = line_counter + 1
    if line_counter==2:
        frame_no.append(int(line))
        print(int(line),file = p_file)
    if (line_counter>9)&(line_counter <= 1009):
        f_list = [float(i) for i in line.split()]
        n_list = [float(i) for i in line2.split()]
        part_id.append(f_list[0])
        currcluster_id.append(f_list[1])
        curr_coord.append(n_list[1])
    if (line_counter == 1009):
        part_index = [0]
        change_count = 0
        line_counter = 0
        for i in range(0,len(part_id)):
            size.append(countX(currcluster_id,part_id[i]))
            curr_particle_cluster_size.append(countX(currcluster_id,currcluster_id[i]))
            print(int(part_id[i]),"  ",curr_particle_cluster_size[i], file = p_file) 
        avg_len.append(np.average(size,weights = size)) #weighted average
        #avg_len.append(np.mean(size))  #regular mean
        max_len.append(np.max(size))
        print(frame_no[iter],"   ",countX(curr_particle_cluster_size,1),file = m_file)
        dimer.append(countX(curr_particle_cluster_size,2)/2)
        trimer.append(countX(curr_particle_cluster_size,3)/3)
        hexamer.append(countX(curr_particle_cluster_size,6)/6)
        if iter >= 1:
            part_index = [0]
            for j in range(0,len(oldcluster_id)):
                if oldcluster_id[j]!=currcluster_id[j]:
                    change_count += 1
            f_elon = 0
            f_frag = 0
            f_depo = 0
            for k in range(0,len(curr_particle_cluster_size)):
                if part_id[k] not in part_index:
                    part_index.append(part_id[k])
                    if curr_particle_cluster_size[k] > old_particle_cluster_size[k]:
                        f_elon += 1
                        for q in range(0,len(curr_particle_cluster_size)):
                            if (currcluster_id[k]==currcluster_id[q]) & (part_id[q] not in part_index):
                                part_index.append(part_id[q])
                    elif (curr_particle_cluster_size[k] < old_particle_cluster_size[k]) & (curr_coord[k] != 0):# & (old_particle_cluster_size[k] > 3):
                        f_frag += 1
                        for q in range(0,len(curr_particle_cluster_size)):
                            if (currcluster_id[k]==currcluster_id[q]) & (part_id[q] not in part_index):
                                part_index.append(part_id[q])
                    elif (curr_particle_cluster_size[k] < old_particle_cluster_size[k]) & (curr_coord[k] == 0):# & (old_particle_cluster_size[k] > 2):
                        f_depo += 1
                        for q in range(0,len(curr_particle_cluster_size)):
                            if (currcluster_id[k]==currcluster_id[q]) & (part_id[q] not in part_index):
                                part_index.append(part_id[q])
                    
                                                
        elongation.append(f_elon)
        fragmentation.append(f_frag)
        depoly.append(f_depo)                    
        migration.append(change_count)
        oldcluster_id = currcluster_id
        old_coord = curr_coord
        old_particle_cluster_size = curr_particle_cluster_size
        part_id = []
        currcluster_id = []
        curr_coord = []
        size = []
        curr_particle_cluster_size =[]
        iter += 1
for j in range(0,len(frame_no)):
    print(frame_no[j]," ",elongation[j],"   ",fragmentation[j],"    ",depoly[j],file=c_file)
c_file.close()
p_file.close()
m_file.close()
print(part_index)        

# %%
fig = plt.figure()
ax = fig.add_axes([0.15,0.15,0.80,0.80])
ax.plot(elongation)
plt.xlabel("Frame")
plt.ylabel("elongation instances")
plt.savefig('figures/elongation_inst.png',dpi = 600)

# %%
fig = plt.figure()
ax = fig.add_axes([0.15,0.15,0.80,0.80])
ax.plot(fragmentation)
plt.xlabel("Frame")
plt.ylabel("fragmentation instances")
plt.savefig('figures/fragmentation_inst.png',dpi = 600)

# %%
fig = plt.figure()
ax = fig.add_axes([0.15,0.15,0.80,0.80])
ax.plot(depoly)
plt.xlabel("Frame")
plt.ylabel("depolymerization instances")
plt.savefig('figures/depoly_inst.png',dpi = 600)


