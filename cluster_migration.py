# %%
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.pyplot as plt
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

# %%
try:
    os.mkdir("figures")
except OSError as error:
    print(error)
try:
    os.mkdir("post_process_files")
except OSError as error:
    print(error)

# %%
def countX(lst, x):
    if lst.count(x):
        return lst.count(x)
    else:
        return 1

# %% [markdown]
# ## Modified Algorithm:

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
#fig = plt.figure()
#ax = fig.add_axes([0.15,0.15,0.80,0.80])
#ax.plot(elongation)
#plt.xlabel("Frame")
#plt.ylabel("elongation instances")
#plt.savefig('figures/elongation_inst.png',dpi = 600)

# %%
#fig = plt.figure()
#ax = fig.add_axes([0.15,0.15,0.80,0.80])
#ax.plot(fragmentation)
#plt.xlabel("Frame")
#plt.ylabel("fragmentation instances")
#plt.savefig('figures/fragmentation_inst.png',dpi = 600)

# %%
#fig = plt.figure()
#ax = fig.add_axes([0.15,0.15,0.80,0.80])
#ax.plot(depoly)
#plt.xlabel("Frame")
#plt.ylabel("depolymerization instances")
#plt.savefig('figures/depoly_inst.png',dpi = 600)

# %% [markdown]
# # Depolymerization algorithm:
# 1. We have data for cluster ID and the coordination number of each particle at different timesteps.
# &nbsp;
# &nbsp;
# 2. We compare the cluster ID of a particle at the current time step to the cluster ID of the same particle at a previous time step.
# &nbsp;
# &nbsp;
# 3. If the cluster ID is different, it means that this particle belongs to a different chain than what was in the previous time step.
# &nbsp;
# &nbsp;
# 4. We identify the particle ID and then go back to the coordination number data of the particle at the current and the previous step.
# &nbsp;
# &nbsp;
# 5. If the coordination number of that particle:
# &nbsp;
# &nbsp;
# ***
# >        a) changes from 2 to 1    :   the particle belonged in the middle of a chain previously and now blongs to  the end of a chain. This could either be migration between chains or fragmentation.
# 
# >        b) changes from 1 to 2    :   the particle was at the end of the chain and now belongs to the middle of a > chain. This could be because of elongation also.
#         
# >        c) changes from 1 to 0    :   the particle was at the end of the chain and now its a free monomer.     therefore if we count the number of such occurences between different timesteps, we have the depolymerization rate.
# ***

# %% [markdown]
# def Block_Average(data,Max_block_Size=2):
#     block_Mean = np.zeros(Max_block_Size)
#     block_Var = np.zeros(Max_block_Size)
#     block_sig = np.zeros(Max_block_Size)
#     block_Err = np.zeros(Max_block_Size)
#     Size = np.arange(1,Max_block_Size+1)
#     for block_Size in Size:
#         Nb = int(len(data)/block_Size)
#         b_mean = np.zeros(Nb)
#         for j in range(1,Nb+1):
#             ibeg = (j-1) * block_Size
#             iend = ibeg + block_Size
#             b_mean[j-1] = np.mean(data[ibeg:iend])
#         block_Mean[block_Size-1] = np.mean(b_mean)
#         block_Var[block_Size-1] = np.var(b_mean)/(Nb-1)
#         block_sig[block_Size-1] = np.sqrt(block_Var[block_Size-1])
#         block_Err[block_Size-1] = block_sig[block_Size-1]/np.sqrt(Nb)
#     return Size,block_Mean,block_sig,block_Err

# %% [markdown]
# ## The code box below is the original algorithm

# %% [markdown]
# 
# frame_no = []
# line_counter = 0
# iter = 0
# part_id = [] #stores the particle id for each frame
# currcluster_id = [] #stores the cluster ID for the current frame
# oldcluster_id=[] #stores the the cluster ID for the previous frame
# curr_coord = [] #stores the coordination number for the current frame 
# old_coord = [] #stores the coordination number for the previous frame
# size = [] #stores the sizes of each cluster ID
# avg_len = [] #average length of the chain
# max_len = [] #maximum length of the chain
# migration = [] #stores the number of times there is a change between the surrent and the previous cluster ID
# curr_particle_cluster_size = [] #stores the chain size that the particle belongs to
# old_particle_cluster_size = []
# dimer = [] #stores the number of dimers in each frame
# trimer = [] #stores the number of trimers in each frame
# hexamer = [] #stores the number of 6-mers in each frame
# elongation = [] #stores elongation instance for each frame
# fragmentation = [] #stores fragmentation instance for each frame
# depoly = [] #stores depolymerization instance for each frame
# f_elon = 0
# f_frag = 0
# f_depo = 0
# iter = 0
# file1 = open("cluster.dat","r")
# file2 = open("coord.dat","r")
# for line,line2 in zip(file1,file2):
#     line_counter = line_counter + 1
#     if line_counter==2:
#         frame_no.append(int(line))
#     if (line_counter>9)&(line_counter <= 1009):
#         f_list = [float(i) for i in line.split()]
#         n_list = [float(i) for i in line2.split()]
#         part_id.append(f_list[0])
#         currcluster_id.append(f_list[1])
#         curr_coord.append(n_list[1])
#     if (line_counter == 1009):
#         part_index = [0]
#         change_count = 0
#         line_counter = 0
#         for i in range(0,len(part_id)):
#             size.append(countX(currcluster_id,part_id[i]))
#             curr_particle_cluster_size.append(countX(currcluster_id,currcluster_id[i])) 
#         avg_len.append(np.mean(size))
#         max_len.append(np.max(size))
#         dimer.append(countX(curr_particle_cluster_size,2))
#         trimer.append(countX(curr_particle_cluster_size,3))
#         hexamer.append(countX(curr_particle_cluster_size,6))
#         if iter >= 1:
#             part_index = [0]
#             for j in range(0,len(oldcluster_id)):
#                 if oldcluster_id[j]!=currcluster_id[j]:
#                     change_count += 1
#             f_elon = 0
#             f_frag = 0
#             f_depo = 0
#             for k in range(0,len(curr_particle_cluster_size)):
#                 if part_id[k] not in part_index:
#                     part_index.append(part_id[k])
#                     if curr_particle_cluster_size[k] > old_particle_cluster_size[k]:
#                         f_elon += 1
#                         for q in range(0,len(curr_particle_cluster_size)):
#                             if (currcluster_id[k]==currcluster_id[q]) & (part_id[q] not in part_index):
#                                 part_index.append(part_id[q])
#                                 
#             part_index = [0]
#             for k in range(0,len(curr_particle_cluster_size)):
#                 if part_id[k] not in part_index:
#                     part_index.append(part_id[k])
#                     if (curr_particle_cluster_size[k] < old_particle_cluster_size[k]) & (curr_coord[k] != 0):
#                         f_frag += 1
#                         for q in range(0,len(curr_particle_cluster_size)):
#                             if (currcluster_id[k]==currcluster_id[q]) & (part_id[q] not in part_index):
#                                 part_index.append(part_id[q])
# 
#             part_index = [0]
#             for k in range(0,len(curr_particle_cluster_size)):
#                 if part_id[k] not in part_index:
#                     part_index.append(part_id[k])
#                     if (curr_particle_cluster_size[k] < old_particle_cluster_size[k]) & (curr_coord[k] == 0):
#                         f_depo += 1
#                         for q in range(0,len(curr_particle_cluster_size)):
#                             if (currcluster_id[k]==currcluster_id[q]) & (part_id[q] not in part_index):
#                                 part_index.append(part_id[q])
#                                                             
#         elongation.append(f_elon)
#         fragmentation.append(f_frag)
#         depoly.append(f_depo)                    
#         migration.append(change_count)
#         oldcluster_id = currcluster_id
#         old_coord = curr_coord
#         old_particle_cluster_size = curr_particle_cluster_size
#         part_id = []
#         currcluster_id = []
#         curr_coord = []
#         size = []
#         curr_particle_cluster_size =[]
#         iter += 1        

# %% [markdown]
# x = [1,1,1,1,1,2,1,1,1,1,1,1,1,3,1,1,1,1,1,1,1,1,1,1,1,4,1,1,1,1,1,1,1,1,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
# print(np.mean(x))
# print(np.average(x,weights=x))

# %% [markdown]
# import pandas as pd
# cdata = pd.read_excel('post_process_files/last_frame_cluster_size.xlsx')
# cdata.head()
# 

# %% [markdown]
# cluster_size = np.arange(1,np.max(cdata.cl_size)+1)
# number_cluster = []
# for i in cluster_size:
#     number_cluster.append((cdata.cl_size.to_numpy().tolist().count(i))/i)

# %% [markdown]
# plt.bar(cluster_size,number_cluster)
# print(np.sum(number_cluster))

# %%



