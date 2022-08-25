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
def countX(lst, x):
    if lst.count(x):
        return lst.count(x)
    else:
        return 1

# %%
def countN(lst, x):
    if lst.count(x):
        return lst.count(x)
    else:
        return 0

# %%

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

file1 = open("cluster.dat","r")
f_file = open("analysis_files/chain_size_distribution_new.dat","w")
print("Number   Length",file = f_file)
for line in file1.readlines()[25223991:]:
    line_counter += 1
    if line_counter==2:
        frame_no = int(line)
    if (line_counter>9):
        f_list = [int(i) for i in line.split()]
        part_id.append(f_list[0])
        currcluster_id.append(f_list[1])
    
for i in range(0,len(part_id)):
    size.append(countX(currcluster_id,part_id[i]))
    curr_particle_cluster_size.append(countX(currcluster_id,currcluster_id[i]))
np.savetxt("temp.txt",curr_particle_cluster_size,fmt = "%d")
sum = 0    
for q in range(1,np.max(curr_particle_cluster_size)+1):
    print(q,"   ",countN(curr_particle_cluster_size,q)/q,file = f_file)
    sum += countN(curr_particle_cluster_size,q)
      
f_file.close()
file1.close()
print(sum)

# %%
data = np.genfromtxt("analysis_files/monomer_conc.dat",skip_header=1)
data[:,1] = (1000-data[:,1])/1000
np.savetxt("temp.txt",data,fmt = "%0.3f")

# %%



