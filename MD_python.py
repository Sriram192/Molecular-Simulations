import numpy as np

# Simulation parameters
num_particles = 100
num_steps = 1000
time_step = 0.001
box_size = 10.0
cutoff = 2.5
verlet_cutoff = cutoff + 0.2

epsilon = 1.0
sigma = 1.0

target_temperature = 1.0
tau = 0.1

pressure = 1.0
target_density = num_particles / box_size**3
target_volume = num_particles / target_density
box_volume = box_size**3
beta = 1.0 / target_temperature
beta_p = beta * pressure

# Initialize particle positions in a simple cubic lattice
num_particles_per_dim = int(np.ceil(num_particles**(1/3)))
spacing = box_size / num_particles_per_dim
positions = np.zeros((num_particles, 3))
index = 0
for x in range(num_particles_per_dim):
    for y in range(num_particles_per_dim):
        for z in range(num_particles_per_dim):
            positions[index] = [x * spacing, y * spacing, z * spacing]
            index += 1
            if index >= num_particles:
                break
        if index >= num_particles:
            break
    if index >= num_particles:
        break

# Initialize particle velocities randomly
velocities = np.random.rand(num_particles, 3)

# Initialize the Verlet list
verlet_list = [[] for _ in range(num_particles)]

# Initialize arrays to store simulation data
potential_energy_array = np.zeros(num_steps)
temperature_array = np.zeros(num_steps)
pressure_array = np.zeros(num_steps)

# Main simulation loop
for step in range(num_steps):
    forces = np.zeros((num_particles, 3))
    potential_energy = 0.0
    
    # Update the Verlet list every few steps (this is a simple approach)
    if step % 10 == 0:
        verlet_list = [[] for _ in range(num_particles)]
        for i in range(num_particles):
            for j in range(i + 1, num_particles):
                r_ij = positions[j] - positions[i]
                r_ij = r_ij - np.round(r_ij / box_size) * box_size
                r_sq = np.sum(r_ij**2)
                if r_sq < verlet_cutoff**2:
                    verlet_list[i].append(j)
                    verlet_list[j].append(i)
    
    # Calculate forces and potential energy for LJ interactions
    for i in range(num_particles):
        for j in verlet_list[i]:
            r_ij = positions[j] - positions[i]
            r_ij = r_ij - np.round(r_ij / box_size) * box_size
            r_sq = np.sum(r_ij**2)
            if r_sq < cutoff**2:
                r = np.sqrt(r_sq)
                lj_energy = 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)
                lj_force = 48 * epsilon * ((sigma**12 / r**13) - (sigma**6 / r**7))
                forces[i] += lj_force * r_ij / r
                forces[j] -= lj_force * r_ij / r
                potential_energy += lj_energy / 2
    
    # Calculate current temperature and update velocities using thermostat
    current_temperature = np.sum(np.sum(velocities**2, axis=1)) / (3 * num_particles)
    scaling_factor = np.sqrt(1 + (time_step / tau) * (target_temperature / current_temperature - 1))
    velocities *= scaling_factor
    
    # Apply volume scaling to maintain pressure using the barostat
    volume_scaling_factor = np.exp(beta_p * (target_volume - box_volume))
    positions *= volume_scaling_factor**(1/3)
    box_size = box_size * volume_scaling_factor**(1/3)
    box_volume = box_size**3
    
    # Update positions and velocities using Verlet integration
    positions += velocities * time_step + 0.5 * forces * time_step**2
    velocities += 0.5 * forces * time_step

    # Apply periodic boundary conditions
    positions %= box_size

    # Store simulation data for each step
    potential_energy_array[step] = potential_energy
    temperature_array[step] = current_temperature
    pressure_array[step] = pressure + (1 / beta_p) * (1 - box_volume / target_volume)

    # Print simulation progress
    if step % 100 == 0:
        print(f"Step: {step}, Potential Energy: {potential_energy}, Current Temperature: {current_temperature}, Current Pressure: {pressure_array[step]}")

print("Simulation complete.")
