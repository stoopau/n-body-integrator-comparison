import numpy as np
import matplotlib.pyplot as plt 

#Order: Sun, Mercury, Venus, Earth, Mars, Jupiter

#Physics equations used for this simulation. All manually derived. : 
#   Velocity: V = V + a*dt
#   Force of Gravity: F = G (m1m2 / r^2)
#   Pairwise difference: dx = x2-x1, dy = y2-y1        
#   Distance between two bodies (pythag): r = sqrt(dx^2 + dy^2)
#   Unit vector (pure direc, mag=1): (dx/r, dy/r)
#   Force vector components: Fx = G m1m2/r^3, Fy = G m1m2/r^3
#   Newtons second law: ax = Fx/m, ay = Fy/m 
#   Velocity update (Euler integration): new = old + a * dt 
#   Position update (Euler integration): new = old + new * dt
#   Some slight energy drift will be seen over longer durations bc Euler integration. 
#   Softening term to avoid division by zero. r_safe = r + epsilon

#Goal: 
#   Upgrade Euler integration to rk4 or leapfrog integration. Done. 


#   Gravitational constant 
G = 6.6743e-11

#   Simulated time passed (seconds). (How much time passes per step)
#   3600/hr, 86400/day, 9.467e7/3 years
#   This multiplied with steps = total time elapsed. 

dt = 86400

#Orbiting Distances. Starting with Sun at (0,0) in Meters. In order. 
x = np.array([0, 5.79e10, 1.08e11, 1.496e11, 2.28e11, 7.78e11])
y = np.array([0,0,0,0,0,0])

#Each Velocity. (Meters/Second)
velocity_x = ([0,0,0,0,0,0])
velocity_y = ([0, 47900, 35000, 29780, 24100, 13070])

#Create the mass array. In order.
mass = np.array([1.989e30, 3.30e23, 4.87e24, 5.97e24, 6.42e23, 1.90e27])


#Loop to simulate. 
#Prepare a 2D numpy array of shape (num_steps (365), num_bodies(6)) for history. 
#Faster than a list where we wait for it to create new slots for appending. No waiting. Slots exist and are empty and ready. 
#History will reflect the x and y histories of each entity, so we can see their orbit path. 

position_old_x = 0
position_old_y = 0 

#   365 steps/year, 1095/3 year
num_steps = 73000
num_bodies = 6

#Empty histories
x_history = np.zeros((num_steps, len(x)))
y_history = np.zeros((num_steps, len(y)))

PE_history = np.zeros(num_steps)
KE_history = np.zeros(num_steps)

#Define energy function for finding potential/kinetic energies. 
def compute_energy(x, y, velocity_x, velocity_y, mass): 
    
    #Pairwise differences. Dx and Dy. Again.
    col_x = x.reshape(-1,1)
    dx = x - col_x
    
    col_y = y.reshape(-1,1)
    dy = y - col_y 

    #Update distance matrix. Safe guard r. Again.
    r = np.sqrt(dx**2 + dy**2)
    r_epsilon = 0.0000000000000001 
    r_safe = r + r_epsilon 
    
    #Mass Matrix. m1m2. Again.
    col_mass = mass.reshape(-1,1)
    m1m2 = col_mass * mass 
    
    #Compute potential energy matrix. 
    PE_matrix = -G * (m1m2/r_safe)
    
    #Cant np.sum the matrix, need conversion. 
    #Using index of i and index of j 
    #np.triu_indices gives row/col indices of upper triangle. 
    #Unique pairs. no double count. 
    #Compute PE from matrix. 
    
    i_index, j_index = np.triu_indices (len(mass), k=1)
    PE = np.sum(PE_matrix[i_index, j_index])
    
    #Compute kinetic energy. 
    KE = np.sum(0.5 * mass * (velocity_x**2 + velocity_y**2))
    
    return PE, KE

for step in range(num_steps): 
    
    x_history[step] = x
    y_history[step] = y
    
    #Pairwise differences. Dx and Dy.
    col_x = x.reshape(-1,1)
    dx = x - col_x
    
    #print("Dx Matrix:")
    #print(dx)
    #print()
    
    col_y = y.reshape(-1,1)
    dy = y - col_y 
    
    #print("Dy Matrix:")
    #print(dy)
    #print()

    #Update distance matrix. Safe guard r. 
    r = np.sqrt(dx**2 + dy**2)
    r_epsilon = 0.0000000000000001 
    r_safe = r + r_epsilon 
    
    #Mass Matrix. m1m2. 
    col_mass = mass.reshape(-1,1)
    m1m2 = col_mass * mass 
    
    #print("M1*M2 Matrix:")
    #print()
    #print(m1m2)
    #print()
    
    #Force Components  
    Fx = G * (m1m2*dx/r_safe**3)
    Fy = G * (m1m2*dy/r_safe**3)
    
    if step == 0:
        print("Sun x:", x[0], "Earth x:", x[3])
        print("Fx on Earth from Sun (Fx[3][0]):", Fx[3][0])
        
    #Total Force
    total_Fx = np.sum(Fx, axis=1)
    total_Fy = np.sum(Fy, axis=1)
    
    #Acceleration
    acceleration_x = total_Fx / mass 
    acceleration_y = total_Fy / mass 
    
    #print("Acceleration X: ")
    #print(acceleration_x)
    #print()
    #print("Acceleration Y: ")
    #print (acceleration_y)
    #print()
    
    #Update velocity 
    velocity_x = velocity_x + acceleration_x * dt
    velocity_y = velocity_y + acceleration_y * dt 
    
    #Update position 
    x = x + velocity_x * dt 
    y = y + velocity_y * dt
        
    #Store the history 
    x_history[step] = x
    y_history[step] = y 
    
    #Extract energies. 
    PE, KE = compute_energy(x, y, velocity_x, velocity_y, mass)
    
    #Store energy history
    KE_history[step] = KE
    PE_history[step] = PE
    
    #I need to visualize it. I dont wanna look at matricies and data. \
    #Need the Sun's history. 
    #Need the Earth's history. 
    #Need the other bodies. 
    # (:, y) = all rows. 
    
    sun_x = x_history[:,0] 
    sun_y = y_history[:,0]
    
    mercury_x = x_history[:,1]
    mercury_y = y_history[:,1]
    
    venus_x = x_history[:,2]
    venus_y = y_history[:,2]
    
    earth_x = x_history[:,3]
    earth_y = y_history[:,3] 
    
    mars_x = x_history[:,4]
    mars_y = y_history[:,4]
    
    jupiter_x = x_history[:,5]
    jupiter_y = y_history[:,5]
    

#Use matplotlib to visually plot the history and show its path. 
plt.plot(sun_x, sun_y, label="Sun")
plt.plot(mercury_x, mercury_y, label="Mercury")
plt.plot(venus_x, venus_y, label="Venus")
plt.plot(earth_x, earth_y, label="Earth")
plt.plot(mars_x, mars_y, label="Mars")
plt.plot(jupiter_x, jupiter_y, label="Jupiter")
    
plt.legend()

#keep x and y axis equal for circular orbit. 
plt.axis('equal')

#Save orbital results
plt.savefig("results/euler_orbits.png")
plt.show()

#Inspect energy drift. 
total_energy = KE_history + PE_history 
plt.plot(total_energy / total_energy[0]-1)

#Save energy results. 
plt.savefig("results/euler_energy_drift.png")
plt.show() 





