import numpy as np
import random as ran
#import scipy.constants as con
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.cm as cm

class Par:
    def __init__(self,x_position,y_position,x_velocity,y_velocity,charge,mass):
        self.x = x_position
        self.y = y_position
        self.vx = x_velocity
        self.vy = y_velocity
        self.q = charge
        self.m = mass

class Grid_Info:
    def __init__(self,x_grid_position,y_grid_position,weight_ij,weight_i1j,weight_i1j1,weight_ij1):
        self.x_i = x_grid_position
        self.y_i = y_grid_position
        self.f_ij = weight_ij
        self.f_i1j = weight_i1j
        self.f_i1j1 = weight_i1j1
        self.f_ij1 = weight_ij1

class Physical_Constants:
    #CGS
    def __init__(self):
        self.c = 3e10  #Speed of light
        self.e = 4.80320427*10**-10 #Elementary charge
        self.m_e = 9.10938215*10**-28 #Electron Mass
        self.pi = 3.14159265359 #pi

con = Physical_Constants()

#Spatial domain variables
x_min = 0 #Minimum position in x domain
x_max = 50 #Maximum position in x domain
y_min = 0 #Minimum position in y domain
y_max = 50 #Maximum position in y domain

nx = 100 #Number of steps taken from x_min to x_max
ny = 100#Number of steps taken from y_min to y_max

dx = (x_max - x_min) / nx #Size of step in x domain
dy = (y_max - y_min) / ny #Size of step in y domain
T_c = nx * ny #Total number of cells

v_min = -(dx/10) #Minimum initial particle speed
v_max = (dx/10) #Maximum initial particle speed

#Populate dimesnional vectors
x = np.linspace(x_min,x_max,num=(nx+1))
y = np.linspace(y_min,y_max,num=(ny+1))

#Meshgrid for plotting
X,Y = np.meshgrid((x),(y))

#Electromagentic field variables
phi = np.zeros(((nx + 1),(ny + 1))) #Electric potential
E_x_n_mhalf = np.zeros(((nx + 1),(ny + 1)))
E_y_n_mhalf = np.zeros(((nx + 1),(ny + 1)))
E_x_n_phalf = np.zeros(((nx + 1),(ny + 1)))
E_y_n_phalf = np.zeros(((nx + 1),(ny + 1)))
J_x = np.zeros(((nx + 1),(ny + 1)))
J_y = np.zeros(((nx + 1),(ny + 1)))
B_z_n = np.zeros(((nx + 1),(ny + 1)))
B_z_n_pone = np.zeros(((nx + 1),(ny + 1)))

#Time Variables
t_min = 0 #Start time
t_max = 2.0 #End time
nt = 200 #Number of time steps
t = np.linspace(t_min,t_max,num=nt) #Populate dimesnional vectors
dt = (abs(t_min) + t_max) / nt #Size of step on the temporal domain

#Particle charge & mass lists
#q_c = [(10*con.e),-(100*con.e)] #particle charge
#m = [(10*con.m_p),(100*con.m_e)] #particle mass
q_c = [-con.e] #particle charge
m = [con.m_e] #particle mass

n_smoothing = 50 #For differencing methods that require smoothing
PPC = .1 #Number of particles per cell

T_PPC = int(T_c * PPC) #Total number of particles in simulation

#Create list of particles properties
Par_x = np.zeros(T_PPC)
Par_y = np.zeros(T_PPC)
Par_vx = np.zeros(T_PPC)
Par_vy = np.zeros(T_PPC)
Par_q = np.zeros(T_PPC)
Par_m = np.zeros(T_PPC)
Par_species_count = np.zeros(len(q_c))

p = []

for i in range(0,T_PPC):

    #TODO: Generalise
    j = ran.randint(0,int(len(q_c)-1)) #Random integer for mass and charge selection
    s = ran.choice((-1,1))

    Par_x = ran.uniform(x_min,x_max) #Particle x coordinate
    Par_y = ran.uniform(y_min,y_max) #Particle y coordinate
    Par_vx = ran.uniform(v_min,v_max) #Particle velocity on the x direction
    Par_vy = ran.uniform(v_min,v_max) #Particle velocity on the y direction
    Par_q = q_c[j] #Particle charge
    Par_m = m[j] #Particle mass

    p.append(Par(Par_x,Par_y,Par_vx,Par_vy,Par_q,Par_m))

print('Initialising particle simulation...')
##########################Calculate initial conditions##########################
print('Calculating Initial Conditions...')
#B_z_n = (((X**2)+(Y**2))**(3/2))*10e-12
B_z_n.fill(1)
q_grid = np.zeros(((nx + 1),(ny + 1))) #Charge density
for i in range(0,T_PPC):
    x_i = 0 #Particle x grid position
    y_i = 0 #Particle y grid position
    hx = 0 #Partcle distance from its grid postion in the x domain
    hy = 0 #Partcle distance from its grid postion in the y domain

    x_i = int(np.floor(p[i].x / dx))
    y_i = int(np.floor(p[i].y / dy))

    hx = p[i].x - x[x_i]
    hy = p[i].y - y[y_i]

    q_grid[x_i][y_i] += (p[i].q/(dx*dy))*(((dx - hx) * (dy - hy)) / (dx * dy))
    q_grid[x_i + 1][y_i] += (p[i].q/(dx*dy))*((hx * (dy - hy)) / (dx * dy))
    q_grid[x_i + 1][y_i + 1] += (p[i].q/(dx*dy))*((hx * hy) / (dx * dy))
    q_grid[x_i][y_i + 1] += (p[i].q/(dx*dy))*(((dx - hx)*hy) / (dx*dy))

phi = (((q_grid)*((dx**2)*(dy**2)) + ((dy**2)*(np.roll(phi,1,axis=0) +
np.roll(phi,-1,axis=0))) + ((dx**2)*(np.roll(phi,1,axis=1) +
np.roll(phi,-1,axis=1)))) / (2*((dx**2) + (dy**2))))

for n in range(0, n_smoothing):
    E_x_n_mhalf = ((np.roll(phi,-1,axis=0)) - np.roll(phi,1,axis=0)) / (2 * dx)
    E_y_n_mhalf = ((np.roll(phi,-1,axis=1)) - np.roll(phi,1,axis=1)) / (2 * dy)
################################################################################

################################Main Loop######################################
print('Initialising Main Loop...')
#Numerical loop responsible for interating the particles through time
for t in range(0,10): #nt

    plt.clf()

    #Print index 't' every 10 time steps to show progress
    #if t % 10 == 0:

    print(t)

    #Perpare index 'i' for while loop responsible for removing particles which
    #move outside of the simulation region
    i = 0
    g = []
    #While loop to remove particles that move outside of simulation area. This
    #is done by calculating grid positions and deleting particles that lie on
    #a grid ID outside the area.
    while i < T_PPC:
        #Open / clear variables for use
        x_i = 0
        y_i = 0

        #Calculate grid locations of particle
        x_i = int(np.floor(p[i].x / dx))
        y_i = int(np.floor(p[i].y / dy))

        #If outside in x domain delete
        if x_i < 0 or x_i > int(nx-1):

            #DEBUG: Print which particle is removed
            print('Out of range: {}'.format(i))

            #Remove particle
            p.pop(int(i))

        else:

            pass

        #If outside in y domain delete
        if y_i < 0 or y_i > int(ny-1):

            #DEBUG: Print which particle is removed
            print('Out of range: {}'.format(i))

            #Remove particle
            p.pop(int(i))

        else:

            pass

        T_PPC = len(p)

        i += 1

    for i in range(0,T_PPC):

        x_i = int(np.floor(p[i].x / dx))
        y_i = int(np.floor(p[i].y / dy))

        #Calculate contribution of each particle's charge to surrounding grid
        #locations in order to calculate current density on grid
        f_ij = (((dx - hx) * (dy - hy)) / (dx * dy))
        f_i1j = ((hx * (dy - hy)) / (dx * dy))
        f_i1j1 = ((hx * hy) / (dx * dy))
        f_ij1 = (((dx - hx)*hy) / (dx*dy))

        J_x[x_i][y_i] += p[i].vx*((p[i].q/(dx*dy))*f_ij)
        J_x[x_i + 1][y_i] += p[i].vx*((p[i].q/(dx*dy))*f_i1j)
        J_x[x_i + 1][y_i + 1] += p[i].vx*((p[i].q/(dx*dy))*f_i1j1)
        J_x[x_i][y_i + 1] += p[i].vx*((p[i].q/(dx*dy))*f_ij1)

        J_y[x_i][y_i] += p[i].vy*((p[i].q/(dx*dy))*f_ij)
        J_y[x_i + 1][y_i] += p[i].vy*((p[i].q/(dx*dy))*f_i1j)
        J_y[x_i + 1][y_i + 1] += p[i].vy*((p[i].q/(dx*dy))*f_i1j1)
        J_y[x_i][y_i + 1] += p[i].vy*((p[i].q/(dx*dy))*f_ij1)

        g.append(Grid_Info(x_i,y_i,f_ij,f_i1j,f_i1j1,f_ij1))

    #Step 1: Advance Electric field ###########################################

    E_x_n_phalf = E_x_n_mhalf + dt*(con.c*((B_z_n - np.roll(B_z_n,1,axis=1))/dy)
     - 4*con.pi*(J_x))
    E_y_n_phalf = E_y_n_mhalf - dt*(con.c*((B_z_n - np.roll(B_z_n,1,axis=0))/dx)
     + 4*con.pi*(J_y))

    #Step 2: Advance Magnetic field ############################################

    B_z_n_pone = B_z_n - con.c*dt*(((E_y_n_phalf-np.roll(E_y_n_phalf,1,axis=0))/dx)
     -((E_x_n_phalf-np.roll(E_x_n_phalf,1,axis=1))/dy))

    #Step 3: Update particle velocities ########################################

    B_norm = B_z_n_pone / con.c

    for i in range(0,T_PPC):

        B = (B_norm[g[i].x_i][g[i].y_i]*g[i].f_ij)+(B_norm[int(g[i].x_i+1)][g[i].y_i]*g[i].f_i1j)+(B_norm[int(g[i].x_i+1)][int(g[i].y_i+1)]*g[i].f_i1j1)+(B_norm[g[i].x_i][int(g[i].y_i+1)]*g[i].f_ij1)
        E_x = (E_x_n_phalf[g[i].x_i][g[i].y_i]*g[i].f_ij)+(E_x_n_phalf[int(g[i].x_i+1)][g[i].y_i]*g[i].f_i1j)+(E_x_n_phalf[int(g[i].x_i+1)][int(g[i].y_i+1)]*g[i].f_i1j1)+(E_x_n_phalf[g[i].x_i][int(g[i].y_i+1)]*g[i].f_ij1)
        E_y = (E_y_n_phalf[g[i].x_i][g[i].y_i]*g[i].f_ij)+(E_y_n_phalf[int(g[i].x_i+1)][g[i].y_i]*g[i].f_i1j)+(E_y_n_phalf[int(g[i].x_i+1)][int(g[i].y_i+1)]*g[i].f_i1j1)+(E_y_n_phalf[g[i].x_i][int(g[i].y_i+1)]*g[i].f_ij1)

        h = (p[i].q*dt)/p[i].m
        K = 1 / (1+((h**2)*(B**2))/4)

        p[i].vx = K*((1-((h**2)*(B**2))/4)*p[i].vx + h*((p[i].vy*B)+E_x) + ((h**2)/2)*(E_x*B))

        p[i].vy = K*((1-((h**2)*(B**2))/4)*p[i].vy + h*(E_y-(p[i].vx*B)) - ((h**2)/2)*(E_y*B))



    #Step 4: Update particle positions #########################################
    for i in range(0,T_PPC):

        p[i].x = float(p[i].x + dt*p[i].vx)
        p[i].y = float(p[i].y + dt*p[i].vy)

    #Reassign variables ########################################################
    E_x_n_mhalf = E_x_n_phalf
    E_y_n_mhalf = E_y_n_phalf
    B_z_n = B_z_n_pone

    #Plotting ##################################################################
    p_x = np.zeros(len(p))
    p_y = np.zeros(len(p))
    p_vx = np.zeros(len(p))
    p_vy = np.zeros(len(p))

    for i in range(0,len(p)):
    	p_x[i] = p[i].x
    	p_y[i] = p[i].y
    	p_vx[i] = p[i].vx
    	p_vy[i] = p[i].vy

    #Debug

    print('vx: min={} max={} avg={}\nvy: min={} max={} avg={}'.format(np.min(p_vx),np.max(p_vx),np.average(p_vx),np.min(p_vy),np.max(p_vy),np.average(p_vy)))

    E_mag_n_mhalf = np.sqrt((E_x_n_mhalf**2)+(E_y_n_mhalf**2))
    E_mag_n_phalf = np.sqrt((E_x_n_phalf**2)+(E_y_n_phalf**2))
    J_mag = np.sqrt((J_x**2)+(J_y**2))

    plt.subplot(2,2,1)
    plt.scatter(p_x,p_y)
    plt.xlim(x_min,x_max)
    plt.ylim(y_min,y_max)
    plt.subplot(2,2,2)
    plt.contourf(X,Y,J_mag)
    plt.colorbar()
    plt.subplot(2,2,3)
    plt.contourf(X,Y,E_mag_n_phalf)
    plt.colorbar()
    plt.subplot(2,2,4)
    plt.contourf(X,Y,B_z_n_pone)
    plt.colorbar()
    plt.draw()
    plt.pause(0.0001)
