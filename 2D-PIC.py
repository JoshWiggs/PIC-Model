import numpy as np
import random as ran
import scipy.constants as con
import matplotlib.pyplot as plt

#Define input parameters
x_min = 0
x_max = 1
y_min = 0
y_max = 1
t_min = 0
t_max = 0.5
v_min = 0
v_max = 10
q_c = 1 #currently unused

nx = 20 #Number of steps taken from y_min to y_max
ny = 20 #Number of steps taken from x_min to x_max
nt = 2 #Number of time steps
n_smoothing = 100 #For differencing methods that require smoothing
PPC = 1 #Number of particles per cell

q_grid= np.zeros(((nx + 1),(ny + 1))) #Charge density
phi = np.zeros(((nx + 1),(ny + 1))) #Electric potential
E_x = np.zeros(((nx + 1),(ny + 1))) #x component of the Electric Field
E_y = np.zeros(((nx + 1),(ny + 1))) #y component of the Electric Field
E_mag = np.zeros((nx,ny)) #magnitude of the Electric field

dx = (abs(x_min) + x_max) / nx #Size of step in x domain
dy = (abs(y_min) + y_max) / ny #Size of step in y domain
dt = (abs(t_min) + t_max) / nt #Size of step on the temporal domain
T_c = nx * ny #Total number of cells
T_PPC = T_c * PPC #Total number of particles in simulation

#Populate dimesnional vectors
x = np.linspace(x_min,x_max,num=nx)
y = np.linspace(y_min,y_max,num=ny)
t = np.linspace(t_min,t_max,num=nt)

#Create list of particles properties
Par_x = np.zeros(T_PPC)
Par_y = np.zeros(T_PPC)
Par_vx = np.zeros(T_PPC)
Par_vy = np.zeros(T_PPC)

for i in range(0,T_PPC):
    Par_x[i] = ran.uniform(x_min,x_max) #Particle x coordinate
    Par_y[i] = ran.uniform(y_min,y_max) #Particle y coordinate
    Par_vx[i] = ran.uniform(v_min,v_max) #Particle velocity on the x direction
    Par_vy[i] = ran.uniform(v_min,v_max) #Particle velocity on the y direction

#Join particle position lists into one array
Par = list(zip(Par_x,Par_y,Par_vx,Par_vy))

for t in range(0,nt):

    q_grid[:][:] = 0
    print(t)

    #Step 1: Compute charge density on grid using bilinear interpolation interpretation
    #(first-order weighting) from particle locations
    g_info = np.zeros((T_PPC,4))
    for i in range(0,T_PPC):

        x_i = 0
        y_i = 0
        hx = 0
        hy = 0

        x_i = np.floor(Par[i][0] / dx)
        y_i = np.floor(Par[i][1] / dy)

        x_i = int(x_i)
        y_i = int(y_i)

        hx = Par[i][0] - x[x_i]
        hy = Par[i][1] - y[y_i]

        #Store grid information relating to each particle for use later
        g_info[i][0] = x_i
        g_info[i][1] = y_i
        g_info[i][2] = hx
        g_info[i][3] = hy

        #TODO: add charge as a particle parameter
        q_grid[x_i][y_i] += ((dx - hx) * (dy - hy) / (dx * dy))
        q_grid[x_i + 1][y_i] += ((hx * (dy - hy)) / (dx * dy))
        q_grid[x_i + 1][y_i + 1] += ((hx * hy) / (dx * dy))
        q_grid[x_i][y_i + 1] += (((dx - hx)*hy) / (dx*dy))

    #Step 2: Compute Electric Potential
    for n in range(0,n_smoothing):
        for i in range(1,nx):
            for j in range(1,ny):

                phi[i][j] = (q_grid[i][j]*((dx**2)*(dy**2)) + ((dy**2)*(phi[i-1][j]
                +phi[i+1][j])) + ((dx**2)*(phi[i][j-1]+phi[i][j+1])) / (2*((dy**2)
                +(dy**2))))

    #Step 3: Compute Electric Field
    for i in range(1,nx):
        for j in range(1,ny):

            E_x[i][j] = ((phi[i + 1][j]) - phi[i - 1][j]) / (2 * dx)

            E_y[i][j] = ((phi[i][j + 1]) - phi[i][j - 1]) / (2 * dy)

            E_mag[i][j] = np.sqrt((E_x[i][j]**2) + (E_y[i][j]**2))

    #Step 4: Move particles
    for i in range(0,T_PPC):

        Par[i] = list(Par[i])

        x_i = int(g_info[i][0])
        y_i = int(g_info[i][1])
        hx = float(g_info[i][2])
        hy = float(g_info[i][3])

        u_x = (E_x[x_i][y_i]*(dx*dy)/((dx-hx)*(dy-hy))) + (E_x[x_i+1][y_i]*((dx*dy)/(hx*(dy-hy)))) + (E_x[x_i+1][y_i+1])*((dx*dy)/(hx*hy)) + (E_x[x_i][y_i+1]*((dx*dy)/(dx-hx)*hy))

        u_y = (E_y[x_i][y_i]*(dx*dy)/((dx-hx)-(dy-hy))) + (E_y[x_i+1][y_i]*((dx*dy)/(hx*(dy-hy)))) + (E_y[x_i+1][y_i+1])*((dx*dy)/(hx*hy)) + (E_y[x_i][y_i+1]*((dx*dy)/(dx-hx)*hy))

        Par[i][2] = u_x
        Par[i][3] = u_y

        Par[i][0] += dt*u_x
        Par[i][1] += dt*u_y

###### PLOTTING #######

new_par_x  = np.zeros(T_PPC)
new_par_y  = np.zeros(T_PPC)

for i in range(0,T_PPC):

    new_par_x[i] = Par[i][0]
    new_par_y[i] = Par[i][1]

X,Y = np.meshgrid(x,y)

plt.subplot(2,1,1)
plt.scatter(Par_x,Par_y)
plt.subplot(2,1,2)
plt.scatter(new_par_x,new_par_y)
plt.show()
