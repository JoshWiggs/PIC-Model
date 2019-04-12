import numpy as np
import random as ran
import scipy.constants as con
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.interpolate as interp
from scipy.spatial import distance
from collections import Counter

#Additional constants
R_j = 7e7 #Jupiter radius [m]

#Options
class options:

    def __init__(self):
            self.DebyeTest = bool(False)
            self.RegenerateParticles = bool(False)
            self.ClosedBoundaries = bool(False) #Not working
            self.JupiterCentre = bool(False)

            self.Electrostatic = bool(True)
            self.Electromagnetic = bool(False)

    class grids:

        def __init__(self):
            self.CartesianGrid = bool(True)
            self.CylindricalGrid = bool(False) #still under development

    class visual:

        def __init__(self):
            self.write_movie = bool(False)

    class boundaries:

        def __init__(self):
            self.CylHardInner = bool(False) #for cylindrical


opt = options()
vis = options.visual()
grid = options.grids()
bound = options.boundaries()

def closest_point(node, nodes):
    nodes = np.asarray(nodes)
    dist_2 = (nodes - node)**2
    return np.argmin(dist_2)

#Ensure one coordinate system has been selected
if grid.CartesianGrid is False and grid.CylindricalGrid is False:

    raise Exception('Please select a coordinate system for the grid.')

elif grid.CartesianGrid is True and grid.CylindricalGrid is True:

    raise Exception('Please select only one coordinate system.')

##### CURRENTLY UNUSED ##########
class maxwell_velocities(stats.rv_continuous):
    def _pdf(self,v):
        return (((con.m_e)/(2*con.pi*con.k*1e6))**(3/2))*(4*con.pi*(v**2))*np.exp(-((con.m_e*(v)**2)/(2*con.k*1e6)))

velocity_pdf = maxwell_velocities(name='velocity_pdf', a=0.0, b=con.c)
#################################

def thermal_velocity(T,Par_m):

    #Open variables
    varience,u1,u2,v = float(0)

    #Calculate varience in eV
    varience = T*1.602e-19/Par_m

    #Generate random uniform numbers between 0 & 1
    u1 = ran.uniform(0.0,1.0)
    u2 = ran.uniform(0.0,1.0)

    #use Box-Muller transform to generate randomly-distributed normal
    v = np.sqrt(-2.0*np.log(u1))*np.cos(2.0*np.pi*u2)*np.sqrt(variance)

    v_x = v*u1
    v_y = v*u2

if vis.write_movie is True:

    #Import modules for gif writing
    import matplotlib
    matplotlib.use('Agg')
    import imageio
    import matplotlib.gridspec as gridspec

#Define input parameters
if grid.CartesianGrid is True:

    x_min = 0 #Minimum position in x domain
    x_max = 5*R_j #Maximum position in x domain
    y_min = 0 #Minimum position in y domain
    y_max = 5*R_j #Maximum position in y domain

    nx = 100 #Number of steps taken from x_min to x_max
    ny = 100#Number of steps taken from y_min to y_max

    dx = (x_max - x_min) / nx #Size of step in x domain
    dy = (y_max - y_min) / ny #Size of step in y domain
    T_c = nx * ny #Total number of cells

    v_min = -R_j #Minimum initial particle speed
    v_max = R_j #Maximum initial particle speed

    #Populate dimesnional vectors
    x = np.linspace(x_min,x_max,num=(nx+1))
    y = np.linspace(y_min,y_max,num=(ny+1))

    #Meshgrid for plotting
    X,Y = np.meshgrid((x/R_j),(y/R_j))

    phi = np.zeros(((nx + 1),(ny + 1))) #Electric potential
    if opt.Electrostatic is True:
        E_x = np.zeros(((nx + 1),(ny + 1))) #x component of the Electric Field
        E_y = np.zeros(((nx + 1),(ny + 1))) #y component of the Electric Field
        E_mag = np.zeros((nx,ny)) #Magnitude of the Electric field

    if opt.Electromagnetic is True:
        E_x_int = np.zeros(((nx + 1),(ny + 1)))
        E_y_int = np.zeros(((nx + 1),(ny + 1)))
        E_x_half = np.zeros(((nx + 1),(ny + 1)))
        E_y_half = np.zeros(((nx + 1),(ny + 1)))
        E_x = np.zeros((2*(nx + 1),2*(ny + 1)))
        E_y = np.zeros((2*(nx + 1),2*(ny + 1)))
        B_z = np.zeros(((nx + 1),(ny + 1)))
        B_z_n = np.zeros(((nx + 1),(ny + 1)))
        B_z_mhalf = np.zeros(((nx + 1),(ny + 1)))

elif grid.CylindricalGrid is True:

    r_min = 1 #Minimum position in the r domain
    r_max = 2 #Maximum position in the r domain
    theta_min = 0 #Minimum position in the theta domain
    theta_max = 2*con.pi #Maximum position in the theta domain

    nr = 50 #Number of steps taken from r_min to r_max
    ntheta = 50 #Number of steps taken from theta_min to theta_max

    dr = (r_max - r_min) / nr #Size of step in r domian
    dtheta = (theta_max - theta_min) / ntheta #Size of step in theta domain
    T_c = nr * ntheta #Total number of cells

    vr_min = -(2*dr) #Minimum intial particle speed in the radial domain
    vr_max = (2*dr) #Maximum intial particle speed in the radial domain
    vtheta_min = -(2*dtheta) #Minimum intial particle speed in the azimuthal domain
    vtheta_max = (2*dtheta) #Maximum intial particle speed in the azimuthal domain

    phi = np.zeros(((nr + 1),(ntheta + 1))) #Electric potential
    E_x = np.zeros(((nr + 1),(ntheta + 1))) #x component of the Electric Field
    E_y = np.zeros(((nr + 1),(ntheta + 1))) #y component of the Electric Field
    E_mag = np.zeros((nr,ntheta)) #magnitude of the Electric field

    #Populate dimesnional vectors
    r = np.linspace(r_min,r_max,num=(nr+1))
    theta = np.linspace(theta_min,theta_max,num=(ntheta+1))

    #Meshgrid for plotting
    R,Theta = np.meshgrid(r,theta)

#Time Variables
t_min = 0 #Start time
t_max = 2.0 #End time
nt = 200 #Number of time steps
t = np.linspace(t_min,t_max,num=nt) #Populate dimesnional vectors
dt = (abs(t_min) + t_max) / nt #Size of step on the temporal domain

#Particle charge & mass lists
q_c = [(10*con.e),-(100*con.e)] #particle charge
m = [(10*con.m_p),(100*con.m_e)] #particle mass
#q_c = [-con.e] #particle charge
#m = [con.m_e] #particle mass
epsilon =  1 #con.epsilon_0 / R_j

n_smoothing = 50 #For differencing methods that require smoothing
PPC = .1 #Number of particles per cell

T_PPC = int(T_c * PPC) #Total number of particles in simulation

if grid.CartesianGrid is True:
    #Create list of particles properties
    Par_x = np.zeros(T_PPC)
    Par_y = np.zeros(T_PPC)
    Par_vx = np.zeros(T_PPC)
    Par_vy = np.zeros(T_PPC)
elif grid.CylindricalGrid is True:
    #Create list of particles properties
    Par_r = np.zeros(T_PPC)
    Par_theta = np.zeros(T_PPC)
    Par_vr = np.zeros(T_PPC)
    Par_vtheta = np.zeros(T_PPC)

Par_q = np.zeros(T_PPC)
Par_m = np.zeros(T_PPC)
Par_species_count = np.zeros(len(q_c))

if vis.write_movie is True:

    #Make a figure, get the axes
    fig = plt.Figure()
    ax = fig.gca() #projection='polar'
    gs = gridspec.GridSpec(2, 1, figure=fig)

    # Set the figure canvas to the Agg backend and get the width/height of the canvas (in pixels)
    fig.canvas = matplotlib.backends.backend_agg.FigureCanvasAgg(fig)
    w,h = fig.canvas.get_renderer().get_canvas_width_height()

    # Setup the imageio writer.
    writer = imageio.get_writer('polar_test.gif', fps=25)

if grid.CartesianGrid is True:

    fig, ax = plt.subplots(2, 1, constrained_layout=False)

if grid.CylindricalGrid is True:

    fig, ax = plt.subplots(1, 2, subplot_kw=dict(projection='polar'), constrained_layout=True)

print('Initialising {} particles...'.format(T_PPC))

for i in range(0,T_PPC):

    #TODO: Generalise
    j = ran.randint(0,int(len(q_c)-1)) #Random integer for mass and charge selection
    s = ran.choice((-1,1))

    if grid.CartesianGrid is True:

        Par_x[i] = ran.uniform(x_min,x_max) #Particle x coordinate
        Par_y[i] = ran.uniform(y_min,y_max) #Particle y coordinate
        Par_vx[i] = ran.uniform(v_min,v_max) #Particle velocity on the x direction
        Par_vy[i] = ran.uniform(v_min,v_max) #Particle velocity on the y direction
        #Par_vx[i] = s*(velocity_pdf.rvs() * 1)  #Particle velocity on the x direction
        #Par_vy[i] = s*(velocity_pdf.rvs() * 1) #Particle velocity on the y direction

    elif grid.CylindricalGrid is True:

        Par_r[i] = ran.uniform(r_min,r_max)
        Par_theta[i] = ran.uniform(theta_min,theta_max)
        Par_vr[i] = ran.uniform(vr_min,vr_max)
        Par_vtheta[i] = ran.uniform(vtheta_min,vtheta_max)

    Par_q[i] = q_c[j] #Particle charge
    Par_m[i] = m[j] #Particle mass

#Join particle position lists into one array
#TODO: Rework this as data types produced seem to make calcultions difficult
if grid.CartesianGrid is True:
    Par = list(zip(Par_x,Par_y,Par_vx,Par_vy,Par_q,Par_m))
elif grid.CylindricalGrid is True:
    Par = list(zip(Par_r,Par_theta,Par_vr,Par_vtheta,Par_q,Par_m))

#If the Debye test is active then generate the Debye test particle and
#append it as the final item in the particle list
if opt.DebyeTest is True:

    if grid.CartesianGrid is True:

        print('Creating Debye particle...')

        D_Par_x = ((x_min + x_max) / 2) + (dx / 2.)
        D_Par_y = ((y_min + y_max) / 2) + (dy / 4.)
        D_Par_vx = 0
        D_Par_vy = 0
        D_Par_q = -(con.e*100)
        D_Par_m = 1.

        D_Par = [D_Par_x,D_Par_y,D_Par_vx,D_Par_vy,D_Par_q,D_Par_m]

        Par.append(D_Par)

    else:

        raise Exception('This test is only for a cartesian coordinate system')

if opt.JupiterCentre is True:

    if grid.CartesianGrid is True:

        print('Creating Jupiter parameters...')

        Jupiter_gridsize_x = int(1/(dx/R_j))
        Jupiter_gridsize_y = int(1/(dy/R_j))
        Jupiter_r = np.sqrt(((((int(1/(dx/R_j)))/2)*(dx/R_j))**2) + ((((int(1/(dy/R_j)))/2)*(dy/R_j))**2))

        grid_min_x = int((nx/2) - (Jupiter_gridsize_x/2))
        grid_max_x = int((nx/2) + (Jupiter_gridsize_x/2))
        grid_min_y = int((ny/2) - (Jupiter_gridsize_y/2))
        grid_max_y = int((ny/2) + (Jupiter_gridsize_y/2))
        for i in range(0,T_PPC):
            Par[i] = list(Par[i])
            if x[grid_min_x]<=float(Par[i][0])<=x[grid_max_x]:
                Par[i][0] = float(Par[i][0]) + R_j
            if y[grid_min_y]<=float(Par[i][1])<=y[grid_max_y]:
                Par[i][1] = float(Par[i][1]) + R_j
    else:

        raise Exception('Jupiter Centre option is only for a cartesian coordinate system')


if grid.CylindricalGrid is True:
    if opt.RegenerateParticles is True:
        print("That's no moon... It's a macro particle" )
        Moon_par_r = (r_min + ((r_max - r_min) / 2.))
        Moon_par_theta = ran.uniform(theta_min,theta_max)
        Moon_par_vr = 0.
        Moon_par_vtheta = (dtheta / 2.)
        Moon_par_q = 0.
        Moon_par_m = 10000

        Moon_Par = [Moon_par_r,Moon_par_theta,Moon_par_vr,Moon_par_vtheta,Moon_par_q,Moon_par_m]

        Par.append(Moon_Par)

        T_PPC = len(Par)

print('Initialising particle simulation...')
if opt.Electromagnetic is True:

    #Calculate initial conditions
    B_z_mhalf = ((X**2)+(Y**2))**(3/2)
    q_grid = np.zeros(((nx + 1),(ny + 1))) #Charge density
    for i in range(0,len(T_PPC)):
        x_i = 0 #Particle x grid position
        y_i = 0 #Particle y grid position
        hx = 0 #Partcle distance from its grid postion in the x domain
        hy = 0 #Partcle distance from its grid postion in the y domain

        x_i = int(np.floor(Par[i][0] / dx))
        y_i = int(np.floor(Par[i][1] / dy))

        hx = Par[i][0] - x[x_i]
        hy = Par[i][1] - y[y_i]

        q_grid[x_i][y_i] += (Par[i][4]/(dx*dy))*(((dx - hx) * (dy - hy)) / (dx * dy))
        q_grid[x_i + 1][y_i] += (Par[i][4]/(dx*dy))*((hx * (dy - hy)) / (dx * dy))
        q_grid[x_i + 1][y_i + 1] += (Par[i][4]/(dx*dy))*((hx * hy) / (dx * dy))
        q_grid[x_i][y_i + 1] += (Par[i][4]/(dx*dy))*(((dx - hx)*hy) / (dx*dy))

    phi = (((q_grid/epsilon)*((dx**2)*(dy**2)) + ((dy**2)*(np.roll(phi,1,axis=0) + np.roll(phi,-1,axis=0))) + ((dx**2)*(np.roll(phi,1,axis=1) + np.roll(phi,-1,axis=1)))) / (2*((dx**2) + (dy**2))))

    for n in range(0, n_smoothing):
        E_x_int = ((np.roll(phi,-1,axis=0)) - np.roll(phi,1,axis=0)) / (2 * dx)
        E_y_int = ((np.roll(phi,-1,axis=1)) - np.roll(phi,1,axis=1)) / (2 * dy)
        E_mag = np.sqrt((E_x**2) + (E_y**2))

    E_x_half = (E_x+np.roll(E_x,1,axis=0))/2
    E_y_half = (E_y+np.roll(E_y,1,axis=1))/2


    for i in range(0,len(E_x)):
        if i % 2 == 0:
            E_x[i] = E_x_int[(i/2)]
            E_y[i] = E_y_int[(i/2)]
        else:
            E_x[i] = E_x_half[int((i/2)-.5)]
            E_y[i] = E_y_half[int((i/2)-.5)]

#Numerical loop responsible for interating the particles through time
for t in range(0,nt): #nt

    #Activate condition only if writing gif
    if vis.write_movie is True:

        #Clear the axes so we can draw the new frame
        ax[0].cla()
        ax[1].cla()

    #Print index 't' every 10 time steps to show progress
    if t % 10 == 0 and t > 0:

        print(t)
        print('Ions: {}, Electrons: {}'.format(Par_species_count[q_c[0]],Par_species_count[q_c[1]]))

    #Perpare index 'i' for while loop responsible for removing particles which
    #move outside of the simulation region
    i = 0

    #While loop to remove particles that move outside of simulation area. This
    #is done by calculating grid positions and deleting particles that lie on
    #a grid ID outside the area.
    while i < T_PPC:

        if grid.CartesianGrid is True:
            #Open / clear variables for use
            x_i = 0
            y_i = 0

            #Calculate grid locations of particle
            x_i = np.floor(Par[i][0] / dx)
            y_i = np.floor(Par[i][1] / dy)

            #TODO: Fixed closed boundary condition
            if opt.ClosedBoundaries is True:

                Par[i] = list(Par[i])

                u_x = 0.
                u_y = 0.

                if x_i == 0:
                    Par[i][2] = u_x
                    if u_x < 0:
                        Par[i][2] = (-1 * u_x)

                if x_i == (nx-1):
                    Par[i][2] = u_x
                    if u_x > 0:
                        Par[i][2] = (-1 * u_x)

                if y_i == 0:
                    Par[i][3] = u_y
                    if u_y < 0:
                        Par[i][3] = (-1 * u_y)

                if y_i == (ny-1):
                    Par[i][3] = u_y
                    if u_y > 0:
                        Par[i][3] = (-1 * u_y)

            else:

                #If outside in x domain delete
                if x_i < 0 or x_i > int(nx-1):

                    #DEBUG: Print which particle is removed
                    print('Out of range: {}'.format(i))

                    #Remove particle
                    Par.pop(int(i))

                else:

                    pass

                #If outside in y domain delete
                if y_i < 0 or y_i > int(ny-1):

                    #DEBUG: Print which particle is removed
                    print('Out of range: {}'.format(i))

                    #Remove particle
                    Par.pop(int(i))

                else:

                    pass

        if grid.CylindricalGrid is True:
            #Open / clear variables for use
            r_i = 0
            theta_k = 0
            r_p = Par[i][0]

            #Calculate grid locations of particle
            r_i = np.floor((Par[i][0] - r_min) / dr)
            theta_k = np.floor(Par[i][1] / dtheta)

            if bound.CylHardInner is True:
                #Closed inner boundary
                if r_p <= (r_min + 0.05):
                    Par[i] = list(Par[i])
                    u_rp = 0.
                    u_rp = float(Par[i][2])
                    if u_rp < 0:
                        Par[i][2] = -1*(u_rp)

            #If outside in radial domain delete
            if r_i < 0 or r_i > int(nr-1):

                #DEBUG: Print which particle is removed
                print('Out of range: {}'.format(i))

                #Remove particle
                Par.pop(int(i))

            else:

                pass

            #Ensure continuous motion in azimuthal domain
            if theta_k < 0:

                par_theta = float(Par[i][1]) + (2*con.pi)
                Par[i][1] = par_theta

            elif theta_k > int(ntheta-1):

                par_theta = float(Par[i][1]) - (2*con.pi)
                Par[i][1] = par_theta

            else:

                pass

        #Determine number of each particle species contained within the
        #simulation area

        List_Par_q = [i[4] for i in Par]

        Par_species_count = Counter(List_Par_q)

        #Recalculate number of particle in simulation area in case current
        #particle in loop has been deleted
        T_PPC = len(Par)

        #Progress to next particle
        i += 1

    #If option to regenerate particles that leave the simulation domain is
    #active then include this step in iteration loop
    if opt.RegenerateParticles is True:

        #Calculate number of particles that require regenerating
        Regenerate_num = int((T_c * PPC) - len(Par))

        #If particles require regenerating then active loop to generate insertion
        #parameters for as many new particles as required
        if Regenerate_num > 0:
            for re in range(0,Regenerate_num):

                #This works in the same way as the loop in the pre-process
                j = ran.randint(0,int(len(q_c)-1))

                if grid.CartesianGrid is True:
                    Regen_par_x = ran.uniform(x_min,x_max)
                    Regen_par_y = ran.uniform(y_min,y_max)
                    Regen_par_vx = ran.uniform(v_min,v_max)
                    Regen_par_vy = ran.uniform(v_min,v_max)
                    Regen_par_q = q_c[j]
                    Regen_par_m = m[j]

                    Regen_par = [Regen_par_x,Regen_par_y,Regen_par_vx,Regen_par_vy,Regen_par_q,Regen_par_m]

                elif grid.CylindricalGrid is True:

                    Moon_pos_r = Par[int(T_PPC - 1)][0]
                    Moon_pos_theta = Par[int(T_PPC - 1)][1]

                    Regen_par_r = Moon_pos_r
                    Regen_par_theta = ran.uniform((Moon_pos_theta - (con.pi / 4)),(Moon_pos_theta + (con.pi / 4)))
                    Regen_par_vr = ran.uniform(vr_min,vr_max)
                    Regen_par_vtheta = ran.uniform(vtheta_min,vtheta_max)
                    Regen_par_q = q_c[j]
                    Regen_par_m = m[j]

                    Regen_par = [Regen_par_r,Regen_par_theta,Regen_par_vr,Regen_par_vtheta,Regen_par_q,Regen_par_m]

                #Insert regenerated particles into particle list preserving the
                #position of the current last element incase Debye test is
                #active
                Par.insert(int(len(Par)-2),Regen_par)

    #DEBUG: Inform that all particles have been checked to ensure they are in
    #the simulation area
    #print('Position check completed')

    #Recalculate the number of particles in the simulation area. THIS IS A
    #SAFETY CHECK AS IT COSTS VERY LITTLE COMPUTATIONALLY
    T_PPC = len(Par)

    #Clear grid data & charge density array ready for re-use
    if grid.CartesianGrid is True:
        g_info = np.zeros((T_PPC,6)) #Grid information
    elif grid.CylindricalGrid is True:
        g_info = np.zeros((T_PPC,6))

    if grid.CartesianGrid is True:
        q_grid = np.zeros(((nx + 1),(ny + 1))) #Charge density
        J_x_half = np.zeros(((nx + 1),(ny + 1)))
        J_y_half = np.zeros(((nx + 1),(ny + 1)))
        J_x_int = np.zeros(((nx + 1),(ny + 1)))
        J_x_int = np.zeros(((nx + 1),(ny + 1)))
        J_x = np.zeros((2*(nx + 1),2*(ny + 1))) #Current density in the x-domain
        J_y = np.zeros((2*(nx + 1),2*(ny + 1))) #Current density in the y-domain
    elif grid.CylindricalGrid is True:
        q_grid = np.zeros(((nr + 1),(ntheta + 1))) #Charge density
        Par_x = np.zeros(T_PPC)
        Par_y = np.zeros(T_PPC)

    if opt.Electrostatic is True:

        for i in range(0,T_PPC):

            #Step 1: Compute charge density on grid using bilinear interpolation interpretation
            #(first-order weighting) from particle locations
            if grid.CartesianGrid is True:
                x_i = 0 #Particle x grid position
                y_i = 0 #Particle y grid position
                hx = 0 #Partcle distance from its grid postion in the x domain
                hy = 0 #Partcle distance from its grid postion in the y domain

                #Calculate grid locations of particle
                x_i = int(np.floor(Par[i][0] / dx))
                y_i = int(np.floor(Par[i][1] / dy))

                #Calculate distance of particle from grid location
                #DEBUG: If a particle outside the simulation area has not been removed
                #       raise error and return particle information
                try:
                    hx = Par[i][0] - x[x_i]
                    hy = Par[i][1] - y[y_i]
                except:
                    raise Exception('ERROR: An error occurred during timestep {} check particle {}'.format(t,i))

                #Store grid information relating to each particle for use later
                g_info[i][0] = x_i
                g_info[i][1] = y_i
                #g_info[i][2] = hx
                #g_info[i][3] = hy

                #Calculate contribution of each particle's charge to surrounding grid
                #locations in order to calculate charge density on grid
                f_ij = (((dx - hx) * (dy - hy)) / (dx * dy))
                f_i1j = ((hx * (dy - hy)) / (dx * dy))
                f_i1j1 = ((hx * hy) / (dx * dy))
                f_ij1 = (((dx - hx)*hy) / (dx*dy))

                g_info[i][2] = f_ij
                g_info[i][3] = f_i1j
                g_info[i][4] = f_i1j1
                g_info[i][5] = f_ij1

                q_grid[x_i][y_i] += (Par[i][4]/(dx*dy))*f_ij
                q_grid[x_i + 1][y_i] += (Par[i][4]/(dx*dy))*f_i1j
                q_grid[x_i + 1][y_i + 1] += (Par[i][4]/(dx*dy))*f_i1j1
                q_grid[x_i][y_i + 1] += (Par[i][4]/(dx*dy))*f_ij1

                #DEBUG: Check area weighting total 1
                f = round((f_ij+f_i1j+f_i1j1+f_ij1),0)
                if f != 1:
                    raise Exception('Error: Particle {} area weighting {}'.format(i,f))

            elif grid.CylindricalGrid is True:
                r_i = 0 #Particle radial grid position
                theta_k = 0 #Particle azimuthal grid position
                hr = 0 #Partcle distance from its grid postion in the radial domain
                htheta = 0 #Partcle distance from its grid postion in the azimuthal domain

                #Calculate grid positions of particles
                r_i = int(np.floor((Par[i][0] - r_min) / dr))
                theta_k = int(np.floor(Par[i][0] / dtheta))

                #Calculate distance of particle from grid location
                #DEBUG: If a particle outside the simulation area has not been removed
                #       raise error and return particle information
                try:
                    hr = Par[i][0] - r[r_i]
                    if hr == 0:
                        hr = 0.0000001
                    htheta = Par[i][1] - theta[theta_k]

                except:
                    raise Exception('ERROR: An error occurred during timestep {} check particle {}'.format(t,i))

                #Store grid information relating to each particle for use later
                g_info[i][0] = r_i
                g_info[i][1] = theta_k
                #g_info[i][2] = hr
                #g_info[i][3] = htheta

                V_ik = (dtheta/2)*(((r[r_i]+(dr/2))**2) - ((r[r_i]-(dr/2))**2))
                area_ik = (((r[r_i+1]**2) - (r[r_i]**2))*(theta[theta_k+1]-theta[theta_k]))

                f_ik = (((r[r_i+1]**2)-(float(Par[i][0])**2))*(theta[theta_k+1]-float(Par[i][1])))/area_ik
                f_i1k = (((float(Par[i][0])**2)-(r[r_i]**2))*(theta[theta_k+1]-float(Par[i][1])))/area_ik
                f_i1k1 = (((float(Par[i][0])**2)-(r[r_i]**2))*(float(Par[i][1])-theta[theta_k]))/area_ik
                f_ik1 = (((r[r_i+1]**2)-(float(Par[i][0])**2))*(float(Par[i][1])-theta[theta_k]))/area_ik

                g_info[i][2] = f_ik
                g_info[i][3] = f_i1k
                g_info[i][4] = f_i1k1
                g_info[i][5] = f_ik1

                q_grid[r_i][theta_k] += ((Par[i][4])*f_ik) / V_ik
                q_grid[r_i + 1][theta_k] += ((Par[i][4])*f_i1k) / V_ik
                q_grid[r_i + 1][theta_k + 1] += ((Par[i][4])*f_i1k1) / V_ik
                q_grid[r_i][theta_k + 1] += ((Par[i][4])*f_ik1) / V_ik

            #Par_r_list = [i[0] for i in Par]
            #Par_theta_list = [i[1] for i in Par]
            #points = list(zip(Par_r_list,Par_theta_list))
            #Par_q_list = [i[4] for i in Par]
            #q_grid = interp.griddata(points,Par_q_list,(R,Theta))

        #Step 2: Compute Electric Potential
        for n in range(0,n_smoothing):

            if grid.CartesianGrid is True:

                """
                for i in range(1,nx):
                    for j in range(1,ny):

                        phi[i][j] = (q_grid[i][j]*((dx**2)*(dy**2)) + ((dy**2)*
                        (phi[i-1][j] +phi[i+1][j])) + ((dx**2)*(phi[i][j-1]+phi[i][j+1]))
                        / (2*((dy**2) +(dy**2))))
                """

                #TODO: need to check this eqn
                phi = (((q_grid/epsilon)*((dx**2)*(dy**2)) + ((dy**2)*(np.roll(phi,1,axis=0)
                + np.roll(phi,-1,axis=0))) + ((dx**2)*(np.roll(phi,1,axis=1)
                + np.roll(phi,-1,axis=1)))) / (2*((dx**2) + (dy**2))))

            elif grid.CylindricalGrid is True:

                phi = ((q_grid*((r+(dr/2)**2)-(r-(dr/2)**2))) + (((2*(r+(dr/2)))/(np.roll(r,-1,axis=0)-r))*np.roll(phi,-1,axis=0)) + (((2*(r-(dr/2)))/(r-np.roll(r,1,axis=0)))*np.roll(phi,1,axis=0)) + (((2*((r+(dr/2))-(r-(dr/2))))/(((dtheta)**2))*r)*(np.roll(phi,-1,axis=1)+np.roll(phi,1,axis=1))))/ (((2*(r+(dr/2)))/(np.roll(r,-1,axis=0)-r)) + ((2*(r-(dr/2)))/(r-np.roll(r,1,axis=0))) + ((4*((r+(dr/2))-(r-(dr/2))))/(((dtheta)**2))*r))
                #phi = ((((q_grid)*((dr**2)*(dtheta**2))) + ((dtheta**2)*(np.roll(phi,1,axis=0) + np.roll(phi,-1,axis=0))) + (((dr**2)/(r**2))*(np.roll(phi,1,axis=1) + np.roll(phi,-1,axis=1))) +((((dtheta**2)*dr)/(2*r))*(np.roll(phi,-1,axis=0)-np.roll(phi,1,axis=0)))) / (2*((dtheta**2)+((dr**2)/r))))

                if bound.CylHardInner is True:

                    phi[0][:] = 0
                    #phi[nr][:] = 0

        #Step 3: Compute Electric Field
        if grid.CartesianGrid is True:

            """
            for i in range(1,nx):
                for j in range(1,ny):

                E_x[i][j] = ((phi[i + 1][j]) - phi[i - 1][j]) / (2 * dx)

                E_y[i][j] = ((phi[i][j + 1]) - phi[i][j - 1]) / (2 * dy)

                E_mag[i][j] = np.sqrt((E_x[i][j]**2) + (E_y[i][j]**2))
            """
            for n in range(0, n_smoothing):

                E_x = ((np.roll(phi,-1,axis=0)) - np.roll(phi,1,axis=0)) / (2 * dx)

                E_y = ((np.roll(phi,-1,axis=1)) - np.roll(phi,1,axis=1)) / (2 * dy)

                E_mag = np.sqrt((E_x**2) + (E_y**2))

        elif grid.CylindricalGrid is True:
            #TODO: Look at half-grid implications in Birdsall & Langdon
            for n in range(0, n_smoothing):

                E_r_half = -(((np.roll(phi,-1,axis=0)) - np.roll(phi,1,axis=0)) / (np.roll(r,-1,axis=0) - r))

                E_theta = -(((np.roll(phi,-1,axis=1)) - np.roll(phi,1,axis=1)) / (r * dtheta))

                E_r = (((r+(dr/2))*np.roll(E_r_half,-1,axis=0)) + ((r-(dr/2))*np.roll(E_r_half,1,axis=0))) / (2*r)

                E_mag = np.sqrt((E_r**2) + (E_theta**2))


        #Step 4: Move particles
        for i in range(0,T_PPC):

            #Safety check for particle information data type
            Par[i] = list(Par[i])

            #Access particle information and the grid information relating to they
            #particle for use in calculating particle motion
            if grid.CartesianGrid is True:
                x_i = int(g_info[i][0])
                y_i = int(g_info[i][1])
                u_x_old = float(Par[i][2])
                u_y_old = float(Par[i][3])
                #hx = float(g_info[i][2])
                #hy = float(g_info[i][3])
                f_ij = float(g_info[i][2])
                f_i1j = float(g_info[i][3])
                f_i1j1 = float(g_info[i][4])
                f_ij1 = float(g_info[i][5])
            elif grid.CylindricalGrid is True:
                r_i = int(g_info[i][0])
                theta_k = int(g_info[i][1])
                u_x_old = float(Par[i][2])
                u_y_old = float(Par[i][3])
                f_ik = float(g_info[i][2])
                f_i1k = float(g_info[i][3])
                f_i1k1 = float(g_info[i][4])
                f_ik1 = float(g_info[i][5])

            #Clear variables for reuse
            E_Par_x = float(0) #Electric field felt by the particle in the x domain
            E_Par_y = float(0) #Electric field felt by the particle in the y domain

            #Calculate the electric field strength at the particle from each of the
            #surrounding grid points using bilinear interpolation interpretation

            if grid.CartesianGrid is True:
                #E_Par_x = (E_x[x_i][y_i]*(dx*dy)/((dx-hx)*(dy-hy))) + (E_x[x_i+1][y_i]*((dx*dy)/(hx*(dy-hy)))) + (E_x[x_i+1][y_i+1])*((dx*dy)/(hx*hy)) + (E_x[x_i][y_i+1]*((dx*dy)/(dx-hx)*hy))
                #E_Par_y = (E_y[x_i][y_i]*(dx*dy)/((dx-hx)-(dy-hy))) + (E_y[x_i+1][y_i]*((dx*dy)/(hx*(dy-hy)))) + (E_y[x_i+1][y_i+1])*((dx*dy)/(hx*hy)) + (E_y[x_i][y_i+1]*((dx*dy)/(dx-hx)*hy))

                E_Par_x = (E_x[x_i][y_i]*f_ij) + (E_x[x_i+1][y_i]*f_i1j) + (E_x[x_i+1][y_i+1]*f_i1j1) + (E_x[x_i][y_i+1]*f_ij1)
                E_Par_y = (E_y[x_i][y_i]*f_ij) + (E_y[x_i+1][y_i]*f_i1j) + (E_y[x_i+1][y_i+1]*f_i1j1) + (E_y[x_i][y_i+1]*f_ij1)

            elif grid.CylindricalGrid is True:
                E_Par_x = (E_r[r_i][theta_k]*f_ik) + (E_r[r_i+1][theta_k]*f_i1k) + (E_r[r_i+1][theta_k+1]*f_i1k1) + (E_r[r_i][theta_k+1]*f_ik1)
                E_Par_y = (E_theta[r_i][theta_k]*f_ik) + (E_theta[r_i+1][theta_k]*f_i1k) + (E_theta[r_i+1][theta_k+1]*f_i1k1) + (E_theta[r_i][theta_k+1]*f_ik1)

            #Using updated electic field calculate updated particle velocities
            u_x = u_x_old - (((dt * Par[i][4]) / Par[i][5]) * E_Par_x)
            u_y = u_y_old - (((dt * Par[i][4]) / Par[i][5]) * E_Par_y)

            ##### If Debye physical test is being conducted then the final particle
            ##### in the list is the stationary Debye test particle. Therefore, the
            ##### velocities are set to zero to prevent motion
            if opt.DebyeTest is True:
                if i == int(T_PPC-1):

                    u_x = 0
                    u_y = 0

            if grid.CylindricalGrid is True:
                if opt.RegenerateParticles is True:
                    if i == int(T_PPC-1):

                        u_x = Moon_par_vr
                        u_y = Moon_par_vtheta

            #Update paricle velocity
            Par[i][2] = u_x
            Par[i][3] = u_y

            #Use particle velocity to calculate new particle position
            Par[i][0] += dt*u_x
            Par[i][1] += dt*u_y

    elif opt.Electromagnetic is True:

        #Step 0: Calculate charge density on all grid points
        # and store bilinear interpolation weights for use later
        if grid.CartesianGrid is True:
            x_i = 0 #Particle x grid position
            y_i = 0 #Particle y grid position
            hx = 0 #Partcle distance from its grid postion in the x domain
            hy = 0 #Partcle distance from its grid postion in the y domain

            #Calculate grid locations of particle
            x_i = int(np.floor(Par[i][0] / dx))
            y_i = int(np.floor(Par[i][1] / dy))

            #Calculate distance of particle from grid location
            #DEBUG: If a particle outside the simulation area has not been removed
            #       raise error and return particle information
            try:
                hx = Par[i][0] - x[x_i]
                hy = Par[i][1] - y[y_i]
            except:
                raise Exception('ERROR: An error occurred during timestep {} check particle {}'.format(t,i))

            #Store grid information relating to each particle for use later
            g_info[i][0] = x_i
            g_info[i][1] = y_i

            #Calculate contribution of each particle's charge to surrounding grid
            #locations in order to calculate charge density on grid
            f_ij = (((dx - hx) * (dy - hy)) / (dx * dy))
            f_i1j = ((hx * (dy - hy)) / (dx * dy))
            f_i1j1 = ((hx * hy) / (dx * dy))
            f_ij1 = (((dx - hx)*hy) / (dx*dy))

            g_info[i][2] = f_ij
            g_info[i][3] = f_i1j
            g_info[i][4] = f_i1j1
            g_info[i][5] = f_ij1

            q_grid[x_i][y_i] += (Par[i][4]/(dx*dy))*f_ij
            q_grid[x_i + 1][y_i] += (Par[i][4]/(dx*dy))*f_i1j
            q_grid[x_i + 1][y_i + 1] += (Par[i][4]/(dx*dy))*f_i1j1
            q_grid[x_i][y_i + 1] += (Par[i][4]/(dx*dy))*f_ij1

            #DEBUG: Check area weighting total 1
            f = round((f_ij+f_i1j+f_i1j1+f_ij1),0)
            if f != 1:
                raise Exception('Error: Particle {} area weighting {}'.format(i,f))

        #Step 1: Calculate B_z^n+1/2 using induction eqn
        for i in range(0,(2*(nx+1)),2):
            for j in range(0,(2*(ny+1)),2):
                B_z = B_z_mhalf - dt * (((E_y[int(i+2)][int(j+1)]-E_y[i][int(j+1)])/dx)-((E_x[int(i+1)][int(j+2)]-E_x[int(i+1)][j])/dy))

        #Step1.5: Calculate B^n

        B_z_n = (B_z +B_z_mhalf)/2

        #Step 2: Calculate u^n+1/2 using eqn for Lorentz force
        for i in range(0,T_PPC):
            if grid.CartesianGrid is True:
                x_i = int(g_info[i][0])
                y_i = int(g_info[i][1])
                u_x_mhalf = float(Par[i][2])
                u_y_mhalf = float(Par[i][3])
                f_ij = float(g_info[i][2])
                f_i1j = float(g_info[i][3])
                f_i1j1 = float(g_info[i][4])
                f_ij1 = float(g_info[i][5])

            #Clear variables for reuse
            E_Par_x = float(0) #Electric field felt by the particle in the x domain
            E_Par_y = float(0) #Electric field felt by the particle in the y domain
            u_x_phalf = float(0)
            u_y_phalf = float(0)
            k = float(0)
            h = float(0)

            E_Par_x = (E_x_int[x_i][y_i]*f_ij) + (E_x_int[x_i+1][y_i]*f_i1j) + (E_x_int[x_i+1][y_i+1]*f_i1j1) + (E_x_int[x_i][y_i+1]*f_ij1)
            E_Par_y = (E_y_int[x_i][y_i]*f_ij) + (E_y_int[x_i+1][y_i]*f_i1j) + (E_y_int[x_i+1][y_i+1]*f_i1j1) + (E_y_int[x_i][y_i+1]*f_ij1)
            B_z_Par = (B_z_n[x_i][y_i]*f_ij) + (B_z_n[x_i+1][y_i]*f_i1j) + (B_z_n[x_i+1][y_i+1]*f_i1j1) + (B_z_n[x_i][y_i+1]*f_ij1)

            h = (Par[i][4]/Par[i][5])*dt
            k = 1/(1+(.25*(h**2)*(B_z_Par**2)))

            Par[i][2] = k*(((1-(.25*(h**2)*(B_z_Par**2)))*u_x_mhalf) + (h*((u_y_mhalf*B_z_Par)+E_Par_x)) + .5*(h**2)*(E_Par_y*B_z_Par))
            Par[i][3] = k*(((1-(.25*(h**2)*(B_z_Par**2)))*u_y_mhalf) + (h*((-1*(u_x_mhalf*B_z_Par))+E_Par_y)) - .5*(h**2)*(E_Par_x*B_z_Par))

            #Step 3: Claculate current density J^n+1/2

            J_x_int[x_i][y_i] += Par[i][2]*((Par[i][4]/(dx*dy))*f_ij)
            J_x_int[x_i + 1][y_i] += Par[i][2]*((Par[i][4]/(dx*dy))*f_i1j)
            J_x_int[x_i + 1][y_i + 1] += Par[i][2]*((Par[i][4]/(dx*dy))*f_i1j1)
            J_x_int[x_i][y_i + 1] += Par[i][2]*((Par[i][4]/(dx*dy))*f_ij1)

            J_y_int[x_i][y_i] += Par[i][3]*((Par[i][4]/(dx*dy))*f_ij)
            J_y_int[x_i + 1][y_i] += Par[i][3]*((Par[i][4]/(dx*dy))*f_i1j)
            J_y_int[x_i + 1][y_i + 1] += Par[i][3]*((Par[i][4]/(dx*dy))*f_i1j1)
            J_y_int[x_i][y_i + 1] += Par[i][3]*((Par[i][4]/(dx*dy))*f_ij1)

            J_x_half = (J_x_int+np.roll(J_x_int,1,axis=0))/2
            J_y_half = (J_y_int+np.roll(J_y_y,1,axis=1))/2

    #Produce visuals for each time step
    if vis.write_movie is False:

        if t % 1==0:

            Par_x_plot = [i[0] for i in Par]
            Par_y_plot = [i[1] for i in Par]
            Par_vx_plot = [i[2] for i in Par]
            Par_vy_plot = [i[3] for i in Par]
            Par_q_plot = [i[4] for i in Par]

            Par_x_plot = np.asarray(Par_x_plot)
            Par_y_plot = np.asarray(Par_y_plot)

            Par_x_plot = (Par_x_plot/R_j)
            Par_y_plot = (Par_y_plot/R_j)

            #phi_plot = np.delete(phi,-1,axis=0)
            #phi_plot = np.delete(phi_plot,-1,axis=1)
            #E_mag_plot = np.delete(E_mag,-1,axis=0)
            #E_mag_plot = np.delete(E_mag_plot,-1,axis=1)

            #q_plot = np.delete(q_grid,-1,axis=0)
            #q_plot = np.delete(q_plot,-1,axis=1)

            if opt.DebyeTest is True:

                D_Par_x_plot = Par_x_plot[-1]
                D_Par_y_plot = Par_y_plot[-1]

                del Par_x_plot[-1]
                del Par_y_plot[-1]
                del Par_vx_plot[-1]
                del Par_vy_plot[-1]
                del Par_q_plot[-1]

            if grid.CartesianGrid is True:
                plt.clf()
                plt.subplot(2,1,1)
                plt.contourf(X,Y,E_mag)
                plt.colorbar(label = '$|E|$')
                plt.title('$t$ = {}, Particle Number = {}'.format(round(t*dt,2),T_PPC))
                #plt.quiver(X,Y,Par_vx_plot,Par_vy_plot,color='white')
                #plt.plot(x,-(phi_plot[49]))
                plt.ylabel('$y$ $(R_j)$')

                plt.subplot(2,1,2)
                plt.scatter(Par_x_plot,Par_y_plot, c=Par_q_plot, cmap=cm.seismic)
                #plt.scatter(Par_x_plot,Par_y_plot, c='blue')
                plt.colorbar(label = '$q$')
                plt.xlabel('$x$ $(R_j)$')
                plt.ylabel('$y$ $(R_j)$')

                if opt.DebyeTest is True:
                    plt.scatter(D_Par_x_plot,D_Par_y_plot)

                plt.xlim((x_min/R_j),(x_max/R_j))
                plt.ylim((y_min/R_j),(y_max/R_j))
                plt.draw()
                plt.pause(0.001)

            elif grid.CylindricalGrid is True:

                if opt.RegenerateParticles is True:

                    Moon_Par_x_plot = Par_x_plot[-1]
                    Moon_Par_y_plot = Par_y_plot[-1]

                    del Par_x_plot[-1]
                    del Par_y_plot[-1]
                    del Par_vx_plot[-1]
                    del Par_vy_plot[-1]
                    del Par_q_plot[-1]

                ax[0].cla()
                ax[1].cla()

                ax[0].scatter(Par_y_plot,Par_x_plot)
                if opt.RegenerateParticles is True:
                    ax[0].scatter(Moon_Par_y_plot,Moon_Par_x_plot,c='green')
                #ax[0].set_title('$t$ = {}'.format(t*dt,2))
                ax[0].set_rlim(0,r_max)

                ax[1].contourf(Theta, R, E_mag)
                ax[1].set_rlim(0,r_max)

                fig.suptitle('$t$ = {}, particles = {}'.format(round(t*dt,2),T_PPC))

                plt.draw()
                plt.pause(0.001)

    else:
        #Write gif
        # Draw the plot and then make sure that the canvas renderer draws it.
        Par_x_plot = [i[0] for i in Par]
        Par_y_plot = [i[1] for i in Par]
        Par_vx_plot = [i[2] for i in Par]
        Par_vy_plot = [i[3] for i in Par]
        Par_q_plot = [i[4] for i in Par]

        #Par_x_plot = np.asarray(Par_x_plot)
        #Par_y_plot = np.asarray(Par_y_plot)

        #Par_x_plot = (Par_x_plot/R_j)
        #Par_y_plot = (Par_y_plot/R_j)

        #E_mag_plot = np.delete(E_mag,-1,axis=0)
        #E_mag_plot = np.delete(E_mag_plot,-1,axis=1)

        if grid.CartesianGrid is True:

            fig.suptitle('$t$ = {}, Particle Number = {}'.format(round(t*dt,2),T_PPC))
            scat = ax[0].contourf(X,Y,E_mag)
            #cbar = plt.colorbar(scat,ax=ax[0])
            #cbar.set_label('$|E|$')
            ax[0].set_ylabel('$y$ $(R_j)$')
            #ax.quiver(X,Y,Par_vx_plot,Par_vy_plot,color='white')
            #ax1 = fig.add_subplot(gs[1,0])
            ax[0].set_xlim((x_min/R_j),(x_max/R_j))
            ax[0].set_ylim((y_min/R_j),(y_max/R_j))
            ax[1].set_xlim((x_min/R_j),(x_max/R_j))
            ax[1].set_ylim((y_min/R_j),(y_max/R_j))
            scat = ax[1].scatter(Par_x_plot,Par_y_plot, c=Par_q_plot, cmap=cm.seismic)
            #cbar = plt.colorbar(scat,ax=ax[1])
            #cbar.set_label('$q$')
            ax[1].set_xlabel('$x$ $(R_j)$')
            ax[1].set_ylabel('$y$ $(R_j)$')
            fig.canvas.draw()

        if grid.CylindricalGrid is True:
            if opt.RegenerateParticles is True:

                Moon_Par_x_plot = Par_x_plot[-1]
                Moon_Par_y_plot = Par_y_plot[-1]

                del Par_x_plot[-1]
                del Par_y_plot[-1]
                del Par_vx_plot[-1]
                del Par_vy_plot[-1]
                del Par_q_plot[-1]

            ax[0].scatter(Par_y_plot,Par_x_plot)
            if opt.RegenerateParticles is True:
                ax[0].scatter(Moon_Par_y_plot,Moon_Par_x_plot,c='green')
            ax[0].set_rlim(0,r_max)

            ax[1].contourf(Theta, R, E_mag)
            ax[1].set_rlim(0,r_max)

            fig.suptitle('$t$ = {}, particles = {}'.format(round(t*dt,2),T_PPC))
            fig.canvas.draw()

        # Grab the canvas buffer as a set of pixels, convert them to a Numpy array,
        # then reshape the 1D set of pixels into an h x w x 3 (RGB) array. Then write
        # to the movie file.
        image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8').reshape((int(h),int(w),3))
        writer.append_data(image)


if vis.write_movie is True:

    #End video
    writer.close()
