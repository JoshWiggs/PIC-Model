import numpy as np
import random as ran
import scipy.constants as con
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#Options
class options:

    def __init__(self):
            self.DebyeTest = bool(True)
            self.RegenerateParticles = bool(False)

    class visual:

        def __init__(self):
            self.write_movie = bool(False)

opt = options()
vis = options.visual()

class maxwell_velocities(stats.rv_continuous):
    def _pdf(self,v):
        return (((con.m_e)/(2*con.pi*con.k*300))**(3/2))*(4*con.pi*(v**2))*np.exp(-((con.m_e*(v)**2)/(2*con.k*300)))

velocity_pdf = maxwell_velocities(name='velocity_pdf', a=0.0, b=con.c)

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
x_min = 0
x_max = 10
y_min = 0
y_max = 10
t_min = 0
t_max = 2.0
v_min = -1
v_max = 1
q_c = [con.e,-(10*con.e)] #particle charge
m = [con.m_p,(10*con.m_e)] #particle mass

nx = 50 #Number of steps taken from y_min to y_max
ny = 50 #Number of steps taken from x_min to x_max
nt = 200 #Number of time steps
n_smoothing = 100 #For differencing methods that require smoothing
PPC = 0.2 #Number of particles per cell

phi = np.zeros(((nx + 1),(ny + 1))) #Electric potential
E_x = np.zeros(((nx + 1),(ny + 1))) #x component of the Electric Field
E_y = np.zeros(((nx + 1),(ny + 1))) #y component of the Electric Field
E_mag = np.zeros((nx,ny)) #magnitude of the Electric field

dx = (abs(x_min) + x_max) / nx #Size of step in x domain
dy = (abs(y_min) + y_max) / ny #Size of step in y domain
dt = (abs(t_min) + t_max) / nt #Size of step on the temporal domain
T_c = nx * ny #Total number of cells
T_PPC = int(T_c * PPC) #Total number of particles in simulation

#Populate dimesnional vectors
x = np.linspace(x_min,x_max,num=nx)
y = np.linspace(y_min,y_max,num=ny)
t = np.linspace(t_min,t_max,num=nt)

#Meshgrid for plotting
X,Y = np.meshgrid(x,y)

#Create list of particles properties
Par_x = np.zeros(T_PPC)
Par_y = np.zeros(T_PPC)
Par_vx = np.zeros(T_PPC)
Par_vy = np.zeros(T_PPC)
Par_q = np.zeros(T_PPC)
Par_m = np.zeros(T_PPC)

if vis.write_movie is True:

    #Make a figure, get the axes
    fig = plt.Figure()
    ax = fig.gca()
    gs = gridspec.GridSpec(2, 1, figure=fig)

    # Set the figure canvas to the Agg backend and get the width/height of the canvas (in pixels)
    fig.canvas = matplotlib.backends.backend_agg.FigureCanvasAgg(fig)
    w,h = fig.canvas.get_renderer().get_canvas_width_height()

    # Setup the imageio writer.
    writer = imageio.get_writer('regen_electron_proton.gif', fps=25)

print('Initialising {} particles...'.format(T_PPC))

for i in range(0,T_PPC):

    #TODO: Generalise
    j = ran.randint(0,int(len(q_c)-1)) #Random integer for mass and charge selection

    Par_x[i] = ran.uniform(x_min,x_max) #Particle x coordinate
    Par_y[i] = ran.uniform(y_min,y_max) #Particle y coordinate
    Par_vx[i] = ran.uniform(v_min,v_max) #Particle velocity on the x direction
    Par_vy[i] = ran.uniform(v_min,v_max) #Particle velocity on the y direction
    #Par_vx[i] = (velocity_pdf.rvs() * 0.001)  #Particle velocity on the x direction
    #Par_vy[i] = (velocity_pdf.rvs() * 0.001) #Particle velocity on the y direction
    Par_q[i] = q_c[j] #Particle charge
    Par_m[i] = m[j] #Particle mass

#Join particle position lists into one array
#TODO: Rework this as data types produced seem to make calcultions difficult
Par = list(zip(Par_x,Par_y,Par_vx,Par_vy,Par_q,Par_m))

#If the Debye test is active then generate the Debye test particle and
#append it as the final item in the particle list
if opt.DebyeTest is True:

    print('Creating Debye particle...')

    D_Par_x = ((x_min + x_max) / 2) + (dx / 2)
    D_Par_y = ((y_min + y_max) / 2) + (dy / 2)
    D_Par_vx = 0.0
    D_Par_vy = 0.0
    D_Par_q = (con.e * 10000)
    D_Par_m = (con.m_p * 1)

    D_Par = [D_Par_x,D_Par_y,D_Par_vx,D_Par_vy,D_Par_q,D_Par_m]

    Par.append(D_Par)

print('Initialising particle simulation...')

#Numerical loop responsible for interating the particles through time
for t in range(0,nt):

    #Activate condition only if writing gif
    if vis.write_movie is True:

        #Clear the axes so we can draw the new frame
        ax.cla()

    #Print index 't' every 10 time steps to show progress
    if t % 10 == 0:

        print(t)

    #Perpare index 'i' for while loop responsible for removing particles which
    #move outside of the simulation region
    i = 0

    #While loop to remove particles that move outside of simulation area. This
    #is done by calculating grid positions and deleting particles that lie on
    #a grid ID outside the area.
    while i < T_PPC:

        #Open / clear variables for use
        x_i = 0
        y_i = 0

        #Calculate grid locations of particle
        x_i = np.floor(Par[i][0] / dx)
        y_i = np.floor(Par[i][1] / dy)

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
            for r in range(0,Regenerate_num):

                #This works in the same way as the loop in the pre-process
                j = ran.randint(0,int(len(q_c)-1))

                Regen_par_x = ran.uniform(x_min,x_max)
                Regen_par_y = ran.uniform(y_min,y_max)
                Regen_par_vx = ran.uniform(v_min,v_max)
                Regen_par_vy = ran.uniform(v_min,v_max)
                Regen_par_q = q_c[j]
                Regen_par_m = m[j]

                Regen_par = [Regen_par_x,Regen_par_y,Regen_par_vx,Regen_par_vy,Regen_par_q,Regen_par_m]

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
    g_info = np.zeros((T_PPC,4)) #Grid information
    q_grid= np.zeros(((nx + 1),(ny + 1))) #Charge density

    #Step 1: Compute charge density on grid using bilinear interpolation interpretation
    #(first-order weighting) from particle locations
    for i in range(0,T_PPC):

        x_i = 0 #Particle x grid position
        y_i = 0 #Particle y grid position
        hx = 0 #Partcle distance from its grid postion in the x domain
        hy = 0 #Partcle distance from its grid postion in the y domain

        #Calculate grid locations of particle
        x_i = int(np.floor(Par[i][0] / dx))
        y_i = int(np.floor(Par[i][1] / dy))

        #Calculate distance of particle from grid location
        #DEBUG: If a particle outside the simulation  area has not been removed
        #       raise error and return particle information
        try:
            hx = Par[i][0] - x[x_i]
            hy = Par[i][1] - y[y_i]
        except:
            raise Exception('ERROR: An error occurred during timestep {} check particle {}'.format(t,i))


        #Store grid information relating to each particle for use later
        g_info[i][0] = x_i
        g_info[i][1] = y_i
        g_info[i][2] = hx
        g_info[i][3] = hy

        #Calculate contribution of each particle's charge to surrounding grid
        #locations in order to calculate charge density on grid
        #TODO: add charge as a particle parameter
        q_grid[x_i][y_i] += Par[i][4]*((dx - hx) * (dy - hy) / (dx * dy))
        q_grid[x_i + 1][y_i] += Par[i][4]*((hx * (dy - hy)) / (dx * dy))
        q_grid[x_i + 1][y_i + 1] += Par[i][4]*((hx * hy) / (dx * dy))
        q_grid[x_i][y_i + 1] += Par[i][4]*(((dx - hx)*hy) / (dx*dy))

    #Step 2: Compute Electric Potential
    for n in range(0,n_smoothing):

        """
        for i in range(1,nx):
            for j in range(1,ny):

                phi[i][j] = (q_grid[i][j]*((dx**2)*(dy**2)) + ((dy**2)*
                (phi[i-1][j] +phi[i+1][j])) + ((dx**2)*(phi[i][j-1]+phi[i][j+1]))
                 / (2*((dy**2) +(dy**2))))
        """

        phi = (q_grid*((dx**2)*(dy**2)) + ((dy**2)*(np.roll(phi,1,axis=0)
        + np.roll(phi,-1,axis=0))) + ((dx**2)*(np.roll(phi,1,axis=1)
        + np.roll(phi,-1,axis=1))) / (2*((dy**2) +(dy**2))))

    #Step 3: Compute Electric Field

    """
    for i in range(1,nx):
        for j in range(1,ny):

            E_x[i][j] = ((phi[i + 1][j]) - phi[i - 1][j]) / (2 * dx)

            E_y[i][j] = ((phi[i][j + 1]) - phi[i][j - 1]) / (2 * dy)

            E_mag[i][j] = np.sqrt((E_x[i][j]**2) + (E_y[i][j]**2))
    """

    E_x = ((np.roll(phi,-1,axis=0)) - np.roll(phi,1,axis=0)) / (2 * dx)

    E_y = ((np.roll(phi,-1,axis=1)) - np.roll(phi,1,axis=1)) / (2 * dy)

    E_mag = np.sqrt((E_x**2) + (E_y**2))
    E_mag = np.delete(E_mag,-1,axis=0)
    E_mag = np.delete(E_mag,-1,axis=1)


    #Step 4: Move particles
    for i in range(0,T_PPC):

        #Safety check for particle information data type
        Par[i] = list(Par[i])

        #Access particle information and the grid information relating to they
        #particle for use in calculating particle motion
        x_i = int(g_info[i][0])
        y_i = int(g_info[i][1])
        u_x_old = float(Par[i][2])
        u_y_old = float(Par[i][3])
        hx = float(g_info[i][2])
        hy = float(g_info[i][3])

        #Clear variables for reuse
        E_Par_x = float(0) #Electric field felt by the particle in the x domain
        E_Par_y = float(0) #Electric field felt by the particle in the y domain

        #Calculate the electric field strength at the particle from each of the
        #surrounding grid points using bilinear interpolation interpretation
        E_Par_x = (E_x[x_i][y_i]*(dx*dy)/((dx-hx)*(dy-hy))) + (E_x[x_i+1][y_i]*((dx*dy)/(hx*(dy-hy)))) + (E_x[x_i+1][y_i+1])*((dx*dy)/(hx*hy)) + (E_x[x_i][y_i+1]*((dx*dy)/(dx-hx)*hy))

        E_Par_y = (E_y[x_i][y_i]*(dx*dy)/((dx-hx)-(dy-hy))) + (E_y[x_i+1][y_i]*((dx*dy)/(hx*(dy-hy)))) + (E_y[x_i+1][y_i+1])*((dx*dy)/(hx*hy)) + (E_y[x_i][y_i+1]*((dx*dy)/(dx-hx)*hy))

        #Using updated electic field calculate updated particle velocities
        #TODO: missing charge and mass from calculation
        u_x = u_x_old - (((dt * Par[i][4]) / Par[i][5]) * E_Par_x)

        u_y = u_y_old - (((dt * Par[i][4]) / Par[i][5]) * E_Par_y)

        ##### If Debye physical test is being conducted then the final particle
        ##### in the list is the stationary Debye test particle. Therefore, the
        ##### velocities are set to zero to prevent motion
        if opt.DebyeTest is True:
            if i == int(T_PPC-1):

                u_x = 0.0
                u_y = 0.0

        #Update paricle velocity
        Par[i][2] = u_x
        Par[i][3] = u_y

        #Use particle velocity to calculate new particle position
        Par[i][0] += dt*u_x
        Par[i][1] += dt*u_y

    #Produce visuals for each time step
    if vis.write_movie is False:

        if t % 1==0:

            Par_x_plot = [i[0] for i in Par]
            Par_y_plot = [i[1] for i in Par]
            Par_vx_plot = [i[2] for i in Par]
            Par_vy_plot = [i[3] for i in Par]
            Par_q_plot = [i[4] for i in Par]

            if opt.DebyeTest is True:
                del Par_x_plot[-1]
                del Par_y_plot[-1]
                del Par_vx_plot[-1]
                del Par_vy_plot[-1]
                del Par_q_plot[-1]

            plt.clf()
            plt.subplot(2,1,1)
            plt.contourf(X,Y,E_mag)
            plt.colorbar(label = '$|E|$')
            plt.title('$t$ = {}, Particle Number = {}'.format(round(t*dt,2),T_PPC))
            #plt.quiver(X,Y,Par_vx_plot,Par_vy_plot,color='white')
            plt.subplot(2,1,2)
            plt.scatter(Par_x_plot,Par_y_plot, c=Par_q_plot, cmap=cm.seismic)
            plt.colorbar(label = '$q$')
            plt.xlim(x_min,x_max)
            plt.ylim(y_min,y_max)
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
        #ax = fig.add_subplot(gs[0,0])
        #ax.contourf(X,Y,E_mag)
        #ax.colorbar(label = '$|E|$')
        ax.set_title('$t$ = {}, Particle Number = {}'.format(round(t*dt,2),T_PPC))
        #ax.quiver(X,Y,Par_vx_plot,Par_vy_plot,color='white')
        #ax1 = fig.add_subplot(gs[1,0])
        ax.scatter(Par_x_plot,Par_y_plot, c=Par_q_plot, cmap=cm.seismic)
        #ax.colorbar(label = '$q$')
        ax.set_xlim(x_min,x_max)
        ax.set_ylim(y_min,y_max)
        fig.canvas.draw()

        # Grab the canvas buffer as a set of pixels, convert them to a Numpy array,
        # then reshape the 1D set of pixels into an h x w x 3 (RGB) array. Then write
        # to the movie file.
        image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8').reshape((int(h),int(w),3))
        writer.append_data(image)

if vis.write_movie is True:

    #End video
    writer.close()

"""
###### PLOTTING #######
new_par_x  = np.zeros(T_PPC)
new_par_y  = np.zeros(T_PPC)

for i in range(0,T_PPC):

    new_par_x[i] = Par[i][0]
    new_par_y[i] = Par[i][1]

plt.clf()
plt.subplot(2,1,1)
plt.scatter(Par_x,Par_y)
plt.subplot(2,1,2)
plt.scatter(new_par_x,new_par_y)
plt.show()
"""
