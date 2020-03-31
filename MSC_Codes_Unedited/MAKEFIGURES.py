#RUN ME IN JUPYTER NOTEBOOK, SEPARATE OUT THE FUNCTIONS FROM THE MAIN CODE AT THE END

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import colorsys
import segyio #Not standard with python https://github.com/equinor/segyio

############################
### ADDITIONAL FUNCTIONS ###
############################

class MidpointNormalize(colors.Normalize):
    """
    Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
    
def RayCalculationProgress(i, count):
    """
    My first attempt at some kind of a loading bar.
    """
    print "Ray density calculation progress: "
    flags = dict.fromkeys(["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"], "no")
    if i/count <= 0.1 and flags["a"] == "no":
        flags.update(dict.fromkeys(["a"], "yes"))
        print "|~~~~~                                                            |"
    if i/count <= 0.2 and flags["b"] == "no":
        flags.update(dict.fromkeys(["b"], "yes"))
        print "|~~~~~*~~~~~                                                      |"
    if i/count <= 0.3 and flags["c"] == "no":
        flags.update(dict.fromkeys(["c"], "yes"))
        print "|~~~~~*~~~~~*~~~~~                                                |"
    if i/count <= 0.4 and flags["d"] == "no":
        flags.update(dict.fromkeys(["d"], "yes"))
        print "|~~~~~*~~~~~*~~~~~*~~~~~                                          |"
    if i/count <= 0.5 and flags["e"] == "no":
        flags.update(dict.fromkeys(["e"], "yes"))
        print "|~~~~~*~~~~~*~~~~~*~~~~~*~~~~~                                    |"
    if i/count <= 0.6 and flags["f"] == "no":
        flags.update(dict.fromkeys(["f"], "yes"))
        print "|~~~~~*~~~~~*~~~~~*~~~~~*~~~~~*~~~~~                              |"
    if i/count <= 0.7 and flags["g"] == "no":
        flags.update(dict.fromkeys(["g"], "yes"))
        print "|~~~~~*~~~~~*~~~~~*~~~~~*~~~~~*~~~~~*~~~~~                        |"
    if i/count <= 0.8 and flags["h"] == "no":
        flags.update(dict.fromkeys(["h"], "yes"))
        print "|~~~~~*~~~~~*~~~~~*~~~~~*~~~~~*~~~~~*~~~~~*~~~~~                  |"
    if i/count <= 0.9 and flags["i"] == "no":
        flags.update(dict.fromkeys(["i"], "yes"))
        print "|~~~~~*~~~~~*~~~~~*~~~~~*~~~~~*~~~~~*~~~~~*~~~~~*~~~~~            |"
    if i/count <= 1.0 and flags["j"] == "no":
        flags.update(dict.fromkeys(["j"], "yes"))
        print "|~~~~~*~~~~~*~~~~~*~~~~~*~~~~~*~~~~~*~~~~~*~~~~~*~~~~~*~~~~~      |"

################################################
### RAY DENSITY SHADING, NOT YET IMPLIMENTED ###
################################################
        
#I'm going to HAVE to correct for bathymetry (I've shifted OBS down to make it work)
def ray_density_plot(raypaths, nodes): #, bathymetry
    """
    Input a TOMO2D ray path file and shade out regions with poor ray coverage.
    """
    
    #Read text file
    raypaths = np.genfromtxt(raypaths)
    
    #Count the number of rays and save indicies
    count = 0
    rayidx = []
    for i in range(np.shape(raypaths)[0]):
        if raypaths[i,0] == 0 and raypaths[i,1] == 0:
            count += 1
            rayidx.append(i)

    print "Number of rays: ", count
    
    #Plot for ray densities
    raydensity = np.zeros((len(zcorr),len(xcorr)))

    #Iterate through each ray
    for i in range(len(rayidx)-1):
        #Remove duplicates; not sure why there's duplicates produced through the ray-tracing process
        raysegmentpaths = np.unique(raypaths[rayidx[i]+1:rayidx[i+1],:],axis=0)
        #Iterate through each point within a raypath
        for j in range(len(raysegmentpaths[:,0])):
            #Find the node closest to the current point
            xnode = min(range(len(xcorr)), key=lambda x: abs(xcorr[x]-raysegmentpaths[j,0]))
            znode = min(range(len(zcorr)), key=lambda z: abs(zcorr[z]-raysegmentpaths[j,1]))

            #Distribute a weight to the surrounding four nodes with a four point stencil

            #Figure out where the other x node is
            if xcorr[xnode] > raysegmentpaths[j,0]:
                altxnode = xnode - 1
            elif xcorr[xnode] < raysegmentpaths[j,0]:
                altxnode = xnode + 1
            else:
                altxnode = "exact"

            #Figure out where the other z node is
            if zcorr[znode] > raysegmentpaths[j,1]:
                altznode = znode - 1
            elif zcorr[znode] < raysegmentpaths[j,1]:
                altznode = znode + 1
            else:
                altznode = "exact"

            #Define the nodes for the four point stencil, scenario considered where ray lies on a node &
            #Calculate contribution weights for each node as a precentage of linear distance from the ray
            #Normalize such that the highest contribution is due to raypaths close to nodes
            #Continually add to 
            if altxnode != "exact" and altznode != "exact":
                #Points
                point = [xnode,znode]
                point1 = [xnode,altznode]
                point2 = [altxnode,znode]
                point3 = [altxnode,altznode]
                #Weights
                weight = sqrt(abs(xcorr[point[0]]-raysegmentpaths[j,0])**2+abs(zcorr[point[1]]-raysegmentpaths[j,1])**2)
                weight1 = sqrt(abs(xcorr[point1[0]]-raysegmentpaths[j,0])**2+abs(zcorr[point1[1]]-raysegmentpaths[j,1])**2)
                weight2 = sqrt(abs(xcorr[point2[0]]-raysegmentpaths[j,0])**2+abs(zcorr[point2[1]]-raysegmentpaths[j,1])**2)
                weight3 = sqrt(abs(xcorr[point3[0]]-raysegmentpaths[j,0])**2+abs(zcorr[point3[1]]-raysegmentpaths[j,1])**2)
                #Normalize weights
                normtotalweight = (1 - weight) + (1 - weight1) + (1 - weight2) + (1 - weight3)
                weight_pct = (1 - weight)/normtotalweight
                weight1_pct = (1 - weight1)/normtotalweight
                weight2_pct = (1 - weight2)/normtotalweight
                weight3_pct = (1 - weight3)/normtotalweight
                #Add to raydensity plot
                raydensity[point[1],point[0]] = raydensity[point[1],point[0]] + weight_pct
                raydensity[point1[1],point1[0]] = raydensity[point1[1],point1[0]] + weight1_pct
                raydensity[point2[1],point2[0]] = raydensity[point2[1],point2[0]] + weight2_pct
                raydensity[point3[1],point3[0]] = raydensity[point3[1],point3[0]] + weight3_pct
                #Flag
                flag = "Uncentered"
            elif altxnode == "exact" and altznode != "exact":
                #Points
                point = [xnode,znode]
                point1 = [xnode,altznode]
                #Weights
                weight = sqrt(abs(xcorr[point[0]]-raysegmentpaths[j,0])**2+abs(zcorr[point[1]]-raysegmentpaths[j,1])**2)
                weight1 = sqrt(abs(xcorr[point1[0]]-raysegmentpaths[j,0])**2+abs(zcorr[point1[1]]-raysegmentpaths[j,1])**2)
                #Normalize weights
                totalweight = weight + weight1 
                weight_pct = weight/totalweight
                weight1_pct = weight1/totalweight
                normtotalweight = (1 - weight) + (1 - weight1)
                weight_pct = (1 - weight)/normtotalweight
                weight1_pct = (1 - weight1)/normtotalweight
                #Add to raydensity plot
                raydensity[point[1],point[0]] = raydensity[point[1],point[0]] + weight_pct
                raydensity[point1[1],point1[0]] = raydensity[point1[1],point1[0]] + weight1_pct
                #Flag
                flag = "Centeredx"
            elif altxnode != "exact" and altznode == "exact":
                #Points
                point = [xnode,znode]
                point1 = [altxnode,znode]
                #Weights
                weight = sqrt(abs(xcorr[point[0]]-raysegmentpaths[j,0])**2+abs(zcorr[point[1]]-raysegmentpaths[j,1])**2)
                weight1 = sqrt(abs(xcorr[point1[0]]-raysegmentpaths[j,0])**2+abs(zcorr[point1[1]]-raysegmentpaths[j,1])**2)
                #Normalize weights
                totalweight = weight + weight1 
                weight_pct = weight/totalweight
                weight1_pct = weight1/totalweight
                normtotalweight = (1 - weight) + (1 - weight1)
                weight_pct = (1 - weight)/normtotalweight
                weight1_pct = (1 - weight1)/normtotalweight
                #Add to raydensity plot
                raydensity[point[1],point[0]] = raydensity[point[1],point[0]] + weight_pct
                raydensity[point1[1],point1[0]] = raydensity[point1[1],point1[0]] + weight1_pct
                #Flag
                flag = "Centeredz"
            else:
                #Point
                point = [xnode,znode]
                #Weight
                weight_pct = 1.0
                #Add to raydensity plot
                raydensity[point[1],point[0]] = weight_pct
                #Flag
                flag = "Centered"
        RayCalculationProgress(i, count)
    
    print "|~~~~~*~~~~~*~~~~~*~~~~~*~~~~~*~~~~~*~~~~~*~~~~~*~~~~~*~~~~~*~~~~~|"
    
    return raydensity

#Convert a GMT cpt file to Matplotlib colormap

def gmtColormap(fileName):
    try:
        f = open(fileName)
    except:
        print "file ",fileName, "not found"
        return None

    lines = f.readlines()

    x = []
    r = []
    g = []
    b = []
    colorModel = "RGB"
    for l in lines:
        ls = l.split() 
        x.append(float(ls[0]))
        r.append(float(ls[1]))
        g.append(float(ls[2]))
        b.append(float(ls[3]))
        xtemp = float(ls[4])
        rtemp = float(ls[5])
        gtemp = float(ls[6])
        btemp = float(ls[7])

    x.append(xtemp)
    r.append(rtemp)
    g.append(gtemp)
    b.append(btemp)

    nTable = len(r)
    x = np.array( x , np.float32)
    r = np.array( r , np.float32)
    g = np.array( g , np.float32)
    b = np.array( b , np.float32)
    if colorModel == "HSV":
        for i in range(r.shape[0]):
            rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
            r[i] = rr ; g[i] = gg ; b[i] = bb
    if colorModel == "HSV":
        for i in range(r.shape[0]):
            rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
            r[i] = rr ; g[i] = gg ; b[i] = bb
    if colorModel == "RGB":
        r = r/255.
        g = g/255.
        b = b/255.
    xNorm = (x - x[0])/(x[-1] - x[0])

    red = []
    blue = []
    green = []
    for i in range(len(x)):
        red.append([xNorm[i],r[i],r[i]])
        green.append([xNorm[i],g[i],g[i]])
        blue.append([xNorm[i],b[i],b[i]])
    colorDict = {"red":red, "green":green, "blue":blue}
    return (colorDict)

#################################################
### PLOT THREE MODELS, PLOT THREE DIFFERENCES ###
#################################################

def plt_3synth_mdl(model_1, model_2, model_3, truemdl, nodes, bathy, src_loc_full, src_loc_sparse, src_loc_usparse,
                   spacing, vrange_vel, vrange_diff, unit="[m]", fontscale=1, aspect=1,
                   shift=[0,0], letterlbl=None, size=[72,48], cmap="jet", savepath=None, savename=None):

    """ 
    This code is intended to generate a high quality image of one model, and the two difference 
    models.
    """
    
    #########################
    ### READ BINARY FILES ###
    #########################

    #Open Model 1
    with open(model_1, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    mesh_1 = np.reshape(data, [nodes[0], nodes[1]])
    
    #Open Model 2
    with open(model_2, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    mesh_2 = np.reshape(data, [nodes[0], nodes[1]])
    
    #Open Model 3
    with open(model_3, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    mesh_3 = np.reshape(data, [nodes[0], nodes[1]])
    
    #Open True Model
    with open(truemdl, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    truemdl = np.reshape(data, [nodes[0], nodes[1]])

    #Read the bathymetry file
    with open(bathy, 'rb') as f:
        bathy = np.fromfile(f, dtype=np.float32)
    if unit=="[km]":
        for i in range(len(bathy)):
            bathy_mesh = bathy/1000.0
    else:
        bathy_mesh = bathy
    
    ###################
    ### Plot Models ###
    ###################

    #Establish units, default is meters
    if unit == "[km]":
        spacing = spacing/1000.
        
    #Subtract the true model from each mesh to formulate residuals
    diff_1 = mesh_1 - truemdl
    diff_2 = mesh_2 - truemdl
    diff_3 = mesh_3 - truemdl
        
    #Define Model Dimensions
    nx_min = 0 * spacing
    nx_max = nodes[0] * spacing
    nz_min = nodes[1] * spacing
    nz_max = 0 * spacing

    #Colorbar range for difference plots
    elev_min, elev_max = vrange_diff[0], vrange_diff[1]
    mid_val=0
    #Colorbar range for velocity model
    elev_min_vel, elev_max_vel = vrange_vel[0], vrange_vel[1]
    
    #Setup subplots and define figure size
    fig = plt.figure(figsize=(size[0], size[1]))
    gs = gridspec.GridSpec(nrows=3, ncols=2,)

    #Plot 1
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=2)
    ax1.scatter(src_loc_full[:,0]*1000.0, src_loc_full[:,1]*1000.0, s=60*fontscale, color="white", marker="v")
    im = ax1.imshow(mesh_1.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min_vel, vmax=elev_max_vel, 
                    cmap=cmap, aspect=aspect)
    #ax1.set_title('Sparse Acquisition',fontweight="bold", size=20.0*fontscale) # Title
    ax1.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    #ax1.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax1.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax1.tick_params(axis='y', which='both', labelsize=16.0*fontscale)
    #ax1.grid()
    
    #Difference Plot 1a
    ax1a = fig.add_subplot(gs[0, 1])
    ax1a.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=2)
    ax1a.scatter(src_loc_full[:,0]*1000.0, src_loc_full[:,1]*1000.0, s=60*fontscale, color="black", marker="v")
    dif = ax1a.imshow(diff_1.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min, vmax=elev_max, 
                     cmap='seismic', norm=MidpointNormalize(midpoint=mid_val), aspect=aspect)
    #ax1a.set_title('Residual Velocity Field',fontweight="bold", size=20.0*fontscale)
    #ax1a.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    #ax1a.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax1a.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax1a.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
    #ax1a.get_yaxis().set_visible(False)
    #ax1a.grid()
    
    #Plot 2
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=2)
    ax2.scatter(src_loc_sparse[:,0]*1000.0, src_loc_sparse[:,1]*1000.0, s=60*fontscale, color="white", marker="v")
    im = ax2.imshow(mesh_2.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min_vel, vmax=elev_max_vel, 
                    cmap=cmap, aspect=aspect)
    #ax2.set_title('Sparse Acquisition: Multiscale Inversion',fontweight="bold", size=20.0*fontscale) # Title
    ax2.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    #ax2.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax2.tick_params(axis='y', which='both', labelsize=16.0*fontscale)
    #ax2.grid()
    
    #Difference Plot 2a
    ax2a = fig.add_subplot(gs[1, 1])
    ax2a.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=2)
    ax2a.scatter(src_loc_sparse[:,0]*1000.0, src_loc_sparse[:,1]*1000.0, s=60*fontscale, color="black", marker="v")
    dif = ax2a.imshow(diff_2.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min, vmax=elev_max, 
                     cmap='seismic', norm=MidpointNormalize(midpoint=mid_val), aspect=aspect)
    #ax2a.set_title('Residual Velocity Field',fontweight="bold", size=20.0*fontscale)
    #ax2a.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    #ax2a.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax2a.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax2a.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
    #ax2a.get_yaxis().set_visible(False)
    #ax2a.grid()
    
    #Plot 3
    ax3 = fig.add_subplot(gs[2, 0])
    ax3.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=2)
    ax3.scatter(src_loc_usparse[:,0]*1000.0, src_loc_usparse[:,1]*1000.0, s=60*fontscale, color="white", marker="v")
    im = ax3.imshow(mesh_3.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min_vel, vmax=elev_max_vel, 
                    cmap=cmap, aspect=aspect)
    #ax3.set_title('Dense Acquisition',fontweight="bold", size=20.0*fontscale) # Title
    ax3.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    ax3.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax3.tick_params(axis='x', which='both', labelsize=16.0*fontscale)
    ax3.tick_params(axis='y', which='both', labelsize=16.0*fontscale)
    #ax3.grid()

    #Difference Plot 3a
    ax3a = fig.add_subplot(gs[2, 1])
    ax3a.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=2)
    ax3a.scatter(src_loc_usparse[:,0]*1000.0, src_loc_usparse[:,1]*1000.0, s=60*fontscale, color="black", marker="v")
    dif = ax3a.imshow(diff_3.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min, vmax=elev_max, 
                     cmap='seismic', norm=MidpointNormalize(midpoint=mid_val), aspect=aspect)
    #ax3a.set_title('Residual Velocity Field',fontweight="bold", size=20.0*fontscale)
    #ax3a.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    ax3a.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax3a.tick_params(axis='x', which='both', labelsize=16.0*fontscale)
    ax3a.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
    #ax3a.get_yaxis().set_visible(False)
    #ax3a.grid()

    #Positioning Adjustments
    fig.subplots_adjust(right=0.8)
    fig.subplots_adjust(hspace=0.5)

    #Colorbar velocities
    ax_cbar = fig.add_axes([0.0550, -0.05, 0.4555, 0.035]) #[0.46, 0.355, 0.0175, 0.3] VERT #[0.053, -0.05, 0.434, 0.035] HORIZ W YAX
    cbar = fig.colorbar(im, cax=ax_cbar, orientation='horizontal')
    cbar.ax.tick_params(labelsize=16.0*fontscale)
    cbar.set_label("Velocity (m/s)", fontsize=18.0*fontscale)

    #Colorbar difference
    ax_cbar1 = fig.add_axes([0.533, -0.05, 0.4555, 0.035]) #[0.9985, 0.355, 0.0175, 0.3] VERT #[0.56, -0.05, 0.434, 0.035] HORIZ W YAX
    cbar1 = fig.colorbar(dif, cax=ax_cbar1, orientation='horizontal')
    cbar1.ax.tick_params(labelsize=16.0*fontscale)
    cbar1.set_label("Velocity Difference (m/s)", fontsize=18.0*fontscale)

    if letterlbl != None:
        #Add text
        fig.text(0.21+shift[0], 0.885, letterlbl[0], fontweight="bold", fontsize=18*fontscale)
        fig.text(0.21+shift[0], 0.602, letterlbl[1], fontweight="bold", fontsize=18*fontscale)
        fig.text(0.21+shift[0], 0.319, letterlbl[2], fontweight="bold", fontsize=18*fontscale)
        
    if savepath != None and savename != None:
        plt.gcf().savefig(savepath+savename+".png", bbox_inches = 'tight', pad_inches = 0.0)
        
    plt.tight_layout(w_pad=1,h_pad=1)

    plt.show()
    
#########################
### PLOT SINGLE MODEL ###
#########################

def plot_onemdl(mdl, nodes, spacing, savepath=None, savename=None, vrange="true", 
                          unit="m", fontscale=1, shift=[0,0], letterlbl=None, size=[72,36], cmap='jet',
                         bathy=None, src_loc=None, Title=None):
    
    """ 
    This code is developed to gereate a high quality image of the starting model and true model for
    a synthetic inversion.
    """

    #########################
    ### READ BINARY FILES ###
    #########################

    #Read true model
    with open(mdl, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    smesh_1 = np.reshape(data, [nodes[0], nodes[1]])
    
    #Read the bathymetry file
    with open(bathy, 'rb') as f:
        bathy = np.fromfile(f, dtype=np.float32)
    if unit=="[km]":
        for i in range(len(bathy)):
            bathy_mesh = bathy/1000.0
    else:
        bathy_mesh = bathy

    ###################
    ### Plot Models ###
    ###################
    
    #Establish units, default is meters
    if unit == "[km]":
        spacing = spacing/1000.0
    
    #Define Model Dimensions
    nx_min = 0 * spacing
    nx_max = nodes[0] * spacing
    nz_min = nodes[1] * spacing
    nz_max = 0 * spacing

    #Define Figure Size
    size = size
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = size[0]
    fig_size[1] = size[1]
    plt.rcParams["figure.figsize"] = fig_size

    #Colorbar range, define limits based on starting or true model, or user specified
    if vrange == "true":
        elev_min, elev_max = min(smesh_1.flatten()), max(smesh_1.flatten())
    else:
        elev_min, elev_max = vrange[0], vrange[1]

    #Setup subplots
    fig = plt.figure() 
    ax1 = fig.add_subplot(111)

    #Plot
    if type(bathy) != str:
        ax1.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=4)
    if type(src_loc) != str:
        ax1.scatter(src_loc[:,0], src_loc[:,1], s=100*fontscale, color="black", marker="v")
    im = ax1.imshow(smesh_1.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min, vmax=elev_max, cmap=cmap,
                   aspect=2.0)
    if Title != None:
        ax1.set_title(Title,fontweight="bold", size=28.0*fontscale) # Title
    #else:
    #    ax1.set_title('Velocity Model',fontweight="bold", size=28.0*fontscale) # Title        
    ax1.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    ax1.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax1.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax1.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)

    #Colorbar
    ax_cbar = fig.add_axes([0.81+shift[1], 0.3, 0.03, 0.40])
    cbar = fig.colorbar(im, cax=ax_cbar, orientation='vertical')
    cbar.ax.tick_params(labelsize=16.0*fontscale)
    cbar.set_label("Velocity (m/s)", fontsize=18.0*fontscale)
    
    if letterlbl != None:
        #Add text
        fig.text(0.086+shift[0], 0.86, letterlbl[0], fontweight="bold", fontsize=18*fontscale)
        
    if savepath != None and savename != None:
        plt.gcf().savefig(savepath+savename+".png", bbox_inches = 'tight', pad_inches = 0.1)

    plt.show()
    
##########################################################
### PLOT TWO MODELS AND THE DIFFERENCE BETWEEN THE TWO ###
##########################################################

def plot_onemdl_diff(model1, model2, nodes, spacing, savepath=None, savename=None, overlay=None, 
                         vrange_diff=None, mid_val=0, unit="m", fontscale=1, shift=[0,0], 
                         letterlbl=None, size=[72,48], aspect=1.0,
                         cmap='jet', bathy=None, src_loc=None):

    """ 
    This code is intended to generate a high quality image of two models, and the difference between
    the two.
    """
    
    #########################
    ### READ BINARY FILES ###
    #########################

    #Open Model 1
    with open(model1, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    smesh_1 = np.reshape(data, [nodes[0], nodes[1]])

    #Open Model 2
    with open(model2, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    smesh_2 = np.reshape(data, [nodes[0], nodes[1]])
    
    #Read the bathymetry file
    with open(bathy, 'rb') as f:
        bathy = np.fromfile(f, dtype=np.float32)
    if unit=="[km]":
        for i in range(len(bathy)):
            bathy_mesh = bathy/1000.0
    else:
        bathy_mesh = bathy

    ###################
    ### Plot Models ###
    ###################

    #Establish units, default is meters
    if unit == "[km]":
        spacing = spacing/1000.0

    #Second model subtracted from first model as difference
    diff = smesh_2 - smesh_1

    #Define Model Dimensions
    nx_min = 0 * spacing
    nx_max = nodes[0] * spacing
    nz_min = nodes[1] * spacing
    nz_max = 0 * spacing

    #Colorbar range for difference plot, define limits based on min/max difference or user specified
    if vrange_diff == None:
        elev_min_diff, elev_max_diff = min(diff.flatten()), max(diff.flatten())
        mid_val=mid_val
    else:
        elev_min_diff, elev_max_diff = vrange_diff[0], vrange_diff[1]
        mid_val=mid_val

    #Define Figure Size
    size = size
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = size[0]
    fig_size[1] = size[1]
    plt.rcParams["figure.figsize"] = fig_size
    
    #Compute Ray Density
    if overlay != None:
        #Make Transparent Colorbar
        color1 = colors.colorConverter.to_rgba('white')
        color2 = colors.colorConverter.to_rgba('black')
        cmap1 = colors.LinearSegmentedColormap.from_list('my_cmap2',[color2,color1],256)
        cmap1._init()
        alphas = np.linspace(0.5, 0, cmap1.N+3) #Change degree of transparency
        cmap1._lut[:,-1] = alphas
        #Ray Density
        raydensity = ray_density_plot(raypaths, nodes)

    #Setup subplots
    fig = plt.figure() 
    ax3 = fig.add_subplot(111)

    #Third Plot Plot
    if type(bathy) != str:
        ax3.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=2)
    if type(src_loc) != str:
        ax3.scatter(src_loc[:,0], src_loc[:,1], s=60*fontscale, color="black", marker="v")
    dif = ax3.imshow(diff.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min_diff, vmax=elev_max_diff, 
                     cmap='seismic', norm=MidpointNormalize(midpoint=mid_val), aspect=aspect)
    ax3.set_title('Starting Model Subtraced from the Inverted Model',fontweight="bold", size=20.0*fontscale) # Title
    ax3.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    ax3.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax3.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax3.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)
    
    if overlay != None:
        #Ray-Density Image
        img2 = ax2.imshow(overlay, extent=[nx_min,nx_max,nz_min,nz_max], interpolation='nearest', 
                          cmap=cmap1, vmin=0, vmax=1)

    #Positioning Adjustments
    fig.subplots_adjust(right=0.8)
    fig.subplots_adjust(hspace=0.5)

    #Colorbar difference
    ax_cbar1 = fig.add_axes([0.7+shift[1], 0.3, 0.03, 0.40])
    cbar1 = fig.colorbar(dif, cax=ax_cbar1, orientation='vertical')
    cbar1.ax.tick_params(labelsize=16.0*fontscale)
    cbar1.set_label("Velocity Difference (m/s)", fontsize=18.0*fontscale)

    if letterlbl != None:
        #Add text
        fig.text(0.21+shift[0], 0.885, letterlbl[0], fontweight="bold", fontsize=18*fontscale)
        fig.text(0.21+shift[0], 0.602, letterlbl[1], fontweight="bold", fontsize=18*fontscale)
        fig.text(0.21+shift[0], 0.319, letterlbl[2], fontweight="bold", fontsize=18*fontscale)
        
    if savepath != None and savename != None:
        plt.gcf().savefig(savepath+savename+".png", bbox_inches = 'tight', pad_inches = 0.1)

    plt.show()
    
    
#######################
### PLOT TWO MODELS ###
#######################

def plot_truemdl_startmdl(truemdl, startmdl, nodes, spacing, savepath=None, savename=None, vrange="true", 
                          unit="m", fontscale=1, shift=[0,0], letterlbl=None, size=[72,36], cmap='jet',
                         bathy=None, src_loc=None):
    
    """ 
    This code is developed to gereate a high quality image of the starting model and true model for
    a synthetic inversion.
    """

    #########################
    ### READ BINARY FILES ###
    #########################

    #Read true model
    with open(truemdl, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    smesh_1 = np.reshape(data, [nodes[0], nodes[1]])

    #Read starting model
    with open(startmdl, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    smesh_2 = np.reshape(data, [nodes[0], nodes[1]])
    
    if type(bathy) == str:
        #Read the bathymetry file
        with open(bathy, 'rb') as f:
            bathy = np.fromfile(f, dtype=np.float32)
        if unit=="[km]":
            for i in range(len(bathy)):
                bathy_mesh = bathy/1000.0
        else:
            bathy_mesh = bathy
    else:
        bathy="None"
        
    #if src_loc == None:
    #    #If the source location is None set to string
    #    src_loc="None"

    ###################
    ### Plot Models ###
    ###################
    
    #Establish units, default is meters
    if unit == "[km]":
        spacing = spacing/1000.0
    
    #Define Model Dimensions
    nx_min = 0 * spacing
    nx_max = nodes[0] * spacing
    nz_min = nodes[1] * spacing
    nz_max = 0 * spacing

    #Define Figure Size
    size = [72,36]
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = size[0]
    fig_size[1] = size[1]
    plt.rcParams["figure.figsize"] = fig_size

    #Colorbar range, define limits based on starting or true model, or user specified
    if vrange == "true":
        elev_min, elev_max = min(smesh_1.flatten()), max(smesh_1.flatten())
    elif vrange == "starting":
        elev_min, elev_max = min(smesh_2.flatten()), max(smesh_2.flatten())
    else:
        elev_min, elev_max = vrange[0], vrange[1]

    #Setup subplots
    fig = plt.figure() 
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    #First Plot
    if type(bathy) != str:
        ax1.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=6)
    if type(src_loc) != str:
        ax1.scatter(src_loc[:,0], src_loc[:,1], s=500, color="black", marker="v")
    im = ax1.imshow(smesh_1.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min, vmax=elev_max, cmap=cmap,
                   aspect=2.0)
    #ax1.set_title('Eastern Mediterranean Model',fontweight="bold", size=20.0*fontscale) # Title
    ax1.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    ax1.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax1.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax1.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)

    #Second Plot
    if type(bathy) != str:
        ax2.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=6)
    if type(src_loc) != str:
        ax2.scatter(src_loc[:,0], src_loc[:,1], s=500, color="black", marker="v")
    im = ax2.imshow(smesh_2.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min, vmax=elev_max, cmap=cmap,
                   aspect=2.0)
    #ax2.set_title('Smoothed Eastern Mediterranean Model',fontweight="bold", size=20.0*fontscale) # Title
    ax2.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    ax2.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax2.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax2.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)

    #Colorbar
    ax_cbar = fig.add_axes([0.81+shift[1], 0.065, 0.03, 0.9])#0.140,,0.725])
    cbar = fig.colorbar(im, cax=ax_cbar, orientation='vertical')
    cbar.ax.tick_params(labelsize=16.0*fontscale)
    cbar.set_label("Velocity (m/s)", fontsize=18.0*fontscale)
    
    if letterlbl != None:
        #Add text
        fig.text(0.086+shift[0], 0.86, letterlbl[0], fontweight="bold", fontsize=18*fontscale)
        fig.text(0.086+shift[0], 0.45, letterlbl[1], fontweight="bold", fontsize=18*fontscale)
        
    if savepath != None and savename != None:
        plt.gcf().savefig(savepath+savename+".png", bbox_inches = 'tight', pad_inches = 0.1)
    
    plt.tight_layout(w_pad=1,h_pad=1)
    plt.show()
    
##########################################################
### PLOT TWO MODELS AND THE DIFFERENCE BETWEEN THE TWO ###
##########################################################

def plot_2mdls_plus_diff(model1, model2, nodes, spacing, savepath=None, savename=None, overlay=None, 
                         vrange_vel="mdl1", vrange_diff=None, mid_val=0, unit="m", fontscale=1, shift=[0,0], 
                         letterlbl=None, size=[72,48], aspect=1.0,
                         cmap='jet', bathy=None, src_loc=None, save_diff_mdl=False):

    """ 
    This code is intended to generate a high quality image of two models, and the difference between
    the two.
    """
    
    #########################
    ### READ BINARY FILES ###
    #########################

    #Open Model 1
    with open(model1, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    smesh_1 = np.reshape(data, [nodes[0], nodes[1]])

    #Open Model 2
    with open(model2, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    smesh_2 = np.reshape(data, [nodes[0], nodes[1]])
    
    #Read the bathymetry file
    with open(bathy, 'rb') as f:
        bathy = np.fromfile(f, dtype=np.float32)
    if unit=="[km]":
        for i in range(len(bathy)):
            bathy_mesh = bathy/1000.0
    else:
        bathy_mesh = bathy

    ###################
    ### Plot Models ###
    ###################

    #Establish units, default is meters
    if unit == "[km]":
        spacing = spacing/1000.0

    #Second model subtracted from first model as difference
    diff = smesh_2 - smesh_1
    
    if save_diff_mdl == True:
        #Manually specify a directory
        direct = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/006_RealDataTests/134_Nov26_EvenItCnstRampDW/_finalmodels/"
        a = np.array(diff,"float32")
        output_file = open(direct+'finalmdl_134_diffstart', 'wb')
        a.tofile(output_file)
        output_file.close()

    #Define Model Dimensions
    nx_min = 0 * spacing
    nx_max = nodes[0] * spacing
    nz_min = nodes[1] * spacing
    nz_max = 0 * spacing

    #Colorbar range for velocity models, define limits based on starting or true model, or user specified
    if vrange_vel == "model1":
        elev_min, elev_max = min(smesh_1.flatten()), max(smesh_1.flatten())
    elif vrange_vel == "model2":
        elev_min, elev_max = min(smesh_2.flatten()), max(smesh_2.flatten())
    else:
        elev_min, elev_max = vrange_vel[0], vrange_vel[1]

    #Colorbar range for difference plot, define limits based on min/max difference or user specified
    if vrange_diff == None:
        elev_min_diff, elev_max_diff = min(diff.flatten()), max(diff.flatten())
        mid_val=mid_val
    else:
        elev_min_diff, elev_max_diff = vrange_diff[0], vrange_diff[1]
        mid_val=mid_val

    #Define Figure Size
    size = size
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = size[0]
    fig_size[1] = size[1]
    plt.rcParams["figure.figsize"] = fig_size
    
    #Compute Ray Density
    if overlay != None:
        #Make Transparent Colorbar
        color1 = colors.colorConverter.to_rgba('white')
        color2 = colors.colorConverter.to_rgba('black')
        cmap1 = colors.LinearSegmentedColormap.from_list('my_cmap2',[color2,color1],256)
        cmap1._init()
        alphas = np.linspace(0.5, 0, cmap1.N+3) #Change degree of transparency
        cmap1._lut[:,-1] = alphas
        #Ray Density
        raydensity = ray_density_plot(raypaths, nodes)

    #Setup subplots
    fig = plt.figure() 
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)

    #First Plot
    if type(bathy) != str:
        ax1.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=2)
    if type(src_loc) != str:
        #Fudge with the plus 5.0 km, I've incorrectly cropped my models, FIX AFTER!
        ax1.scatter(src_loc[:,0], src_loc[:,1], s=60*fontscale, color="black", marker="v")
    im = ax1.imshow(smesh_1.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min, vmax=elev_max, cmap=cmap,
                   aspect=aspect)
    ax1.set_title('Starting Model',fontweight="bold", size=20.0*fontscale) # Title
    ax1.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    ax1.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax1.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax1.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)
    
    if overlay != None:
        #Ray-Density Image
        img2 = ax1.imshow(overlay, extent=[nx_min,nx_max,nz_min,nz_max], interpolation='nearest', 
                          cmap=cmap1, vmin=0, vmax=1)

    #Second Plot
    if type(bathy) != str:
        ax2.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=2)
    if type(src_loc) != str:
        ax2.scatter(src_loc[:,0], src_loc[:,1], s=60*fontscale, color="black", marker="v")
    im = ax2.imshow(smesh_2.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min, vmax=elev_max, cmap=cmap,
                   aspect=aspect)
    ax2.set_title('Inverted Model',fontweight="bold", size=20.0*fontscale) # Title
    ax2.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    ax2.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax2.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax2.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)
    
    if overlay != None:
        #Ray-Density Image
        img2 = ax2.imshow(overlay, extent=[nx_min,nx_max,nz_min,nz_max], interpolation='nearest', 
                          cmap=cmap1, vmin=0, vmax=1)

    #Third Plot Plot
    if type(bathy) != str:
        ax3.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=2)
    if type(src_loc) != str:
        ax3.scatter(src_loc[:,0], src_loc[:,1], s=60*fontscale, color="black", marker="v")
    dif = ax3.imshow(diff.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min_diff, vmax=elev_max_diff, 
                     cmap='seismic', norm=MidpointNormalize(midpoint=mid_val), aspect=aspect)
    ax3.set_title('Starting Model Subtraced from the Inverted Model',fontweight="bold", size=20.0*fontscale) # Title
    ax3.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    ax3.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax3.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax3.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)
    
    if overlay != None:
        #Ray-Density Image
        img2 = ax2.imshow(overlay, extent=[nx_min,nx_max,nz_min,nz_max], interpolation='nearest', 
                          cmap=cmap1, vmin=0, vmax=1)

    #Positioning Adjustments
    fig.subplots_adjust(right=0.8)
    fig.subplots_adjust(hspace=0.5)

    #Colorbar velocities
    ax_cbar = fig.add_axes([0.7+shift[1], 0.4075, 0.02, 0.4725])
    cbar = fig.colorbar(im, cax=ax_cbar, orientation='vertical')
    cbar.ax.tick_params(labelsize=16.0*fontscale)
    cbar.set_label("Velocity (m/s)", fontsize=18.0*fontscale)

    #Colorbar difference
    ax_cbar1 = fig.add_axes([0.7+shift[1], 0.1225, 0.02, 0.1925])
    cbar1 = fig.colorbar(dif, cax=ax_cbar1, orientation='vertical')
    cbar1.ax.tick_params(labelsize=16.0*fontscale)
    cbar1.set_label("Velocity Difference (m/s)", fontsize=18.0*fontscale)

    if letterlbl != None:
        #Add text
        fig.text(0.21+shift[0], 0.885, letterlbl[0], fontweight="bold", fontsize=18*fontscale)
        fig.text(0.21+shift[0], 0.602, letterlbl[1], fontweight="bold", fontsize=18*fontscale)
        fig.text(0.21+shift[0], 0.319, letterlbl[2], fontweight="bold", fontsize=18*fontscale)
        
    if savepath != None and savename != None:
        plt.gcf().savefig(savepath+savename+".png", bbox_inches = 'tight', pad_inches = 0.1)

    plt.show()
    
##################################################################
### PLOT MODEL, PLOT DIFFERENCE STARTING, PLOT DIFFERENE EXACT ###
##################################################################

def plot_mdl_plus_2diff(model, startmdl, truemdl, nodes, spacing, savepath=None, savename=None, 
                        vrange_vel=None, vrange_diff="diff1", mid_val=0, unit="m", fontscale=1, 
                        shift=[0,0], letterlbl=None, size=[72,48], cmap="jet", bathy=None, src_loc=None):

    """ 
    This code is intended to generate a high quality image of one model, and the two difference 
    models.
    """
    
    #########################
    ### READ BINARY FILES ###
    #########################

    #Open Model 1
    with open(model, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    smesh = np.reshape(data, [nodes[0], nodes[1]])

    #Open Model 2
    with open(startmdl, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    startmdl = np.reshape(data, [nodes[0], nodes[1]])
    
    #Open Model 2
    with open(truemdl, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    truemdl = np.reshape(data, [nodes[0], nodes[1]])

    #Read the bathymetry file
    with open(bathy, 'rb') as f:
        bathy = np.fromfile(f, dtype=np.float32)
    if unit=="[km]":
        for i in range(len(bathy)):
            bathy_mesh = bathy/1000.0
    else:
        bathy_mesh = bathy
    
    ###################
    ### Plot Models ###
    ###################

    #Establish units, default is meters
    if unit == "[km]":
        spacing = spacing/1000.
        
    #Starting model subtracted from inverted model as first difference
    diff_1 = smesh - startmdl
    #True model subtracted from inverted model as first difference
    diff_2 = smesh - truemdl
        
    #Define Model Dimensions
    nx_min = 0 * spacing
    nx_max = nodes[0] * spacing
    nz_min = nodes[1] * spacing
    nz_max = 0 * spacing

    #Colorbar range for difference plot, define limits based on starting or true model, or user specified
    if vrange_diff == "diff1":
        elev_min, elev_max = min(diff_1.flatten()), max(diff_1.flatten())
    elif vrange_diff == "diff1":
        elev_min, elev_max = min(diff_2.flatten()), max(diff_2.flatten())
    else:
        elev_min, elev_max = vrange_diff[0], vrange_diff[1]

    #Colorbar range for velocity model, define limits based on min/max difference or user specified
    if vrange_vel == None:
        elev_min_vel, elev_max_vel = min(smesh.flatten()), max(smesh.flatten())
        mid_val=mid_val
    else:
        elev_min_vel, elev_max_vel = vrange_vel[0], vrange_vel[1]
        mid_val=mid_val

    #Define Figure Size
    size = size
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = size[0]
    fig_size[1] = size[1]
    plt.rcParams["figure.figsize"] = fig_size

    #Setup subplots
    fig = plt.figure() 
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)

    #First Plot
    if type(bathy) != str:
        ax1.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=2)
    if type(src_loc) != str:
        ax1.scatter(src_loc[:,0], src_loc[:,1], s=60*fontscale, color="black", marker="v")
    im = ax1.imshow(smesh.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min_vel, vmax=elev_max_vel, 
                    cmap=cmap, aspect=2.0)
    #ax1.set_title('Inverted Model',fontweight="bold", size=20.0*fontscale) # Title
    ax1.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    ax1.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax1.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax1.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)

    #Second Plot
    if type(bathy) != str:
        ax2.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=2)
    if type(src_loc) != str:
        ax2.scatter(src_loc[:,0], src_loc[:,1], s=60*fontscale, color="black", marker="v")
    dif = ax2.imshow(diff_1.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min, vmax=elev_max, 
                     cmap='seismic', norm=MidpointNormalize(midpoint=mid_val), aspect=2.0)
    #ax2.set_title('Starting Model Subtracted from the Inverted Model',fontweight="bold", size=20.0*fontscale) # Title
    ax2.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    ax2.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax2.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax2.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)

    #Third Plot Plot
    if type(bathy) != str:
        ax3.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=2)
    if type(src_loc) != str:
        #Adding 5 km is a fudge for now, I messed up with the OBS locations
        ax3.scatter(src_loc[:,0], src_loc[:,1], s=60*fontscale, color="black", marker="v")
    dif = ax3.imshow(diff_2.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min, vmax=elev_max, 
                     cmap='seismic', norm=MidpointNormalize(midpoint=mid_val), aspect=2.0)
    #ax3.set_title('Detrended Starting Model Subtracted from the Inverted Model',fontweight="bold", 
                  #size=20.0*fontscale) # Title
    ax3.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    ax3.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax3.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax3.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)

    #Positioning Adjustments
    fig.subplots_adjust(right=0.8)
    fig.subplots_adjust(hspace=0.5)

    #Colorbar velocities
    ax_cbar = fig.add_axes([0.7+shift[1], 0.691, 0.02, 0.19])
    cbar = fig.colorbar(im, cax=ax_cbar, orientation='vertical')
    cbar.ax.tick_params(labelsize=16.0*fontscale)
    cbar.set_label("Velocity (m/s)", fontsize=18.0*fontscale)

    #Colorbar difference
    ax_cbar1 = fig.add_axes([0.7+shift[1], 0.1225, 0.02, 0.4725])
    cbar1 = fig.colorbar(dif, cax=ax_cbar1, orientation='vertical')
    cbar1.ax.tick_params(labelsize=16.0*fontscale)
    cbar1.set_label("Velocity Difference (m/s)", fontsize=18.0*fontscale)

    if letterlbl != None:
        #Add text
        fig.text(0.21+shift[0], 0.885, letterlbl[0], fontweight="bold", fontsize=18*fontscale)
        fig.text(0.21+shift[0], 0.602, letterlbl[1], fontweight="bold", fontsize=18*fontscale)
        fig.text(0.21+shift[0], 0.319, letterlbl[2], fontweight="bold", fontsize=18*fontscale)
        
    if savepath != None and savename != None:
        plt.gcf().savefig(savepath+savename+".png", bbox_inches = 'tight', pad_inches = 0.1)

    plt.show()
    
#######################################################
### PLOT TRUE MODEL, STARTING MODEL, INVERTED MODEL ###
#######################################################

def plot_true_start_inv_mdls(truemdl, startmdl, invmdl, spacing, savepath=None, savename=None,
                             vrange_vel="mdltrue", unit="m", fontscale=1, shift=[0,0], letterlbl=None, 
                             size=[72,48], cmap="jet", bathy=None, src_loc=None):

    """ savepath
    This code is intended to generate a publication grade image of a true model, starting model, 
    and inverted model.
    """

    #########################
    ### READ BINARY FILES ###
    #########################

    #Read True Model
    with open(truemdl, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    smesh_1 = np.reshape(data, [nodes[0], nodes[1]])

    #Read Starting Model
    with open(startmdl, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    smesh_2 = np.reshape(data, [nodes[0], nodes[1]])

    #Read Inverted Model
    with open(invmdl, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    smesh_3 = np.reshape(data, [nodes[0], nodes[1]])
    
    #Read the bathymetry file
    with open(bathy, 'rb') as f:
        bathy = np.fromfile(f, dtype=np.float32)
    if unit=="[km]":
        for i in range(len(bathy)):
            bathy_mesh = bathy/1000.0
    else:
        bathy_mesh = bathy

    ###################
    ### Plot Models ###
    ###################

    #Establish units, default is meters
    if unit == "[km]":
        spacing = spacing/1000. 

    #Define Model Dimensions
    nx_min = 0 * spacing
    nx_max = nodes[0] * spacing
    nz_min = nodes[1] * spacing
    nz_max = 0 * spacing

    #Define Figure Size
    size = size
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = size[0]
    fig_size[1] = size[1]
    plt.rcParams["figure.figsize"] = fig_size

    #Colorbar range for velocity models, define limits based on starting or true model, or user specified
    if vrange_vel == "mdltrue":
        elev_min, elev_max = min(smesh_1.flatten()), max(smesh_1.flatten())
    elif vrange_vel == "mdlstart":
        elev_min, elev_max = min(smesh_2.flatten()), max(smesh_2.flatten())
    elif vrange_vel == "mdlinv":
        elev_min, elev_max = min(smesh_3.flatten()), max(smesh_3.flatten())
    else:
        elev_min, elev_max = vrange_vel[0], vrange_vel[1]

    #Setup subplots
    fig = plt.figure() 
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)

    #First Plot
    if type(bathy) != str:
        ax1.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=2)
    if type(src_loc) != str:
        ax1.scatter(src_loc[:,0], src_loc[:,1], s=60*fontscale, color="black", marker="v")
    im = ax1.imshow(smesh_1.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min, vmax=elev_max, 
                    cmap='jet')
    ax1.set_title('True Model',fontweight="bold", size=20.0*fontscale) # Title
    ax1.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    #ax1.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax1.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax1.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)

    #Second Plot
    if type(bathy) != str:
        ax2.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=2)
    if type(src_loc) != str:
        ax2.scatter(src_loc[:,0], src_loc[:,1], s=60*fontscale, color="black", marker="v")
    im = ax2.imshow(smesh_2.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min, vmax=elev_max, 
                    cmap='jet')
    ax2.set_title('Starting Model; Sparse Dataset',fontweight="bold", size=20.0*fontscale) # Title
    ax2.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    ax2.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax2.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax2.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)

    #Third Plot Plot
    if type(bathy) != str:
        ax3.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=2)
    if type(src_loc) != str:
        ax3.scatter(src_loc[:,0], src_loc[:,1], s=60*fontscale, color="black", marker="v")
    im = ax3.imshow(smesh_2.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min, vmax=elev_max, 
                    cmap='jet')
    ax3.set_title('Inverted Model',fontweight="bold", size=20.0*fontscale) # Title
    ax3.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    ax3.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax3.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax3.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)

    #Positioning Adjustments
    fig.subplots_adjust(right=0.8)
    fig.subplots_adjust(hspace=0.5)

    #Colorbar velocities
    ax_cbar = fig.add_axes([0.81+shift[1], 0.1275, 0.02, 0.7490])
    cbar = fig.colorbar(im, cax=ax_cbar, orientation='vertical')
    cbar.ax.tick_params(labelsize=16.0*fontscale)
    cbar.set_label("Velocity (m/s)", fontsize=18.0*fontscale)

    if letterlbl != None:
        #Add text
        fig.text(0.21+shift[0], 0.885, letterlbl[0], fontweight="bold", fontsize=18)
        fig.text(0.21+shift[0], 0.602, letterlbl[1], fontweight="bold", fontsize=18)
        fig.text(0.21+shift[0], 0.319, letterlbl[2], fontweight="bold", fontsize=18)
        
    if savepath != None and savename != None:
        plt.gcf().savefig(savepath+savename+".png", bbox_inches = 'tight', pad_inches = 0.1)

    plt.show()
    
##################################################################
### PLOT MODEL, PLOT DIFFERENCE STARTING, PLOT DIFFERENE EXACT ###
##################################################################

def plot_mdl_plus_2diff_manu(model, startmdl, truemdl, nodes, spacing, savepath=None, savename=None, 
                        vrange_vel=None, vrange_diff="diff1", mid_val=0, unit="m", fontscale=1, 
                        shift=[0,0], letterlbl=None, size=[72,48], cmap="jet", bathy=None, src_loc=None):

    """ 
    This code is intended to generate a high quality image of one model, and the two difference 
    models.
    """
    
    #########################
    ### READ BINARY FILES ###
    #########################

    #Open Model 1
    with open(model, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    smesh = np.reshape(data, [nodes[0], nodes[1]])

    #Open Model 2
    with open(startmdl, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    startmdl = np.reshape(data, [nodes[0], nodes[1]])
    
    #Open Model 2
    with open(truemdl, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    truemdl = np.reshape(data, [nodes[0], nodes[1]])

    #Read the bathymetry file
    with open(bathy, 'rb') as f:
        bathy = np.fromfile(f, dtype=np.float32)
    if unit=="[km]":
        for i in range(len(bathy)):
            bathy_mesh = bathy/1000.0
    else:
        bathy_mesh = bathy
    
    ###################
    ### Plot Models ###
    ###################

    #Establish units, default is meters
    if unit == "[km]":
        spacing = spacing/1000.
        
    #Starting model subtracted from inverted model as first difference
    diff_1 = smesh - startmdl
    #True model subtracted from inverted model as first difference
    diff_2 = smesh - truemdl
        
    #Define Model Dimensions
    nx_min = 0 * spacing
    nx_max = nodes[0] * spacing
    nz_min = nodes[1] * spacing
    nz_max = 0 * spacing

    #Colorbar range for difference plot, define limits based on starting or true model, or user specified
    if vrange_diff == "diff1":
        elev_min, elev_max = min(diff_1.flatten()), max(diff_1.flatten())
    elif vrange_diff == "diff1":
        elev_min, elev_max = min(diff_2.flatten()), max(diff_2.flatten())
    else:
        elev_min, elev_max = vrange_diff[0], vrange_diff[1]

    #Colorbar range for velocity model, define limits based on min/max difference or user specified
    if vrange_vel == None:
        elev_min_vel, elev_max_vel = min(smesh.flatten()), max(smesh.flatten())
        mid_val=mid_val
    else:
        elev_min_vel, elev_max_vel = vrange_vel[0], vrange_vel[1]
        mid_val=mid_val

    #Define Figure Size
    size = size
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = size[0]
    fig_size[1] = size[1]
    plt.rcParams["figure.figsize"] = fig_size

    #Setup subplots
    fig = plt.figure() 
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    
    #White water column
    for i in range(len(smesh.T[0,:])):
        for j in range(len(smesh.T[:,0])):
            if j * spacing < bathy_mesh[i]:
                smesh.T[j,i] = np.nan            

    #First Plot
    if type(bathy) != str:
        ax1.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=2)
    if type(src_loc) != str:
        ax1.scatter(src_loc[:,0], src_loc[:,1], s=60*fontscale, color="black", marker="v")
    im = ax1.imshow(smesh.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min_vel, vmax=elev_max_vel, 
                    cmap=cmap, aspect=2.0)
    #ax1.set_title('Inverted Model',fontweight="bold", size=20.0*fontscale) # Title
    ax1.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    #ax1.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax1.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax1.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)

    #Second Plot
    if type(bathy) != str:
        ax2.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=2)
    if type(src_loc) != str:
        ax2.scatter(src_loc[:,0], src_loc[:,1], s=60*fontscale, color="black", marker="v")
    dif = ax2.imshow(diff_1.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min, vmax=elev_max, 
                     cmap='seismic', norm=MidpointNormalize(midpoint=mid_val), aspect=2.0)
    #ax2.set_title('Starting Model Subtracted from the Inverted Model',fontweight="bold", size=20.0*fontscale) # Title
    ax2.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    #ax2.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax2.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax2.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)

    #Third Plot Plot
    if type(bathy) != str:
        ax3.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=2)
    if type(src_loc) != str:
        #Adding 5 km is a fudge for now, I messed up with the OBS locations
        ax3.scatter(src_loc[:,0], src_loc[:,1], s=60*fontscale, color="black", marker="v")
    dif = ax3.imshow(diff_2.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min, vmax=elev_max, 
                     cmap='seismic', norm=MidpointNormalize(midpoint=mid_val), aspect=2.0)
    #ax3.set_title('Detrended Starting Model Subtracted from the Inverted Model',fontweight="bold", 
                  #size=20.0*fontscale) # Title
    ax3.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    ax3.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax3.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax3.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)

    #Positioning Adjustments
    fig.subplots_adjust(right=0.8)
    fig.subplots_adjust(hspace=0.5)

    #Colorbar velocities
    ax_cbar = fig.add_axes([0.9725, 0.692, 0.0225, 0.296]) #[0.8225, 0.692, 0.0225, 0.296] 30km crop #[0.9725, 0.692, 0.0225, 0.296] 20 km crop
    cbar = fig.colorbar(im, cax=ax_cbar, orientation='vertical')
    cbar.ax.tick_params(labelsize=16.0*fontscale)
    cbar.set_label("Velocity (m/s)", fontsize=18.0*fontscale)

    #Colorbar difference
    ax_cbar1 = fig.add_axes([0.9725, 0.04, 0.0225, 0.624]) #[0.8225, 0.04, 0.0225, 0.624] 30km crop # 20 km crop
    cbar1 = fig.colorbar(dif, cax=ax_cbar1, orientation='vertical')
    cbar1.ax.tick_params(labelsize=16.0*fontscale)
    cbar1.set_label("Velocity Difference (m/s)", fontsize=18.0*fontscale)

    if letterlbl != None:
        #Add text
        fig.text(0.21+shift[0], 0.885, letterlbl[0], fontweight="bold", fontsize=18*fontscale)
        fig.text(0.21+shift[0], 0.602, letterlbl[1], fontweight="bold", fontsize=18*fontscale)
        fig.text(0.21+shift[0], 0.319, letterlbl[2], fontweight="bold", fontsize=18*fontscale)
        
    if savepath != None and savename != None:
        plt.gcf().savefig(savepath+savename+".png", bbox_inches = 'tight', pad_inches = 0.1)
        
    plt.tight_layout(w_pad=1,h_pad=1)
    
    plt.show()
    
##################################################################
### PLOT MODEL, PLOT DIFFERENCE STARTING, PLOT DIFFERENE EXACT ###
##################################################################

def plot_mdl_plus_3diff(model, model2, startmdl, truemdl, nodes, spacing, 
                        vrange_vel=None, vrange_diff="diff1", vrange_diff2=None,
                        mid_val=0, unit="m", fontscale=1, 
                        shift=[0,0], size=[72,48], cmap="jet", bathy=None, src_loc=None):

    """ 
    This code is intended to generate a high quality image of one model, and the two difference 
    models.
    """
    
    #########################
    ### READ BINARY FILES ###
    #########################

    #Open Model 1
    with open(model, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    smesh = np.reshape(data, [nodes[0], nodes[1]])
    
    #Open Model 2
    with open(model2, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    smesh2 = np.reshape(data, [nodes[0], nodes[1]])

    #Open Model 3
    with open(startmdl, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    startmdl = np.reshape(data, [nodes[0], nodes[1]])
    
    #Open Model 4
    with open(truemdl, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    truemdl = np.reshape(data, [nodes[0], nodes[1]])

    #Read the bathymetry file
    with open(bathy, 'rb') as f:
        bathy = np.fromfile(f, dtype=np.float32)
    if unit=="[km]":
        for i in range(len(bathy)):
            bathy_mesh = bathy/1000.0
    else:
        bathy_mesh = bathy
    
    ###################
    ### Plot Models ###
    ###################

    #Establish units, default is meters
    if unit == "[km]":
        spacing = spacing/1000.
        
    #Starting model subtracted from inverted model as first difference
    diff_1 = smesh - startmdl
    #True model subtracted from inverted model as first difference
    diff_2 = smesh - truemdl
    #External model subtracted from inverted model as thrid difference
    diff_3 = smesh - smesh2
        
    #Define Model Dimensions
    nx_min = 0 * spacing
    nx_max = nodes[0] * spacing
    nz_min = nodes[1] * spacing
    nz_max = 0 * spacing

    #Colorbar range for difference plot, define limits based on starting or true model, or user specified
    if vrange_diff == "diff1":
        elev_min, elev_max = min(diff_1.flatten()), max(diff_1.flatten())
    elif vrange_diff == "diff1":
        elev_min, elev_max = min(diff_2.flatten()), max(diff_2.flatten())
    else:
        elev_min, elev_max = vrange_diff[0], vrange_diff[1]

    #Colorbar range for velocity model, define limits based on min/max difference or user specified
    if vrange_vel == None:
        elev_min_vel, elev_max_vel = min(smesh.flatten()), max(smesh.flatten())
        mid_val=mid_val
    else:
        elev_min_vel, elev_max_vel = vrange_vel[0], vrange_vel[1]
        mid_val=mid_val

    #Define Figure Size
    size = size
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = size[0]
    fig_size[1] = size[1]
    plt.rcParams["figure.figsize"] = fig_size

    #Setup subplots
    fig = plt.figure() 
    ax1 = fig.add_subplot(411)
    ax2 = fig.add_subplot(412)
    ax3 = fig.add_subplot(413)
    ax4 = fig.add_subplot(414)
    
    #White water column
    for i in range(len(smesh.T[0,:])):
        for j in range(len(smesh.T[:,0])):
            if j * spacing < bathy_mesh[i]:
                smesh.T[j,i] = np.nan            

    #First Plot
    if type(bathy) != str:
        ax1.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=2)
    if type(src_loc) != str:
        ax1.scatter(src_loc[:,0], src_loc[:,1], s=60*fontscale, color="black", marker="v")
    im = ax1.imshow(smesh.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min_vel, vmax=elev_max_vel, 
                    cmap=cmap, aspect=2.0)
    #ax1.set_title('Inverted Model',fontweight="bold", size=20.0*fontscale) # Title
    ax1.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    #ax1.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax1.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax1.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)

    #Second Plot
    if type(bathy) != str:
        ax2.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=2)
    if type(src_loc) != str:
        ax2.scatter(src_loc[:,0], src_loc[:,1], s=60*fontscale, color="black", marker="v")
    dif = ax2.imshow(diff_1.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min, vmax=elev_max, 
                     cmap='seismic', norm=MidpointNormalize(midpoint=mid_val), aspect=2.0)
    #ax2.set_title('Starting Model Subtracted from the Inverted Model',fontweight="bold", size=20.0*fontscale) # Title
    ax2.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    #ax2.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax2.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax2.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)

    #Third Plot Plot
    if type(bathy) != str:
        ax3.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=2)
    if type(src_loc) != str:
        ax3.scatter(src_loc[:,0], src_loc[:,1], s=60*fontscale, color="black", marker="v")
    dif = ax3.imshow(diff_2.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min, vmax=elev_max, 
                     cmap='seismic', norm=MidpointNormalize(midpoint=mid_val), aspect=2.0)
    #ax3.set_title('Detrended Starting Model Subtracted from the Inverted Model',fontweight="bold", 
                  #size=20.0*fontscale) # Title
    ax3.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    ax3.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax3.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax3.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)
    
    #Fourth Plot
    if type(bathy) != str:
        ax4.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=2)
    if type(src_loc) != str:
        ax4.scatter(src_loc[:,0], src_loc[:,1], s=60*fontscale, color="black", marker="v")
    dif_2 = ax4.imshow(diff_3.T, extent=[nx_min,nx_max,nz_min,nz_max], vmin=vrange_diff2[0], vmax=vrange_diff2[1], 
                     cmap='seismic', norm=MidpointNormalize(midpoint=mid_val), aspect=2.0)
    #ax3.set_title('Detrended Starting Model Subtracted from the Inverted Model',fontweight="bold", 
                  #size=20.0*fontscale) # Title
    ax4.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    ax4.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax4.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax4.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)

    #Positioning Adjustments
    fig.subplots_adjust(right=0.8)
    fig.subplots_adjust(hspace=0.5)

    #Colorbar velocities
    ax_cbar = fig.add_axes([0.66, 0.79, 0.0175, 0.1975]) 
    cbar = fig.colorbar(im, cax=ax_cbar, orientation='vertical')
    cbar.ax.tick_params(labelsize=16.0*fontscale)
    cbar.set_label("Velocity (m/s)", fontsize=18.0*fontscale)

    #Colorbar difference
    ax_cbar1 = fig.add_axes([0.66, 0.295, 0.0175, 0.445])
    cbar1 = fig.colorbar(dif, cax=ax_cbar1, orientation='vertical')
    cbar1.ax.tick_params(labelsize=16.0*fontscale)
    cbar1.set_label("Velocity Difference (m/s)", fontsize=18.0*fontscale)
    
    #Colorbar difference2
    ax_cbar2 = fig.add_axes([0.66, 0.045, 0.0175, 0.1975])
    cbar2 = fig.colorbar(dif_2, cax=ax_cbar2, orientation='vertical')
    cbar2.ax.tick_params(labelsize=16.0*fontscale)
    cbar2.set_label("Velocity Difference (m/s)", fontsize=18.0*fontscale)
        
    plt.tight_layout(w_pad=1,h_pad=1)
    
    plt.show()
    
#########################
### PLOT SINGLE MODEL ###
#########################

def plot_petrelmdl_manu(smesh, nodes, spacing, vrange="true", mid_val=0, unit="m", fontscale=1, size=[72,36], 
                        cmap='jet', bathy=None, src_loc=None, aspect=1.0):
    
    """ 
    This code is developed to read in a segy file that has been altered in Petrel.
    """

    ###################
    ### Plot Models ###
    ###################
    
    #Read the bathymetry file
    #with open(bathy, 'rb') as f:
    #    bathy = np.fromfile(f, dtype=np.float32)
    #if unit=="[km]":
    #    for i in range(len(bathy)):
    #        bathy_mesh = bathy/1000.0
    #else:
    #    bathy_mesh = bathy
    
    #Establish units, default is meters
    if unit == "[km]":
        spacing[0] = spacing[0]/1000.0
        spacing[1] = spacing[1]/1000.0
    
    #Define Model Dimensions
    nx_min = 0 * spacing[0] + 12.922
    nx_max = nodes[0] * spacing[0] + 12.922 
    nz_min = nodes[1] * spacing[1]
    nz_max = 0 * spacing[1]
    
    #Define Figure Size
    size = size
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = size[0]
    fig_size[1] = size[1]
    plt.rcParams["figure.figsize"] = fig_size

    #Colorbar range, define limits based on starting or true model, or user specified
    elev_min, elev_max = vrange[0], vrange[1]

    #Setup subplots
    fig = plt.figure() 
    ax1 = fig.add_subplot(111)
    
    #White water column (velocities)
    #for i in range(len(smesh[0,:])):
    #    for j in range(len(smesh[:,0])):
    #        if smesh[j,i] == 0:
    #            smesh[j,i] = np.nan    

    #Plot
    if type(bathy) != str:
        ax1.plot(np.arange(len(bathy_mesh))*spacing[0], bathy_mesh, color='black', linewidth=2)
    if type(src_loc) != str:
        ax1.scatter(src_loc[:,0], src_loc[:,1], s=100*fontscale, color="black", marker="v")
    im = ax1.imshow(smesh, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min, vmax=elev_max, cmap=cmap,
                   aspect=aspect)
    ax1.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    ax1.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax1.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax1.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)

    #Colorbar
    ax_cbar = fig.add_axes([0.912, 0.365, 0.025, 0.275])
    cbar = fig.colorbar(im, cax=ax_cbar, orientation='vertical')
    cbar.ax.tick_params(labelsize=16.0*fontscale)
    cbar.set_label("Velocity (m/s)", fontsize=18.0*fontscale)
    
    plt.show()
    
from scipy import interpolate

##################################################################
### PLOT MODEL, PLOT DIFFERENCE STARTING, PLOT DIFFERENE EXACT ###
##################################################################

def plot_profile_2diff_manu(profile1, profile2, profile3, nodes, spacing, decimate_z=1.0, crop_seismic=1.0,
                        alpha=0.5, vrange_vel=[0,9000], vrange_profile=[-10,10], vrange_diff=[-1000,1000],
                        mid_val=0, unit="m", fontscale=1, size=[72,48], cmap="jet", bathy=None):
  
    ##############################
    ### Decimate & Crop Models ###
    ##############################
    
    #Read the bathymetry file
    if bathy != None:
        with open(bathy, 'rb') as f:
            bathy = np.fromfile(f, dtype=np.float32)
        if unit=="[km]":
            for i in range(len(bathy)):
                bathy_mesh = bathy/1000.0
        else:
            bathy_mesh = bathy    
    #Decimate z
    profile_1 = np.zeros((int((len(profile1[:,0])-1)/decimate_z+1),len(profile1[0,:])))
    profile_2 = np.zeros((int((len(profile2[:,0])-1)/decimate_z+1),len(profile2[0,:])))
    profile_3 = np.zeros((int((len(profile3[:,0])-1)/decimate_z+1),len(profile3[0,:])))
    for i in range(len(profile1[:,0])):
        profile_1[i/2,:] = profile1[i,:]
        profile_2[i/2,:] = profile2[i,:]
        profile_3[i/2,:] = profile3[i,:]
    spacing[1] = spacing[1]*decimate_z
    nodes[1] = (nodes[1]-1.0)/decimate_z+1
    
    #Establish units, default is meters
    if unit == "[km]":
        spacing[0] = spacing[0]/1000.0
        spacing[1] = spacing[1]/1000.0
    
    #Define Model Dimensions
    nx_min = 0 * spacing[0] + 12.922
    nx_max = nodes[0] * spacing[0] + 12.922 
    nz_min = nodes[1] * spacing[1]
    nz_max = 0 * spacing[1]
    
    #Crop z for profile 1
    profile_1_crop = profile_1[0:int(crop_seismic*nodes[1]),:]
    nz_min_crop = crop_seismic*nz_min
        
    ###################
    ### Plot Models ###
    ###################

    #Define Figure Size
    size = size
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = size[0]
    fig_size[1] = size[1]
    plt.rcParams["figure.figsize"] = fig_size

    #Setup subplots
    fig = plt.figure() 
    ax1 = fig.add_subplot(411)
    ax2 = fig.add_subplot(412)
    ax3 = fig.add_subplot(413)
    ax4 = fig.add_subplot(414)
    
    #Interpolate the bathymetric profile TEMP
    def upscale_toy2dac_bathy(bathy, nx):

        #Upscale model
        xarange = np.arange(0,3397,1)
        xarange_update = np.arange(0, 3396, 1)
        
        #1D Interpolation for bathymetry
        tmp_bathy = interpolate.interp1d(xarange,bathy)
        bathy_update = tmp_bathy(xarange_update)
        
        #Output scaled to meteres for TOY2DAC - doing so prior results in a memory error
        return bathy_update
    
    bathy_upscale = upscale_toy2dac_bathy(bathy,len(bathy_mesh))
    
    print len(profile_2[:,0])
    print len(profile_2[0,:])
    print spacing[1]
    
    #White water column
    for i in range(len(profile_2[0,:])):
        for j in range(len(profile_2[:,0])):
            if j * spacing[1] < bathy_upscale[i]:
                profile_2[j,i] = np.nan
                
    print profile_2[:,0]
        
    #First Plot
    im = ax1.imshow(profile_1_crop, extent=[nx_min,nx_max,nz_min_crop,nz_max], vmin=vrange_profile[0], 
                    vmax=vrange_profile[1], cmap="gray", aspect=3.0/crop_seismic)
    #ax1.set_title('Inverted Model',fontweight="bold", size=20.0*fontscale) # Title
    ax1.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    #ax1.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax1.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax1.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)

    #Second Plot
    vel = ax2.imshow(profile_2, extent=[nx_min,nx_max,nz_min,nz_max], vmin=vrange_vel[0], vmax=vrange_vel[1], 
                     cmap=cmap, aspect=3.0)
    ax2.imshow(profile_1, extent=[nx_min,nx_max,nz_min,nz_max], vmin=vrange_profile[0], 
                    vmax=vrange_profile[1], cmap="gray", aspect=3.0, alpha=alpha)
    #ax2.set_title('Starting Model Subtracted from the Inverted Model',fontweight="bold", size=20.0*fontscale) # Title
    ax2.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    #ax2.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax2.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax2.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)

    #Third Plot Plot
    dif = ax3.imshow(profile_3, extent=[nx_min,nx_max,nz_min,nz_max], vmin=vrange_diff[0], vmax=vrange_diff[1], 
                     cmap='seismic', norm=MidpointNormalize(midpoint=mid_val), aspect=3.0)
    ax3.imshow(profile_1, extent=[nx_min,nx_max,nz_min,nz_max], vmin=vrange_profile[0], 
                    vmax=vrange_profile[1], cmap="gray", aspect=3.0, alpha=alpha)
    #ax3.set_title('Detrended Starting Model Subtracted from the Inverted Model',fontweight="bold", 
                  #size=20.0*fontscale) # Title
    ax3.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    ax3.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax3.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax3.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)
    
    #Fourth Plot
    dif = ax4.imshow(profile_3, extent=[nx_min,nx_max,nz_min,nz_max], vmin=vrange_diff[0], vmax=vrange_diff[1], 
                     cmap='seismic', norm=MidpointNormalize(midpoint=mid_val), aspect=3.0)
    #ax3.set_title('Detrended Starting Model Subtracted from the Inverted Model',fontweight="bold", 
                  #size=20.0*fontscale) # Title
    ax4.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    ax4.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax4.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax4.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)

    #Positioning Adjustments
    fig.subplots_adjust(right=0.8)
    fig.subplots_adjust(hspace=0.5)

    #Colorbar velocities
    ax_cbar = fig.add_axes([0.705, 0.78, 0.0175, 0.21]) 
    cbar = fig.colorbar(im, cax=ax_cbar, orientation='vertical')
    cbar.ax.tick_params(labelsize=16.0*fontscale)
    cbar.set_label("Velocity (m/s)", fontsize=18.0*fontscale)

    #Colorbar difference
    ax_cbar1 = fig.add_axes([0.705, 0.5315, 0.0175, 0.21])
    cbar1 = fig.colorbar(vel, cax=ax_cbar1, orientation='vertical')
    cbar1.ax.tick_params(labelsize=16.0*fontscale)
    cbar1.set_label("Velocity Difference (m/s)", fontsize=18.0*fontscale)
    
    #Colorbar difference2
    ax_cbar2 = fig.add_axes([0.705, 0.035, 0.0175, 0.4585])
    cbar2 = fig.colorbar(dif, cax=ax_cbar2, orientation='vertical')
    cbar2.ax.tick_params(labelsize=16.0*fontscale)
    cbar2.set_label("Velocity Difference (m/s)", fontsize=18.0*fontscale)
        
    plt.tight_layout(w_pad=1,h_pad=1)
    
    plt.show()
    
#############
### FILES ###
#############

#True Model (Marmousi)
#dirr = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/marmousi/run_marmousi_template_base/" 
#model = "vp_Marmousi_exact"
#mdltrue = dirr + model 

#Starting Model (Marmousi)
#dirr = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/marmousi/run_marmousi_template_base/" 
#model = "vp_Marmousi_init"
#mdlstart = dirr + model

#True Model (EMSynth)
#dirr = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/000_DATA/3_SynthEmed/4_SyntheticEMNov/model_50m/" 
#model = "vp"
#mdltrue = dirr + model 

#Starting Model (EMSynth)
#dirr = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/000_DATA/3_SynthEmed/4_SyntheticEMNov/model_50m/" 
#model = "vp_smooth"
#mdlstart = dirr + model

#Starting Model (EM Real)
dirr = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/000_DATA/2_RealData/17_110StartMdl_THEEND/0_models/model_50m/"
model = "vp_smooth_20crop"
mdlstart = dirr + model

#Detrended Model (EM Real)
dirr = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/000_DATA/2_RealData/17_110StartMdl_THEEND/0_models/model_50m/"
model = "vp_smooth_detrended_preupsc_20km"
mdldetrend = dirr + model

#Inverted Model One
dirr = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/009_Thesis_Inversions/2_EMReal/0_Manu_Results_Final_043/_finalmodels/" 
model = "model043_crop_20km"
invmdl1 = dirr + model

#Inverted Model Two
#dirr = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/000_DATA/2_RealData/17_110StartMdl_THEEND/0_models/model_50m/"
#model = "vp_smooth_detrended_preupsc_20km"
#detrendmdl = dirr + model

#Ray path file
#dirr = "/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/00_TOMO2DCONVERTED/2_fwdmdlraypaths/" 
#model = "Marmousi_FullAcqui"
#rays = dirr + model

#Petrel Model
segy_dir = "/home/chris/0_GlobeClaritas/0_Projects/Ref_EMed_Line1/COMMON/DATA/00_Reflection/3_savedSEGY/"
segy_file = "Velocity_043_CroppedAboveMultiple.sgy"
filename = segy_dir+segy_file
segyoffset = []
with segyio.open(filename) as f:
    profile = np.zeros((len(f.trace[0]),len(f.offsets)))
    obsoffset = f.offsets
    for o in range(len(obsoffset)):
        profile[:,o] = f.trace[len(obsoffset)-o-1]

#Bathymetry file (optional)
#bathy = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/000_DATA/2_RealData/17_110StartMdl_THEEND/0_models/model_50m/fbathy"
bathy = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/008_2020RealDataTests/043_MiscTestsFGProg1_SmthCnstL1/_finalmodels/fbathy_FirstLastOBSCrop"
#bathy = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/marmousi/run_marmousi_template_base/fbathy"

# Make a colormap from GMT file (optional)
cpt = gmtColormap('/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/rgb.VDL.NASA.cpt')#rgb.celw.Jan26.cpt')#
cpt_convert = colors.LinearSegmentedColormap('cpt', cpt)

#Import x,z OBS locations (optional)
dirr = "/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/00_TOMO2DCONVERTED/2_EmedSynthetic/0_FinalSourceGeoms/"
f = "Synthetic_SourceGeom_Final_True"
#dirr = "/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/00_TOMO2DCONVERTED/1_Marmousi/0_Starting_True/"
#f = "source_geom_full"
src_loc = np.genfromtxt(dirr+f)

###################
### DEFINITIONS ###
###################

#Specify number of nodes [x,z]
#nodes = [681,141] #marmousi model
#nodes = [5201,701] #kims model
#nodes = [3801,401] #20km z 
#nodes = [3801,501] #25km z 
#nodes = [3801,601] #synthetic no crop
#nodes = [1401, 601] #synthetic hecataeus crop
#nodes = [4751,501] #20 km z crop real
#nodes = [4751,751] #no crop real
#nodes = [801,141] #TOMO2D Model
nodes = [3396,15001] #Petrel Model

#Specify homogeneous node spacing in meters
#spacing = 50.0 #50 for kims & synthetic, 40 for real, 25 for marmousi
spacing = [50.005891,1.0]

#(OPTIONAL) Specifiy minimum and maximum velocities for velocity colorbar [min,max]
vrange = [1450.0,7250.0]#[1450, 7250.0]

#(OPTIONAL) Specifiy minimum and maximum velocities for difference colorbar [min,max,midval]
drange = [-1000.0,1000.0]

#(OPTIONAL) Specifiy minimum and maximum velocities for difference colorbar [min,max,midval]
drange2 = [-1000.0,1000.0]

#(OPTIONAL) Shift the colorbar shift[1], or the letters shift[0]
#shift=[0, -0.02] #Two Plot (marmousi)
#shift=[0.12,0.05] #Two Plot (EMSynth)
#shift=[-0.24,0.0] #EM Synth Hec Crop, single
#shift=[0.0,0.10] #Two Plot (EMSynthCROP)
#shift=[0.01,0.1] #Single Plot (Kim's EM Model)
#shift = [-0.09,0.09] #Three Plot (marmousi)
#shift = [-0.085,0.01] #Three Plot crop 25km
#shift = [-0.085,0.065] #Three Plot crop 20km
#shift = [0.02,-0.03] #Three Plot (EM), 2.0 V.E.

#(OPTIONAL) Letter Lables
lbl = ["(A)","(B)","(C)"]

#Path to save file, NONE if you don't want to save it. Both savepath and savename have to be defined.
savepath = None

#Name to save file, .png extension will be automatically used.
savename = None

################
### PLOTTING ###
################

#Make Ray Density Image (NEEDS WORK, TOO MEMORY INTENSIVE)
#raydensity = ray_density_plot(rays, nodes)

#Plot a single model
#plot_onemdl(mdlstart, nodes, spacing, vrange=vrange, savepath=savepath, savename=savename,
#                      unit="[km]", fontscale=1.5, shift=shift, size=[24,12], cmap=cpt_convert,
#                      bathy=bathy, src_loc=src_loc, Title=None)

#Read in peterl model 
#plot_petrelmdl_manu(profile, nodes, spacing, vrange=vrange, unit="[km]", fontscale=1.5, size=[36,18], 
#                        cmap=cpt_convert, bathy="bathy", src_loc="src_loc", aspect=3.0)

#Plot one difference model
#plot_onemdl_diff(mdlstart, invmdl1, nodes, spacing, savepath=savepath, savename=savename, 
#                         vrange_diff=drange, unit="[km]", fontscale=3.0, shift=shift, 
#                         letterlbl=None, size=[72,48], aspect=2.0,
#                         bathy=bathy, src_loc=src_loc)

#Plot the true model and starting model
#plot_truemdl_startmdl(mdltrue, mdlstart, nodes, spacing, vrange=vrange, savepath=savepath, savename=savename,
#                      unit="[km]", fontscale=4.0, shift=shift, letterlbl=None, size=[36,24], cmap=cpt_convert, #Pres 4.0 FS
#                      bathy=bathy, src_loc="none")

#Plot two difference models 
#plot_mdl_plus_2diff(invmdl1, mdlstart, mdltrue, nodes, spacing, savepath=savepath, savename=savename, 
#                    vrange_vel=vrange, vrange_diff=drange, unit="[km]", fontscale=30, shift=shift, 
#                    letterlbl=lbl)

#Plot two models, then the difference between the two of them (mdl2 - mdl1)
#plot_2mdls_plus_diff(mdlstart, invmdl1, nodes, spacing, savepath=savepath, savename=savename,
#                     vrange_vel=vrange, vrange_diff=drange,unit="[km]", fontscale=1.75,
#                     shift=shift,  size=[36,24], aspect=2.0,
#                     cmap=cpt_convert, bathy=bathy, src_loc=src_loc, save_diff_mdl=True)

"""SOMETHING FOR REAL DATA"""
#Plot one model plus two difference models 
#plot_mdl_plus_2diff(invmdl1, mdlstart, mdldetrend, nodes, spacing, savepath=savepath, savename=savename,
#                    vrange_vel=vrange, vrange_diff=drange, unit="[km]", fontscale=1.25, shift=shift, #Pres 1.75 FS
#                    letterlbl=lbl, size=[36,24], cmap=cpt_convert, bathy=bathy, src_loc=src_loc)

#Plot one model plus two difference models 
#plot_mdl_plus_2diff_manu(invmdl1, mdlstart, mdldetrend, nodes, spacing, savepath=savepath, savename=savename,
#                    vrange_vel=vrange, vrange_diff=drange, unit="[km]", fontscale=2.6, shift=shift, #OTHER PREV 1.5
#                    letterlbl=None, size=[36,24], cmap=cpt_convert, bathy=bathy, src_loc=src_loc)


#Plot one model plus three difference models 
#plot_mdl_plus_3diff(invmdl1, invmdl2, mdlstart, mdltrue, nodes, spacing,
#                    vrange_vel=vrange, vrange_diff=drange, vrange_diff2=drange2, unit="[km]", 
#                    fontscale=1.75, shift=shift, #OTHER PREV 1.5
#                    size=[48,24], cmap=cpt_convert, bathy=bathy, src_loc=src_loc)

#Plot True, Starting, and Inverted Models (Not yet Implimented)
#plot_true_start_inv_mdls(mdltrue, mdlstart, mdl, nodes, spacing, size, Marmousi_8thAcqui_NoSmth  
#                         unit="km", fontscale=1.2, shift=shift, letterlbl=lbl)

"""Final Interp Image for Manuscript File"""

#Plotting seismic profile, velocity model and detrended difference with seimsic overlays
plot_profile_2diff_manu(profile1, profile2, profile3, nodes, spacing, decimate_z=2.0, crop_seismic=2.0/3.0,
                        alpha=0.6, vrange_vel=[1450,7250], vrange_profile=[-1,1], 
                        vrange_diff=[-1000,1000], mid_val=0, unit="m", fontscale=1.2, size=[48,24], 
                        cmap=cpt_convert, bathy=bathy)

#Plot three models with velocity model in left column and difference model in right column
plt_3synth_mdl(model_1, model_2, model_3, truemdl, nodes, bathy, src_loc_full, src_loc_sparse, src_loc_usparse,
               spacing, vrange, drange, unit="[m]", fontscale=1.75, aspect=2.0,
               shift=[0,0], letterlbl=None, size=[36,24], cmap="jet")

"""Having issues setting this up as a function, just use this for now"""

### Program to read in segy profile for plotting in Python, not formalized at all

import segyio #Not standard with python https://github.com/equinor/segyio

#QC Figure Size
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 24
fig_size[1] = 12
plt.rcParams["figure.figsize"] = fig_size

"""Profile One"""

segy_dir = "/home/chris/0_GlobeClaritas/0_Projects/Ref_EMed_Line1/COMMON/DATA/00_Reflection/3_savedSEGY/"
segy_file = "Profile07_043_CropBelowMultiple.sgy"

filename = segy_dir+segy_file
segyoffset = []

with segyio.open(filename) as f:
    profile1 = np.zeros((len(f.trace[0]),len(f.offsets)))
    obsoffset = f.offsets
    for o in range(len(obsoffset)):
        profile1[:,o] = f.trace[len(obsoffset)-o-1]

#print np.shape(profile_1)
#plt.imshow(profile_1,aspect="auto",cmap="seismic", vmin=-10, vmax=10)
#plt.show()

"""Profile Two"""

segy_file = "Velocity_043_CroppedAboveMultiple.sgy"

filename = segy_dir+segy_file
segyoffset = []

with segyio.open(filename) as f:
    profile2 = np.zeros((len(f.trace[0]),len(f.offsets)))
    obsoffset = f.offsets
    for o in range(len(obsoffset)):
        profile2[:,o] = f.trace[len(obsoffset)-o-1]

#print np.shape(profile_2)
#plt.imshow(profile_2,aspect="auto",cmap="jet", vmin=1450, vmax=7000)
#plt.show()

"""Profile Three"""

segy_file = "DetrendedVelUpdate_043_CropAboveMultiple.sgy"

filename = segy_dir+segy_file
segyoffset = []

with segyio.open(filename) as f:
    profile3 = np.zeros((len(f.trace[0]),len(f.offsets)))
    obsoffset = f.offsets
    for o in range(len(obsoffset)):
        profile3[:,o] = f.trace[len(obsoffset)-o-1]

#print np.shape(profile_3)
#plt.imshow(profile_3,aspect="auto",cmap="seismic", vmin=-1000, vmax=1000)
#plt.show()

#X-length = 169,770 - spacing = 50.005891
#Z-length = 15,000 - spacing = 1.0