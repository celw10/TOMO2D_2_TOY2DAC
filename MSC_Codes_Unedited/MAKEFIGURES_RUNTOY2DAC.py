#RUN ME IN JUPYTER NOTEBOOK, SEPARATE OUT THE FUNCTIONS FROM THE MAIN CODE AT THE END

import numpy as np
import os
import imageio
import scipy.misc as spm
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.colors as colors
import colorsys

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

#Set the colormap to center at desired value
class make_midpoint(colors.Normalize):
    """
    Modified from: http://chris35wills.github.io/matplotlib_diverging_colorbar/
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)
    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
    
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

#########################
### PLOT SINGLE MODEL ###
#########################
nodes = [0,0]
def plot_onemdl(mdl, nodes, spacing, savepath=None, savename=None, vrange="true", 
                        unit="[km]", fontscale=2.0, shift=[0,0], letterlbl=None, size=[72,36], cmap='jet',
                        bathy=None, src_loc=None, zcrop=nodes[1],Title="Model"):
    
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
    #Optinal model cropping
    smesh_1 = smesh_1[:,0:zcrop].T
    
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
    nz_min = zcrop * spacing
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
        ax1.scatter(src_loc[:,0], src_loc[:,1], s=50*fontscale, color="black", marker="v")
    im = ax1.imshow(smesh_1, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min, vmax=elev_max, cmap=cmap,
                   aspect=2.0)
    ax1.set_title(Title,fontweight="bold", size=20.0*fontscale) # Title
    ax1.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    ax1.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax1.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax1.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)

    #Positioning Adjustments
    fig.subplots_adjust(right=0.8)
    fig.subplots_adjust(hspace=0.5)
    
    #Colorbar
    ax_cbar = fig.add_axes([0.7+shift[1], 0.3, 0.03, 0.40])
    cbar = fig.colorbar(im, cax=ax_cbar, orientation='vertical')
    cbar.ax.tick_params(labelsize=16.0*fontscale)
    cbar.set_label("Velocity (m/s)", fontsize=18.0*fontscale)
    
    #Add a number to keep figures in order
    if letterlbl < 10:
        idx = "00"+str(letterlbl)+"_"
    elif letterlbl < 100:
        idx = "0"+str(letterlbl)+"_"
    else:
        idx = str(letterlbl)+"_"
        
    #Optionally add the inversion group number
    if letterlbl != None:
        fig.text(0.086+shift[0], 0.74, letterlbl, fontweight="bold", fontsize=18*fontscale)
    
    #Save the figure
    if savepath != None and savename != None:
        plt.gcf().savefig(savepath+idx+savename+".png", bbox_inches = 'tight', pad_inches = 0.1)

    #Delete current figure
    fig.clear(); del fig
    plt.close()
    
##########################################################
### PLOT TWO MODELS AND THE DIFFERENCE BETWEEN THE TWO ###
##########################################################
nodes=[0,0]
def plot_onemdl_diff(model1, model2, nodes, spacing, savepath=None, savename=None, overlay=None, 
                         vrange_diff=None, mid_val=0, unit="m", fontscale=1, shift=[0,0], 
                         letterlbl=None, size=[72,48], aspect=1.0,
                         cmap='jet', bathy=None, src_loc=None, title="Difference Model",
                         zcrop=nodes[1]):

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
    #Optinal model cropping
    smesh_1 = smesh_1[:,0:zcrop].T
    
    #Open Model 2
    with open(model2, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32)
    smesh_2 = np.reshape(data, [nodes[0], nodes[1]])
    #Optinal model cropping
    smesh_2 = smesh_2[:,0:zcrop].T
    
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
    nz_min = zcrop * spacing
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
        ax3.scatter(src_loc[:,0], src_loc[:,1], s=50*fontscale, color="black", marker="v")
    dif = ax3.imshow(diff, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min_diff, vmax=elev_max_diff, 
                     cmap='seismic', norm=MidpointNormalize(midpoint=mid_val), aspect=aspect)
    ax3.set_title(title,fontweight="bold", size=20.0*fontscale) # Title
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
    
    #Add a number to keep figures in order
    if letterlbl < 10:
        idx = "00"+str(letterlbl)+"_"
    elif letterlbl < 100:
        idx = "0"+str(letterlbl)+"_"
    else:
        idx = str(letterlbl)+"_"
        
    #Optionally add the inversion group number
    if letterlbl != None:
        fig.text(0.086+shift[0], 0.74, letterlbl, fontweight="bold", fontsize=18*fontscale)
    
    #Save the figure
    if savepath != None and savename != None:
        plt.gcf().savefig(savepath+idx+savename+".png", bbox_inches = 'tight', pad_inches = 0.1)

    #Delete current figure
    fig.clear(); del fig
    plt.close()
    
#! /usr/bin/env python

################################################################
# Version 2.0 
# Christopher Williams, Memorial University of Newfoundland
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# What is my purpose? Why am I here???
# -> MAKEFIGURES_V*** is intented to interface with 
#    RUNTOY2DAC_V***, and in particular the INVNAMES.out file
# -> MAKEFIGURES_V*** produces multiple high quality images of  
#    the desired inversion output from RUNTOY2DAC_V***.
# -> As of now MAKEFIGURES_V*** cannot be run on TORNGAT.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################################################################

####################
### DEFINE FLAGS ###
####################

#Plot the inverted model, 0=off, 1=on
plt_inv_mdl = [1,0]

#Plot gradient model, 0=off, 1=on
plt_grad = [0,0]

#Plot the difference between inverted_model-starting_model, 0=off, 1=on
diff_inv_start = [1,0]

#Plot the difference between inverted_model-true_model, 0=off, 1=on
diff_inv_true = [1,0]

#Make extra directories, inversion stats, gifs, forward models, and final models. Adds consistency, 0 off 1 on
extra_dirr = 1

##########################
### DEFINE DIRECTORIES ###
##########################

#Where to look for output (unaltered from RUNTOY2DAC_V***)
relative_path = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/009_Thesis_Inversions/1_EMSynthetic/3_Appendix_ManuParamsPreNLCG/"

#####################
### DEFINE MODELS ###
#####################

#Define starting model (EM)
#if diff_inv_start[0] == 1:
#    file_start = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/000_DATA/2_RealData/17_110StartMdl_THEEND/0_models/model_50m/vp_smooth"

#Define starting model (Marmousi)
#if diff_inv_start[0] == 1:
#    file_start = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/marmousi/run_marmousi_template/vp_Marmousi_init"

#Define starting model (EM Synthetic)
if diff_inv_start[0] == 1:
    file_start = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/000_DATA/3_SynthEmed/4_SyntheticEMNov/model_50m/vp_smooth"
    
#Define true model (Marmousi)
#if diff_inv_true[0] == 1:
#    file_true = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/marmousi/run_marmousi_template/vp_Marmousi_exact"

#Define true model (EM Synthetic)
if diff_inv_true[0] == 1:
    file_true = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/000_DATA/3_SynthEmed/4_SyntheticEMNov/model_50m/vp"

#####################
### DEFINE EXTRAS ###
#####################

#Bathymetry file (optional)
bathy = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/000_DATA/2_RealData/17_110StartMdl_THEEND/0_models/model_50m/fbathy"
#bathy = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/marmousi/run_marmousi_template_base/fbathy"
# Make a colormap from GMT file (optional)
#cpt = gmtColormap('/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/rgb.celw.Jan26.cpt')#rgb.kw.Nov09.cpt')
#cpt_convert = colors.LinearSegmentedColormap('cpt', cpt)

#Import x,z OBS locations Eastern Mediterranean
dirr = "/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/00_TOMO2DCONVERTED/2_EmedSynthetic/0_FinalSourceGeoms/"
f = "Synthetic_SourceGeom_Final_True"
src_loc = np.genfromtxt(dirr+f)

#Import x,z OBS locations Marmousi
#dirr = "/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/00_TOMO2DCONVERTED/1_Marmousi/0_Starting_True/"
#f = "source_geom_full"
#src_loc = np.genfromtxt(dirr+f)

##################################
### DEFINE PLOTTING PARAMETERS ###
##################################

#Specify number of nodes [x,z] 
#nodes = [4751,751] #(Real model 40m)
#nodes = [681,141] #(Marmousi)
nodes = [3801,601] #(EM Synthetic 50m)
#nodes = [1401,601] #(EM Synthetic, HR Crop, 50m)

#Specify homogeneous node spacing in meters
spacing = 50.00 #40.00, 25.00

#Specify maximum figure size [x,z]
size = [24,16]#[24,12], em, 16 for marmousi

#Specifiy minimum and maximum velocities for velocity colorbar [min,max]
vrange = [1500.0, 8000.0]

#Specifiy minimum and maximum velocities for difference colorbar [min,max,midval]
drange = [-800.0,800.0,0]

#Shift (laterally), the colorbar or letters
shift = [0.01, 0.11]

#Scaler for charachter sizes in images
fontscale=1.5

####################
### MAKE FIGURES ###
####################

#Read INVNAMES.out
f = open(relative_path+"INVNAMES.out","r")
model_names = f.read().splitlines()
f.close()

#Make a directory to store velocity model plots in
if plt_inv_mdl[0] == 1:
    os.mkdir(relative_path+"3_mdl_figs/")
    
#Make a directory to store the gradient plots in
if plt_grad[0] == 1:
    os.mkdir(relative_path+"4_gradient_figs/")
    
#Make a directory to store starting model difference plots in
if diff_inv_start[0] == 1:
    #Read the inverted model
    with open(file_start, "rb") as f:
        data = np.fromfile(f, dtype=np.float32)
    mdl_start = np.reshape(data,nodes)
    #Make the directory
    os.mkdir(relative_path+"5_mdl_start_diff_figs/")
    
#Make a directory to store starting model difference plots in
if diff_inv_true[0] == 1:
    #Read the inverted model
    with open(file_true, "rb") as f:
        data = np.fromfile(f, dtype=np.float32)
    mdl_true = np.reshape(data,nodes)
    #Make the directory
    os.mkdir(relative_path+"6_mdl_true_diff_figs/")
    
#Make Extra Directories (keeps everything clean)
if extra_dirr == 1:
    os.mkdir(relative_path+"7_invstats")
    os.mkdir(relative_path+"8_fwdmdl")
    os.mkdir(relative_path+"9_gifs")
    os.mkdir(relative_path+"_finalmodels")
        
#Iterate through every model and make figures for each
for idx in range(len(model_names)):
    
    #Inverted model location
    file_inverted = relative_path+"0_inverted_models/"+model_names[idx]
    
    if plt_grad[0] == 1:
        #Gradient model location
        file_gradient = relative_path+"2_gradient_optput/"+model_names[idx]
    
    #Plot inverted model
    if plt_inv_mdl[0] == 1:
        plot_onemdl(file_inverted, nodes, spacing, vrange=vrange, savepath=relative_path+"3_mdl_figs/", 
                    savename=str(model_names[idx]), unit="[km]", fontscale=fontscale, shift=shift, letterlbl=idx,
                    size=size, cmap="jet", bathy=bathy, src_loc=src_loc, zcrop=601,
                    Title="Velocity: Update "+str(idx+1))

    #Plot gradient
    if plt_grad[0] == 1:
        plot_onemdl(file_gradient, nodes, spacing, vrange=[-1,1], savepath=relative_path+"4_gradient_figs/", 
                    savename=str(model_names[idx]), unit="[km]", fontscale=fontscale, shift=shift, letterlbl=idx,
                    size=size, cmap="bwr", bathy=bathy, src_loc=src_loc, zcrop=601, 
                    Title="Gradient: Update "+str(idx+1))
    
    #Plot difference from inverted-starting
    if diff_inv_start[0] == 1:
        plot_onemdl_diff(file_start, file_inverted, nodes, spacing, savepath=relative_path+"5_mdl_start_diff_figs/", 
                         savename=str(model_names[idx]), vrange_diff=drange, unit="[km]", fontscale=fontscale, 
                         shift=shift, letterlbl=idx, size=size, aspect=2.0, bathy=bathy, src_loc=src_loc,
                         title="Starting Model Subtraced from the Inverted Model",
                         zcrop=601)
    
    #Plot difference from inverted-true
    if diff_inv_true[0] == 1:
        plot_onemdl_diff(file_true, file_inverted, nodes, spacing, savepath=relative_path+"6_mdl_true_diff_figs/", 
                         savename=str(model_names[idx]), vrange_diff=drange, unit="[km]", fontscale=fontscale, 
                         shift=shift, letterlbl=idx, size=size, aspect=2.0, bathy=bathy, src_loc=src_loc,
                         title="Starting Model Subtraced from the True Model",
                         zcrop=601)