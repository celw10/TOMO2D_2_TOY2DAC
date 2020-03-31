#RUN ME IN JUPYTER NOTEBOOK

"""
This code is still quite ad hoc. but I'm trying to clean it up a bit from an earlier version in V1.0
Also, fixing a bug with the RMS velcity computation
"""


import utm #From https://github.com/Turbo87/utm
import numpy as np
import matplotlib.pyplot as plt

filename = "/Profile07_crsstack.sgy"

#CDP location file This is from UTM zone 36
f = "/home/chris/0_GlobeClaritas/0_Projects/Ref_EMed_Line1/COMMON/DATA/00_Reflection/Profile07_CDP_X_Y.txt" 
CDP_xy = np.genfromtxt(f, skip_header=1)
#Organized as, [CDP, UTMx, UTMy]
CDP_xy = CDP_xy[:,2:5]
maxfold = len(CDP_xy)

#Convert to latitude-longitude
Prof_lat_long = np.zeros((2,maxfold))
for i in range(maxfold):
    Prof_lat_long[0,i] = utm.to_latlon(CDP_xy[i,1], CDP_xy[i,2], 36, "U")[0]
    Prof_lat_long[1,i] = utm.to_latlon(CDP_xy[i,1], CDP_xy[i,2], 36, "U")[1]

#Import latitude-longitude OBS locations
f = "/home/chris/Documents/1_data/EasternMed/Chris_Williams/locations/obs.pqrsdz" 
obs = np.genfromtxt(f)
del_obs = [7,9,18,19,20] #OBS removed from inverte
OBS_xy = np.delete(obs,del_obs,axis=0)
OBS_xy = OBS_xy[:,2:4].T

plt.plot(OBS_xy[1], OBS_xy[0], "r", linewidth=3)
plt.plot(Prof_lat_long[0], Prof_lat_long[1], "b")
plt.show()

#Want to know the minimum and Maximum CDP covered by a velocity model
#We will sort by longitude
min_obs_x = min(OBS_xy[1,:])
max_obs_x = max(OBS_xy[1,:])

#Find Nearest
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

min_CDP = find_nearest(Prof_lat_long[0,:],min_obs_x)
max_CDP = find_nearest(Prof_lat_long[0,:],max_obs_x)

print "Minimum and maximum CDP for reflection profile onto refraction: ", min_CDP, max_CDP

#Define Figure Size
size = [12,12]
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = size[0]
fig_size[1] = size[1]
plt.rcParams["figure.figsize"] = fig_size

#Make a promax style velocity file that can be converted to a .nmo file in globe claritas.

from scipy import interpolate
import time

#############
### FLAGS ###
#############

#Output an interval (0), RMS velocity (1) file, or velocity difference model (2)
#Note: A velocity difference model will still require a interval velocity model to compute traveltimes
#Everything must be the same size
type_out = 2

############################
### SETUP VELOCITY MODEL ###
############################

#Define nodes for inported model x,z
nodes = [3397,601] #OBS CROP [3397,601]

#Define the model spacing (constant in x and z)
space = 50.0 #in meters

#Define a velocity difference model if that is to be output
if type_out == 2:
    model_diff = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/008_2020RealDataTests/043_MiscTestsFGProg1_SmthCnstL1/_finalmodels/Mdl043_FirstLasOBSCrop_DiffStart"
    with open(model_diff, "rb") as f:
        data = np.fromfile(f, dtype=np.float32)
    mdl_diff = np.reshape(data,nodes).T
    
#Import velocity model
#NOTE: This has to be cropped from the first and last OBS
model = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/008_2020RealDataTests/043_MiscTestsFGProg1_SmthCnstL1/_finalmodels/Mdl043_FirstLastOBSCrop"
with open(model, "rb") as f:
    data = np.fromfile(f, dtype=np.float32)
mdl_interval = np.reshape(data,nodes).T

###################
### MAKE MODELS ###
###################

print 
print "Computing interval traveltimes from input..."
print

#Make the iterval traveltime matrix
mdl_interval_time = np.zeros((nodes[1],nodes[0]))
for i in range(nodes[0]):
    for j in range(nodes[1]):
        mdl_interval_time[j,i] = space/mdl_interval[j,i]

#Make the traveltime matrix
mdl_time = np.zeros((nodes[1],nodes[0]))
for i in range(nodes[0]):
    for j in range(nodes[1]):
        mdl_time[j,i] = np.sum(mdl_interval_time[0:j,i])*2 #Two way traveltime

#Velocity file out
if type_out == 1:
    print 
    print "Computing RMS velocities from input..."
    print
    #Compute RMS layer velocities, lengthy...
    tt=time.time()
    mdl_rms = np.zeros((nodes[1],nodes[0]))
    for i in range(nodes[0]):
        for j in range(nodes[1]):
            intvelsquare_inttime = np.zeros(j+1)
            for k in range(j+1):
                intvelsquare_inttime[k] = (mdl_interval[j,i]**2)*mdl_interval_time[j,i]
            mdl_rms[j,i]=np.sqrt(np.sum(intvelsquare_inttime)/np.sum(mdl_interval_time[0:j+1,i]))
    #Define out
    vel_mdl_out = mdl_rms
    print
    print "RMS velocity model generation time: ", time.time()-tt
    print
    
elif type_out == 2:    
    print 
    print "Outputting a velocity update model (FWI model - starting model) from input..."
    print
    vel_mdl_out = mdl_diff 

else:
    print 
    print "Outputting an interval velocity file..."
    print
    #Define out
    vel_mdl_out = mdl_interval 
    
#QC Plots

#Define Figure Size
size = size
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 24
fig_size[1] = 6
plt.rcParams["figure.figsize"] = fig_size

#Plot the Inverval Traveltimes with the initial model structure 
plt.imshow(mdl_interval_time,"jet")
plt.title('Interval Time: Model structure')
plt.xlabel('Node x')
plt.ylabel('Node z')
cbar = plt.colorbar()
cbar.set_label("Interval Time [s]")
plt.show()

#Plot the Traveltimes with the initial model structure 
plt.imshow(mdl_time,"jet")
plt.title('Total Time: Model structure')
plt.xlabel('Node x')
plt.ylabel('Node z')
cbar = plt.colorbar()
cbar.set_label("Time [s]")
plt.show()

#Plot the Inverval Velocities with the initial model structure
plt.imshow(mdl_interval,"jet")
plt.title('Inverval Velocity: Model Structure')
plt.xlabel('Node x')
plt.ylabel('Node z')
cbar = plt.colorbar()
cbar.set_label("Inverval Velocity (m/s)")
plt.show()

#Plot the Output Velocity model with the initial model structure 
plt.imshow(vel_mdl_out,"jet")
plt.title('Velocity Model Out: Model structure')
plt.xlabel('Node x')
plt.ylabel('Node z')
cbar = plt.colorbar()
cbar.set_label("Time [s]")
plt.show()

#############
### FLAGS ###
#############

#Optional Gaussian smoothing for traveltime model, 0 = no, [smthx, smthz]=yes
smth_times = [10,10]

#Optional Gaussian smoothing for velocity (RMS or Interval) model, 0 = no, [smthx, smthz]=yes
smth_vels = 0 #[25,10]

##############################
### CONVERT NODES TO CDP's ###
##############################

#Convert nodes to CDP's
nx, nz = nodes[0], nodes[1]
num_cdpx = abs(max_CDP - min_CDP)
xarange = np.arange(0,nx,1)
zarange = np.arange(0,nz,1)
xarange_update = np.arange(0, nx, float(nx)/float(num_cdpx))

print 
print "Converting model nodes to CDP's based on defined geometry above..."
print "Total number of CDP's: ", num_cdpx
print

#Rectangular Bivariate Spline interpolation, kx- ky- degrees to the bivariate spline
#Velocity Model
interp_spline = interpolate.RectBivariateSpline(zarange, xarange, vel_mdl_out)
mesh = interp_spline(zarange, xarange_update)
#Traveltimes
interp_spline = interpolate.RectBivariateSpline(zarange, xarange, mdl_time)
times = interp_spline(zarange, xarange_update)

#########################################
### OPTIONAL SMOOTHING FOR VELOCITIES ###
#########################################

from scipy.ndimage import gaussian_filter

#Optional Model Smoothing, Velocities
if smth_vels != 0:
    print
    print "Smoothing in the x- and z-dimension defined as: ", smth_vels[0], " and ", smth_vels[1], " respectfully."
    print
    smoothing = [smth_vels[1], smth_vels[0]]
    mesh = gaussian_filter(mesh, smoothing, mode='reflect')

#Optional Model Smoothing, Traveltimes
if smth_times != 0:
    print
    print "Smoothing in the x- and z-dimension defined as: ", smth_times[0], " and ", smth_times[1], " respectfully."
    print
    smoothing = [smth_times[1], smth_times[0]]
    times = gaussian_filter(times, smoothing, mode='reflect')

################
### QC PLOTS ###
################

#Plot the Output Velocity model with the initial model structure 
plt.imshow(vel_mdl_out,"jet")
plt.title('Velocity Model Out: Model structure')
plt.xlabel('Node x')
plt.ylabel('Node z')
cbar = plt.colorbar()
cbar.set_label("Time [s]")
plt.show()

#Plot the Output Velocity model with the CDP model structure 
plt.imshow(mesh,"jet")
plt.title('Velocity Model Out: CDP structure')
plt.xlabel('CDP x')
plt.ylabel('Node z')
cbar = plt.colorbar()
cbar.set_label("Time [s]")
plt.show()

#Plot the Traveltimes with the initial model structure 
plt.imshow(mdl_time,"jet")
plt.title('Total Time: Model structure')
plt.xlabel('Node x')
plt.ylabel('Node z')
cbar = plt.colorbar()
cbar.set_label("Time [s]")
plt.show()

#Plot the Traveltimes with the DCP model structure 
plt.imshow(times,"jet")
plt.title('Total Time: CDP structure')
plt.xlabel('CDP x')
plt.ylabel('Node z')
cbar = plt.colorbar()
cbar.set_label("Time [s]")
plt.show()

###########################################
### STRUCTURE PROMAX RMS VELOCITY FILE ###
##########################################

# The first column will be CDP, the second column will be traveltimes, and the third will be RMS velocities
""" 
STRUCTURE PROMAX VEL FILE: 
Col 1 -> CDP
Col 2 -> Traveltimes ()
Col 3 -> RMS Velocities (m/s)
"""

print 
print "Constructing PROMAX Velocity File..."
print 

PROMAX_vel_file = np.zeros((nodes[1],num_cdpx,3))

#This is just storing everything in a 3D Matrix
for i in range(num_cdpx):
    for j in range(nodes[1]):
        #Assign CDP Information
        PROMAX_vel_file[j,i,0] = min_CDP - i #Make sure the minimum CDP is the leftmost portion of the model
        #Assign Velocity Information
        PROMAX_vel_file[j,i,2] = mesh[j,i]
        #Assign Traveltime Information
        PROMAX_vel_file[j,i,1] = times[j,i]*1000.0
        
######################
### FINAL QC PLOTS ###
######################

print
print "### FINAL QC ###"
print

#Plot the Output Velocity model with the initial model structure 
plt.imshow(PROMAX_vel_file[:,:,0],"jet")
plt.title('CDPs')
plt.xlabel('CDP x')
plt.ylabel('Node z')
cbar = plt.colorbar()
cbar.set_label("Time [s]")
plt.show()

#Plot the Traveltimes with the initial model structure 
plt.imshow(PROMAX_vel_file[:,:,1],"jet")
plt.title('Traveltime')
plt.xlabel('CDP x')
plt.ylabel('Node z')
cbar = plt.colorbar()
cbar.set_label("Time [ms]")
plt.show()

#Plot the Traveltimes with the DCP model structure 
plt.imshow(PROMAX_vel_file[:,:,2],"jet")
plt.title('Velocities')
plt.xlabel('CDP x')
plt.ylabel('Node z')
cbar = plt.colorbar()
cbar.set_label("Velocity [m/s]")
plt.show()

#Retain every 10th row, and retain every 12th CDP,
#I think I'm off by one CDP, I thought I had 119, but having 120 is close enough
#idx_lst = np.arange((max_CDP-max_CDP), (min_CDP-max_CDP+1),12)-1 
#idx_lst[0] = idx_lst[0]+1
#Velocities for TD conversions we want everything relatively smooth
#mdl = model[1::10,idx_lst,:]

#Retain every 5th row, and retain every 6th CDP (for RMS)
#What about 5th row and 3rd column for interval?
idx_lst = np.arange((max_CDP-max_CDP), (min_CDP-max_CDP+1),6)-1#CDP SPACING, 6 
idx_lst[0] = idx_lst[0]+1
mdl = PROMAX_vel_file[1::5,idx_lst,:]#Z-SPACING,5

#Reshape to length n with three columns for CDP, TT, RMS Velocity, This is PROMAX format, use the PROMAX to NMO conversion in Globe Clarias
mdl = mdl.reshape(((len(mdl[:,0,0])*len(mdl[0,:,0])), 3))

#NOTE the first column is actually 120 after cropping, last CDP is still 3515

#Try this
out = np.sort(mdl.view('f8,f8,f8'), order=['f0'], axis=0).view(np.float)

from astropy.io import ascii
ascii.write(out, "/home/chris/0_GlobeClaritas/0_Projects/Ref_EMed_Line1/COMMON/DATA/00_Reflection/2_output_models/2020_043/difference_starting/mdl043_10SmthTTs_DS.ascii")
#Write the output to an ascii file, I bring this ascii file into excell and format it in acending CDP, then I save it to a .csv file, before changing it to a PROMAX .vel file. Then this file is converted to a Claritas.nmo file using the PROMAX to NMO feature. 

