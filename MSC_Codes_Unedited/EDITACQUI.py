#RUN ME IN JUPYTER NOTEBOOK, SEPARATE OUT THE FUNCTIONS FROM THE MAIN CODE AT THE END

### MAIN FUNCTION TO READ THE tx.in FB PICKS

def read_fbpicks(picks, shotloc, phase, unc, tt_window, SR):
    """
    |_________________________________________________________________________________
    | PURPOSE:
    |-> Function to read and reorder a tx.in first-break pick file to a phython array.
    |-> Optional manipulation of the picks in tx.in file.
    |_________________________________________________________________________________
    | INPUT:
    |-> picks: tx.in file preferably read in using numpy genfromtxt
    |-> shotloc: a file of length=nshots, with the x-dimension model positions of the shots
    |-> phase: a list of int that represent different desired phases to be retained in this analysis. 
    |   Refracted and ground waves are reccomended. 
    |-> unc: maximum uncertainty to be taken from tx.in file
    |-> tt_window: minimum and maximum traveltimes to read as tt_window=[mintt, maxtt]
    |-> SR: sample rate from segy data
    |_________________________________________________________________________________
    """
    
    #Find breaks in tx.in file, get rid of intermediate points
    for i in range(0, len(picks)-len(shotloc), 1): 
        if picks[i,1] == 1.0 and picks[i,3] == 0:
            picks = np.delete(picks,i,0)
    
    total = len(picks)
    
    print "Original length of tx.in file: ", total
    
    #Remove undesired phases
    picks, phasedeleted = desired_phase(picks, phase)
    
    #Remove high pick uncertainties
    picks, uncdeleted = desired_unc(picks, unc)

    #Retain picks in desired traveltime window
    picks, offsetdeleted = desired_ttwindow(picks, tt_window)
    
    print "Remaining first-break picks in tx.in file: ", len(picks) - len(shotloc) - 1
    print
    
    #Ensure everything adds up
    assert ((phasedeleted)+(uncdeleted)+(offsetdeleted)+(len(picks))) == total
    
    #Return the index for each shot, and number of picks for each shot
    pick_numbers, indices_start = num_picks(picks, shotloc)
    
    #Make a organized list with all data
    OBS = [[] for _ in range(len(shotloc))]
    for i in np.arange(len(shotloc)):
        if i == len(shotloc) - 1:
            OBS[i] = picks[indices_start[i]:len(picks)]
        else:
            OBS[i] = picks[indices_start[i]:indices_start[i+1]]

    #Change traveltime to number of time samples, discard pick type and uncertainty
    OBS_sg = []
    for i in np.arange(len(shotloc)):
        OBS_sg.append(OBS[i][:,0:2])
        OBS_sg[i][:,1] = OBS_sg[i][:,1]/SR
            
    return OBS_sg, pick_numbers

### SUPPORTING FUNCTIONS FOR...
### read_fbpicks

def desired_phase(picks, phase):
    """Remove phases from tx.in first-break pick file."""
    
    i, totalphase, tot = 0, 0, len(picks)
    while i < tot: 
        index = tot - i - 1
        if picks[index][3] in phase:
            totalphase = totalphase+1
        else:
            picks = np.delete(picks, index, 0)
        i = i + 1
        
    phasedeleted = tot - totalphase
        
    print "Undesirable phase picks removed: ", phasedeleted
    
    return picks, phasedeleted

def desired_unc(picks, unc):
    """Remove uncertainties from tx.in first-break pick file."""
    
    i, totalunc, tot = 0, 0, len(picks)
    while i < tot: 
        index = tot - i - 1
        if picks[index][2] <= unc: 
            totalunc = totalunc + 1
        else:
            picks = np.delete(picks, index, 0)
        i = i + 1
    
    uncdeleted = tot - totalunc
    
    print "Undesirable uncertainties removed: ", uncdeleted
    
    return picks, uncdeleted
    
def desired_ttwindow(picks, tt_window):
    """Remove picks outside a desired traveltime window."""

    i, totaloffset, tot = 0, 0, len(picks)
    while i < tot: 
        index = tot - i - 1
        if tt_window[0] <= picks[index][1] <= tt_window[1]: 
            totaloffset = totaloffset + 1
        else:
            if picks[index][1] < 0:
                totaloffset = totaloffset + 1
            else:
                picks = np.delete(picks, index, 0)
        i = i + 1
    
    offsetdeleted = tot - totaloffset
    
    print "Undesirable seismic arrivals removed: ", offsetdeleted
    
    return picks, offsetdeleted

def num_picks(picks, shotloc, breakidx=-1.000):
    """Save the number of first-break picks per shot"""
    
    indices_start, pick_numbers = [], []

    #Find the index corresponding to the start of the shot
    for i in range(0, len(picks), 1): 
        if picks[i,0] in shotloc and picks[i,3] == 0 and picks[i,1] == breakidx: 
            indices_start.append(i)

    #Append the total number of picks
    for i in range(0, len(indices_start), 1): 
        if i < (len(indices_start)-1):
            x = indices_start[i+1] - indices_start[i]
            pick_numbers.append(x)
        else:
            x = len(picks[:,0]) - indices_start[i]
            pick_numbers.append(x)
            
    return pick_numbers, indices_start

#READ IN GATHERS AND OFFSETS FROM SEGY FILES

def extract_segy(fileloc, numobs, maxfold, nsamples):
    """
    |_________________________________________________________________________________
    | PURPOSE:
    |-> Obtain offset information directly from SEGY headers.
    |-> Organize gathers into python array trace by trace.
    |_________________________________________________________________________________
    | IMPORTANT:
    |-> The seismic data must all be in the same folder and properly named.
    |-> The only thing that may change in the name of each file is a number,
    |-> This number must start a 1 and continue to the total number of OBS.
    |-> A consistent lable prior to and after the number is allowed.
    |_________________________________________________________________________________
    """
    
    maxfold_qc = 0
    obsoffset = [[] for _ in range(numobs)]
    gathers = np.zeros((numobs,maxfold,int(nsamples)))
    
    for i in range(0,numobs,1):
        data = fileloc[1] + str(i+1) + fileloc[2]
        filename = fileloc[0] + data
        with segyio.open(filename) as f:
            #Save all offsets
            obsoffset[i] = f.offsets/1000.0 #Want this in meters
            #QC the input maximum fold
            if len(obsoffset[i]) > maxfold_qc:
                maxfold_qc = len(obsoffset[i])
            #Save each trace over all offsets
            for o in range(len(obsoffset[i])):
                gathers[i,o,:] = f.trace[o]
    
    #Make sure the user input maxfold correctly
    assert maxfold == maxfold_qc
    #Makesure the user input nsamples correctly
    assert nsamples == len(gathers[i,o,:])
    
    return gathers, obsoffset

### FUNCTIONS TO ADJUST DATASET AND CROP MODEL

def max_offset(OBS, offset_max):
    """Delete data above the maximum desirable offset
    We want to delete everything, but we need to take note of how much we've deleted to the LEFT (of the OBS location)
    of the model. At this point these points are negative
    I've added a minimum offset, so the nearest data with complicated waveforms due to complicated bathymetry
    can be removed as well."""

    print "WARNING: Maximum offset has been restricted to: ", offset_max
    
    #Use max/min crop here
    indicies = [[] for _ in range(len(OBS))]
    indicies_all = [[] for _ in range(len(OBS))]
    for i in np.arange(len(OBS)):
        for j in np.arange(len(OBS[i][:,0])):
            if abs(OBS[i][j,0]) >= offset_max:
                indicies[i].append(j)
    
    #Information to shift the first break data
    indicies_shift = []
    for i in range(len(OBS)):
        tmp = 0
        for j in indicies[i]:
            if OBS[i][j,0] <= 0:
                tmp += 1
        indicies_shift.append(tmp)

    #Remove undesired offsets
    total_all=0
    for i in np.arange(len(OBS)):
        OBS[i] = np.delete(OBS[i],indicies[i],0)
        total_all += len(indicies[i])
    
    print "Traces removed, max offsets: ", total_all
    print
    
    return OBS, indicies_shift

def shift_data(OBS, shotloc, bulk_shift):
    """Shift entire dataset to have a origin of zero"""
    
    if bulk_shift != 0:
        print "WARNING: Bulk shift of ", bulk_shift, " applied to the acquisition file"
        print
        
    #Shift w.r.t. OBS position on model
    for i in range(len(OBS)):
        OBS[i][:,0] = OBS[i][:,0] + shotloc[i] + bulk_shift
        
    return OBS

def crop_model(OBS, min_max):
    """Specify the minimum and maximum extent of the model in the x dimension.
    To be applied after the data has been shifted to a consistent origin in shift_data above.
    min_max input as min_max = [xminimum, xmaximum]."""
    
    print "WARNING: maximum xcoordinate has been restricted to ", min_max[1]

    indicies_max = [[] for _ in range(len(OBS))]
    for i in np.arange(len(OBS)):
        for j in np.arange(len(OBS[i][:,0])):
            if OBS[i][j,0] >= min_max[1]:
                indicies_max[i].append(j)

    total=0
    for i in np.arange(len(OBS)):
        OBS[i] = np.delete(OBS[i],indicies_max[i],0)
        total += len(indicies_max[i])

    print "Traces removed: ", total
    print
        
    print "WARNING: minimum xcoordinate has been restricted to ", min_max[0]

    #Save indicies of all locations greater than the maximum offset
    indicies_min = [[] for _ in range(len(OBS))]
    for i in np.arange(len(OBS)):
        for j in np.arange(len(OBS[i][:,0])):
            if OBS[i][j,0] <= min_max[0]:
                indicies_min[i].append(j)

    total=0
    for i in np.arange(len(OBS)):
        OBS[i] = np.delete(OBS[i],indicies_min[i],0)
        total += len(indicies_min[i])

    print "Traces removed: ", total
    print
        
    return OBS, indicies_min

### Now we need to read in the data from SGY files

def indices_from_segy(OBS, offset):
    """
    |_________________________________________________________________________________
    | PURPOSE:
    |-> Obtain min/max indicies for segy data that relate to the min/max first break pick.
    |_________________________________________________________________________________
    | IMPORTANT:
    |-> The segy data cannot be shifted, offsets are with respect to the OBS location at this time.
    |-> If using plotsec, make sure that the shift value is zero when producing tx.in file.
    |_________________________________________________________________________________
    """
    
    #Setup lists
    minoffsets_ps = [[] for _ in range(len(OBS))]
    maxoffsets_ps = [[] for _ in range(len(OBS))]
    minoffsets_hdr = [[] for _ in range(len(OBS))]
    maxoffsets_hdr = [[] for _ in range(len(OBS))]

    #Store the maximum and minimum offsets for first break picks and from headers
    for i in np.arange(len(OBS)):
        minoffsets_ps[i] = min(OBS[i][:,0])
        maxoffsets_ps[i] = max(OBS[i][:,0])
        minoffsets_hdr[i] = min(offset[i])
        maxoffsets_hdr[i] = max(offset[i])

    #Save the maximum and minimum indicies of first break picks w.r.t. segy header offsets
    minimum_idx, maximum_idx, a, b = [], [], 0, 0
    for i in np.arange(len(offset)):
        for j in np.arange(len(offset[i])):
            if minoffsets_ps[i] == offset[i][j]:
                minimum_idx.append(j)
            if maxoffsets_ps[i] == offset[i][j]:
                maximum_idx.append(j)
                
    return minimum_idx, maximum_idx

def data_cleanup(OBS, obsoffset, idx, reverse=True):
    """
    |_________________________________________________________________________________
    | PURPOSE:
    |-> Remove duplicate picks, i.e. two separate phases picked on the same channel.
    |-> Interpolate a first break pick for any channels between the minimum and maximum first break pick 
    |   where there is none present.
    |_________________________________________________________________________________
    | IMPORTANT:
    |-> idx is the min_idx, and max_idx where idx=[min_idx, max_idx]
    |-> reverse will flip the numbering of OBS, where the last OBS will now be indexed as the first.
    |_________________________________________________________________________________
    """
    
    #Sort data by offset
    tmp = []
    for i in range(len(OBS)):
        tmp.append(OBS[i][np.argsort(OBS[i][:,0])])
        
    #Only want unique values
    tmp0 = []
    for i in range(len(OBS)):
        tmp0.append(tmp[i][np.unique(tmp[i][:,0], return_index=True, axis=0)[1]])
        
    #Interpolate first breaks over min/max traces
    tmp1 = []
    for i in range(len(tmp0)):
        f = interpolate.interp1d(tmp0[i][:,0],tmp0[i][:,1])
        tmp1.append(np.vstack((obsoffset[i][idx[0][i]:idx[1][i]],
                              f(obsoffset[i][idx[0][i]:idx[1][i]]))).T)
        
    #Reverse shot order
    if reverse == True:
        OBS = []
        for i in range(len(tmp1)):
            it = len(tmp1) - 1 - i
            OBS.append(tmp1[it])
    else:
        OBS = tmp1
    
    #Output pick numbers, for now... #### SEE THIS
    pick_num = []
    for i in range(len(OBS)):
        pick_num.append(len(OBS[i]))
        
    return OBS, pick_num

def gen_toy2dac_acqi(data, obs_x, obs_z, bulk_shift, min_max, del_obs):
    """
    |_________________________________________________________________________________
    | PURPOSE:
    |-> Properly orgainize offset informtion in to a toy2dac geometry file, or "acqui" file
    |-> Delete anything outside of the minimum and maximum range here
    |-> Delete any undesirable OBS here as well
    |-> SEPT 2019 update, I really don't need two columns of data, just one. I've reduced to one here.
    |_________________________________________________________________________________
    """
    
    final_lst = [[] for _ in range(len(data))]
    for i in range(len(data)):
        #num = len(data) - 1
        final_lst[i] = np.zeros((len(data[i]), 5))#.astype(object) 

    for j in range(0,len(final_lst),1):
        for i in range(0, len(data[j]), 1):
            if i == 0: #Write the source line
                final_lst[j][i,0] = N(obs_z[j]*1000.0-1.0,9) #Source Depth in meters obs_z[j]*1000.0
                final_lst[j][i,1] = N((obs_x[j]*1000.0)+(bulk_shift*1000.0),9)
                final_lst[j][i,2] = 0
                final_lst[j][i,3] = 0
                final_lst[j][i,4] = 0
            else: #Write receiver lines
                final_lst[j][i,0] = N(0.005*1000.0,9) #Receiver Depth in meters
                final_lst[j][i,1] = N(data[j][i-1]*1000.0,9)
                final_lst[j][i,2] = 0 
                final_lst[j][i,3] = 0
                final_lst[j][i,4] = 1
                
    #Remove empty columns and report the number of OBS removed
    idx=[]
    for j in range(len(final_lst)):
        if len(final_lst[j]) == 0:
            idx.append(j)
        elif final_lst[j][0,1] < min_max[0]*1000.0:
            idx.append(j)
        elif final_lst[j][0,1] > min_max[1]*1000.0:
            idx.append(j)
        elif j in del_obs:
            idx.append(j)
    output = np.delete(final_lst,idx)
    
    #Ensure zero coordinate origin
    for i in range(len(final_lst)):
        final_lst[i][:,1] = final_lst[i][:,1]-(min_max[0]*1000.0)
        
    count = 0
    for i in range(len(output)):
        count += len(output[i][:,0])

    print "Orignal number of OBS: ", len(data)
    print "Remaining OBS based on maximum and minimum locations: ", len(data)-len(idx)
    print "Total number of traces for all remaining OBS: ", count
    print
                       
    return output, idx

"""
Chris Williams, September 2019
-> Function to upscale toy2dac geometry based on REAL geometry. It will retain the real data, but generate synthetic data inbetween real OBS.
PARAMS
upscale: Factor by which to increase the number of shots between existing shots (upscale*orig_size)-upscale
shot spacing: Shot spacing for generated data (or receiver spacing)
shot space sd noise: Gaussian noise is applied to synthetic shots to mimic real data, define the standard deviation of this noise
strd dev offset: Standard deviation to apply to offset constraints. Will be constrained by min/max offset and model size. Will be output * max offset.
minoffset/max offset: Minimum and maximum offsets to generate. Anything generated larger or smaller will be downsized to these bounds
model size: The size of the model, all in meters.
"""

def edit_acqui(acqui, upscale, shot_spacing, shot_space_sd_noise, strd_dev_offset, 
               min_offset, max_offset, model_size):
    
    print "Upscaling acquisition file by interpolating between true OBS."
    print "Upscaling factor: ",upscale

    #Upscale source geometry
    src = []
    for i in range(len(acqui)):
        src.append(acqui[i][0,1])
    xcorr = np.arange(0, len(src)*upscale, upscale)
    xcorr_update = np.arange((len(src)-1)*upscale)
    tmp = interpolate.interp1d(xcorr,src, kind="cubic")
    src_upscale = tmp(xcorr_update)
    
    print "Number of soruces: ", len(xcorr_update)

    #Project OBS down to seafloor
    b_file = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/000_DATA/3_SynthEmed/0_models/fbathy"
    with open(b_file, 'rb') as f:
        bathy = np.fromfile(f, dtype=np.float32)
    #Generate x coordinates to accompany bathymetry
    spacing, bcorr, src_upscale_z = 50.0, [], []
    for i in range(len(bathy)):
        bcorr.append(spacing*i)
    #Place "shots or OBS" on seafloor, use numpy interpolate
    for i in range(len(src_upscale)):
        src_upscale_z.append(np.interp(src_upscale[i], bcorr, bathy))

    #Generate list of offset ranges for each source
    min_offsets, max_offsets, count, idx = [], [], 0, 0
    for i in range(len(src_upscale)+1):
        if i in xcorr:
            min_offsets.append(min(acqui[idx][1::,1]-acqui[idx][0,1]))
            max_offsets.append(max(acqui[idx][1::,1]-acqui[idx][0,1]))
            idx+=1
        else:
            norm_min_offset = -abs(np.random.normal(0,strd_dev_offset,1)[0])*max_offset
            #Set to minimum offset if too low
            if abs(norm_min_offset) < min_offset:
                norm_min_offset = -min_offset
            #Set to maximum offset if too high
            if abs(norm_min_offset) > max_offset:
                norm_min_offset = -max_offset
            #Set to zero if minimum offset is too high
            if norm_min_offset + src_upscale[i] < 0:
                norm_min_offset = -abs(src_upscale[i]-50.0)

            norm_max_offset = abs(np.random.normal(0,strd_dev_offset,1)[0])*max_offset
            #Set to minimum offset if too low
            if norm_max_offset < min_offset:
                norm_max_offset = min_offset
            #Set to maximum offset if too high
            if norm_max_offset > max_offset:
                norm_max_offset = max_offset
            #Set to maximum model space if maximum offset is too high
            if norm_max_offset + src_upscale[i] > model_size:
                norm_max_offset = model_size-src_upscale[i]-50.0
            min_offsets.append(round(norm_min_offset,0))
            max_offsets.append(round(norm_max_offset,0))
            
    #Make offsets or shot locations
    shot_points = [[] for _ in range(len(src_upscale))]
    for i in range(len(shot_points)):
        shot_points[i] = np.arange(round(min_offsets[i],0),round(max_offsets[i],0),shot_spacing)
        count += len(shot_points[i])
        
    print "Total number of shots: ",count
    print

    #Add noise to shot locations
    for i in range(len(shot_points)):
        for j in range(len(shot_points[i])):
            shot_points[i] += np.random.normal(0,shot_space_sd_noise,1)[0]
        plt.plot(shot_points[i])
    plt.show()
            
    #Needs to be in KM, deal with obsured rounding
    for i in range(len(shot_points)):
        src_upscale[i] = round(src_upscale[i]/1000.,3)
        src_upscale_z[i] = round(src_upscale_z[i]/1000.,3)
        for j in range(len(shot_points[i])):
            shot_points[i][j] = round(shot_points[i][j]/1000. + src_upscale[i],3)
            
    return shot_points, src_upscale, src_upscale_z

###TEST WHAT I'VE DONE SO FAR###

import numpy as np
import segyio #This doesn't come standard with python
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.colors as mcolors
import struct

from scipy import interpolate
from sympy import N
from math import e

#############
### FILES ###
#############

###############
"""REAL DATA"""
###############

#Full path to tx.in file
f = "/home/chris/Documents/1_data/EasternMed/Chris_Williams/first_breaks/pick9.2a"
picks = np.genfromtxt(f) 

#OBS location file
f = "/home/chris/Documents/1_data/EasternMed/Chris_Williams/locations/obs.pqrsdz" 
obs = np.genfromtxt(f)
obs = np.sort(obs[:,0], axis=None)
for i in range(0,len(obs),1):
    obs[i] = (round(obs[i], 3))
    
#OBS depth file
f = "/home/chris/Documents/1_data/EasternMed/Chris_Williams/locations/OBS_Bathy.dat" 
bathy = np.genfromtxt(f)
obs_z = bathy[:,2]

#SEGY seismic data location and extension
segydir = '/home/chris/0_GlobeClaritas/0_Projects/Ref_EMed_Line1/COMMON/DATA/4_newfilt_all_Aug02_2019/'
name_start = ''
name_end = '_VERT_newfilt_celw.sgy'
fileloc = [segydir, name_start, name_end]

#################
### VARIABLES ###
#################

###############
"""REAL DATA"""
###############

#Define variables for read_fbpicks
phase = [0,1,3,6]
unc = 1
tt_window = [0,20.0]

#Define variables pertainig to segy data
sr = 0.004 #sample rate
tmax = 40.0 #maximum time
nsamples = tmax/sr +1 #number of time samples
maxfold = 1821 #maximum fold, maximum numbsetup_flager of traces

#Define variables for model cropping
offset_max = 100
offset_min = 0
bulk_shift = 0 # Not implimented
min_max = [30, 220]

#Indicies of OBS that have already been cropped out, blank if none: PYTHON INDEXING, OBS number - 1
del_obs = [8, 10]

#Define variables for acquisition upscaling
upscale = 4
shot_spacing = 100
shot_space_sd_noise = 3
strd_dev_offset = 0.5
minimum, maximum, size = 10000, 100000, 190000

################
### GEOMETRY ###
################

"""REAL DATA"""
#Read the first-beak picks    
OBS, pick_numbers = read_fbpicks(picks, obs, phase, unc, tt_window, sr)

#################
### SEGY DATA ###
#################

"""REAL DATA"""
#Extract obsoffsets from the header and the gathers from sgy file
gathers, obsoffsets = extract_segy(fileloc, len(obs), maxfold, nsamples)
#Indexes for first break picks in segy data
min_idx, max_idx = indices_from_segy(OBS, obsoffsets)

################
### CLEAN UP ###
################

"""REAL DATA"""
#Ensure uniqueness and interpolate gaps in dataset.
OBS, pick_num_data = data_cleanup(OBS, obsoffsets, [min_idx, max_idx])

####################
### MODEL EXTENT ###
####################

"""REAL DATA"""
#Crop the data and geometry together to ensure consistency.
#Set the maximum offset
OBS, indicesdel = max_offset(OBS, offset_max)
#Shift the data with the OBS model locations and a bulk shift
OBS = shift_data(OBS, obs, bulk_shift)
#Set the maximum and minimum values for model cropping
OBS, indicesdel_crop = crop_model(OBS, min_max)

#########################
### OUTPUT ACQUI FILE ###
#########################

OBS_edit = [[] for _ in range(len(OBS))]
for i in range(len(OBS)):
    OBS_edit[i]=OBS[i][:,0]

"""REAL DATA"""
#Output final version of acqui file
acqui, del_idx = gen_toy2dac_acqi(OBS_edit, obs, obs_z, bulk_shift, min_max, del_obs)

#print acqui

#Edit the acqui file (upscale)
#EVERY TIME YOU RUN THIS, IT WILL CHANGE, MAKE SURE YOU SAVE ACQUI FILES
shot_points, src_upscale, src_upscale_z = edit_acqui(acqui, upscale, shot_spacing, shot_space_sd_noise, 
                                                     strd_dev_offset, minimum, maximum, size)

del_obs_upscale = []
min_max = [0, 190000]
#Output edited version of acqui file
acqui, del_idx = gen_toy2dac_acqi(shot_points, src_upscale, src_upscale_z, bulk_shift, min_max, del_obs_upscale)

##################
### SAVE FILES ###
##################

#Save the full acqui file
for i in range(len(acqui)-1):
    acqui[0] = np.concatenate((acqui[0], acqui[i+1]), axis=0)
np.savetxt("/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/000_DATA/3_SynthEmed/3_models_fixcrop/2_acquidense/"+\
           "acqui_100", acqui[0], fmt='%9s', delimiter='\t')
#USE Edit_Acqui_TOY2DAC to make data 

#Plot Histogram of Acquisition File

offsets = [[] for _ in range(len(acqui))]
obs_loc = []

for i in range(len(acqui)):
    obs_loc.append(acqui[i][0,1])
    for j in range(1,len(acqui[i])):
        offsets[i].append(acqui[i][j,1]-obs_loc[i])
         
#Look at a histogram of our offset distribution
DATA, maxo = [], 0
off_1, tgt1 = 0, 15000
off_2, tgt2 = 0, 30000
off_3, tgt3 = 0, 45000
off_4, tgt4 = 0, 60000
off_5, tgt5 = 0, 100000
for i in range(len(offsets)):
    for j in range(len(offsets[i])):
        DATA.append(abs(offsets[i][j]))
        if abs(offsets[i][j]) > maxo:
            maxo = abs(offsets[i][j])
        if abs(offsets[i][j]) < tgt1:
            off_1 += 1
        if abs(offsets[i][j]) < tgt2:
            off_2 += 1
        if abs(offsets[i][j]) < tgt3:
            off_3 += 1
        if abs(offsets[i][j]) < tgt4:
            off_4 += 1
        if abs(offsets[i][j]) < tgt5:
            off_5 += 1
        
print "Maximum offset: ", maxo
print "Offsets group 1: ", off_1
print "Offsets group 2: ", off_2
print "Offsets group 3: ", off_3
print "Offsets group 4: ", off_4
print "Offsets group 5: ", off_5

fig, ax = plt.subplots()
        
n,bins,patches=ax.hist(DATA, 48, facecolor='blue')

fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 24
fig_size[1] = 24
plt.rcParams["figure.figsize"] = fig_size

ax.set_xlabel("Offset [m]",fontsize=24)
ax.set_ylabel("Number of offsets in bin",fontsize=24)
ax.set_title("Pick 9.2a offset distribution", fontsize=24)
ax.set_xticks([15000, 30000, 45000, 60000, 100000], minor=False)
#ax.set_xticks([10000, 30000, 50000, 70000, 90000], minor=True)
ax.xaxis.grid(True, which='major')
#ax.xaxis.grid(True, which='minor')
#plt.xticks(fontsize=24) 
#plt.yticks(fontsize=24)

plt.show()