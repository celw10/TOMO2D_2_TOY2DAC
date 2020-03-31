#RUN ME IN JUPYTER NOTEBOOK, SEPARATE OUT THE FUNCTIONS FROM THE MAIN CODE AT THE END

### MAIN FUNCTION TO READ THE tx.in FB PICKS

def read_fbpicks(picks, shotloc, phase, unc, tt_window, SR=None, mode="TOMO2D"):
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
        if SR != None:
            OBS_sg.append(OBS[i][:,0:2])
            OBS_sg[i][:,1] = OBS_sg[i][:,1]/SR
        else:
            OBS_sg.append(OBS[i][:,0:2])
            OBS_sg[i][:,1] = OBS_sg[i][:,1]
            
    if mode == "TOY2DAC":
        #We only want time samples and position for TOY2DAC
        return OBS_sg, pick_numbers
    
    elif mode == "TOMO2D":
        #We need arrival time, position, uncertainty, and phase for TOMO2D
        return OBS, pick_numbers

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
    
def crop_model(OBS, min_max):
    """Specify the minimum and maximum extent of the model in the x dimension.
    To be applied after the data has been shifted to a consistent origin in shift_data above.
    min_max input as min_max = [xminimum, xmaximum]."""
    
    print "WARNING: maximum xcoordinate has been restricted to ", min_max[1]

    indicies_max = [[] for _ in range(len(OBS))]
    for i in np.arange(len(OBS)):
        for j in np.arange(len(OBS[i][:,0])):
            if OBS[i][j,0] >= min_max[1] and OBS[i][j,3] != 0: #We don't want to remove any sources... yet...
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
            if OBS[i][j,0] <= min_max[0] and OBS[i][j,3] != 0: #We don't want to remove any sources... yet...
                indicies_min[i].append(j)

    total=0
    for i in np.arange(len(OBS)):
        OBS[i] = np.delete(OBS[i],indicies_min[i],0)
        total += len(indicies_min[i])

    print "Traces removed: ", total
    print
        
    return OBS

def gen_tomo2d_acqui(data, obs_x, obs_z, min_max, del_obs, typ="inv"):
    """
    |_________________________________________________________________________________
    | PURPOSE:
    |-> Properly orgainize offset informtion in to a tomo2d geometry file
    |-> Accounting for position, arrival time, uncertainty, and phase
    |-> Delete anything outside of the minimum and maximum range here
    |-> Delete any undesirable OBS here as well
    |_________________________________________________________________________________
    """
    
    final_lst = [[] for _ in range(len(data)+1)]
    for i in range(len(data)+1):
        if i==0:
            #Add the row on top with num sources - i.e. header line
            final_lst[i] = np.zeros((1,6)).astype(object)
        else:
            final_lst[i] = np.zeros((len(data[i-1][:,0]), 6)).astype(object)
    
    #Write out the data
    for j in range(1,len(final_lst),1):
        for i in range(0, len(data[j-1][:,0]), 1):
            if i == 0: 
                #Write the source line
                final_lst[j][i,0] = "s" 
                final_lst[j][i,1] = obs_x[j-1] #x-position of OBS
                final_lst[j][i,2] = obs_z[j-1] #z-position of OBS
                final_lst[j][i,3] = None
                final_lst[j][i,4] = len(data[j-1][:,0])-1 #number of picks for OBS
                final_lst[j][i,5] = None
            elif typ == "geom": 
                #Write receiver lines for tomo2d geometry file
                final_lst[j][i,0] = "r"
                final_lst[j][i,1] = round(data[j-1][i,0], 3) #x-position of receiver
                final_lst[j][i,2] = 0.005 #z-position of receiver
                #NOTE: This is depended on how you define the phases of your data in tx.in file
                if data[j-1][i,3] == 1 or data[j-1][i,3] == 3 or data[j-1][i,3] == 6:
                    #Refractions
                    final_lst[j][i,3] = 0 #phase of arrival, reflection of refraction
                elif data[j-1][i,3] == 2 or data[j-1][i,3] == 4:
                    #Reflections
                    final_lst[j][i,3] = 1 #phase of arrival, reflection of refraction
                final_lst[j][i,4] = 0
                final_lst[j][i,5] = 0
            elif typ == "inv": 
                #Write receiver lines for tomo2d geometry file
                final_lst[j][i,0] = "r"
                final_lst[j][i,1] = round(data[j-1][i,0], 3) #x-position of receiver
                final_lst[j][i,2] = 0.005 #z-position of receiver
                #NOTE: This is depended on how you define the phases of your data in tx.in file
                if data[j-1][i,3] == 1 or data[j-1][i,3] == 3 or data[j-1][i,3] == 6:
                    #Refractions
                    final_lst[j][i,3] = 0 #phase of arrival, reflection of refraction
                elif data[j-1][i,3] == 2 or data[j-1][i,3] == 4:
                    #Reflections
                    final_lst[j][i,3] = 1 #phase of arrival, reflection of refraction
                final_lst[j][i,4] = round(data[j-1][i,1], 3) #first arrival time of wavefield
                final_lst[j][i,5] = round(data[j-1][i,2], 3) #uncertainty quantification
                
    #Remove empty columns and report the number of OBS removed
    #NOTE: No python indexing here as I added a column on top for the header
    idx=[]
    for j in range(1, len(final_lst), 1):
        if len(final_lst[j]) == 0:
            idx.append(j)
        elif final_lst[j][0,1] < min_max[0]:
            idx.append(j)
        elif final_lst[j][0,1] > min_max[1]:
            idx.append(j)
        elif j in del_obs:
            idx.append(j)

    #Add the header to the first row
    final_lst[0][0,:] = None
    final_lst[0][0,0] = len(final_lst)-len(idx)-1
    
    #Remove undesirable OBSdel_obs
    output = np.delete(final_lst,idx)
        
    count = 0
    for i in range(len(output)):
        count += len(output[i][:,0])

    print "Orignal number of OBS: ", len(data)
    print "Remaining OBS based on maximum and minimum locations: ", len(data)-len(idx)
    print "Total number of traces for all remaining OBS: ", count
    print
                       
    return output

### RUN VERSION 1.1 ###

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.colors as mcolors

#############
### FILES ###
#############

#Full path to tx.in file
f = "/home/chris/Documents/1_data/EasternMed/Chris_Williams/first_breaks/pick9.final.zerophase" 
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

#################
### VARIABLES ###
#################

#Define variables for read_fbpicks
phase = [0,1,3,6]
unc = 1
tt_window = [0,20.0]

#Define variables for model cropping
min_max = [30, 220]

#Indicies for OBS to be removed, NO PYTHON INDEXING HERE!
del_obs = [1, 2, 3]

#Indicate if a inversion or geometry file is being produced, inv or geom
typ = "geom"

####################
### FIRST BREAKS ###
####################

#Read the first-beak picks    
OBS, pick_numbers = read_fbpicks(picks, obs, phase, unc, tt_window)

print OBS[0]

####################
### MODEL EXTENT ###
####################

#Set the maximum and minimum values for model cropping
OBS = crop_model(OBS, min_max)

#########################
### OUTPUT ACQUI FILE ###
#########################

#Future me, who sees this mess, I fudged it to output a unstructured geometry file

#Output final version of tomo2d geometry file
acqui = gen_tomo2d_acqui(OBS, obs, obs_z, min_max, del_obs, typ=typ)

#Remove undesirable OBSdel_obs
print obs,obs_z
del_obs = [0, 1, 2, 11, 13]
OBS = np.delete(OBS,del_obs)
obs = np.delete(obs,del_obs)
obs_z = np.delete(obs_z,del_obs)
print np.shape(OBS), np.shape(obs), np.shape(obs_z)

"""
This is pretty cool. I can make an acquisition file of the picks that fail cycle-skipping
"""
#Optional, only want to keep the points that didn't pass cycle skipping, for some visual entertainment
#bad_idx = np.load("/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/110_pick9.final.zerophase_108run_Oct29/3_modeled_traveltimes/good_idx_2.5Hz.npy")  
#acqui_bad_idx = []#[[] for _ in range(len(bad_idx))]

for i in range(len(bad_idx)):
    q=0
    mat = np.zeros((len(bad_idx[i]),4))
    for j in range(len(OBS[i])):
        if j in bad_idx[i]:
            mat[q,0],mat[q,1],mat[q,2],mat[q,3] = OBS[i][j][0],OBS[i][j][1],OBS[i][j][2],OBS[i][j][3]
            q+=1
    acqui_bad_idx.append(mat)
Output final version of tomo2d geometry file
acqui_bad = gen_tomo2d_acqui(acqui_bad_idx, obs, obs_z, [-30,240], [], typ=typ)
        
############
### SAVE ###
############

#Save everyting into one file (FOR MODIFIED ACQUI)
for i in range(1,len(acqui_bad),1):
    acqui_bad[0] = np.concatenate((acqui_bad[0], acqui_bad[i]), axis=0)

#Save everyting into one file (FOR ORIG ACQUI)
for i in range(1,len(acqui),1):
    acqui[0] = np.concatenate((acqui[0], acqui[i]), axis=0)

if typ=="inv":
    np.savetxt("/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/96_pick9.final.zerophase_94run_Oct28/picks_inv/pg.0_20.all",
               acqui[0], fmt='%5s', delimiter='\t') #Give a name
if typ=="geom":
    np.savetxt("/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/110_pick9.final.zerophase_108run_Oct29/picks_fwd/pg.cycle.skip.good.2.5Hz",
               acqui_bad[0], fmt='%5s', delimiter='\t') #Give a name