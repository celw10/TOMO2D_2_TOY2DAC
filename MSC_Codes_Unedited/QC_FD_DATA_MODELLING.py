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
def gen_toy2dac_acqi(data, obs_x, obs_z, bulk_shift, min_max, del_obs, sr):
    """
    |_________________________________________________________________________________
    | PURPOSE:
    |-> Properly orgainize offset informtion in to a toy2dac geometry file, or "acqui" file
    |-> Delete anything outside of the minimum and maximum range here
    |-> Delete any undesirable OBS here as well
    |_________________________________________________________________________________
    """

    final_lst = [[] for _ in range(len(data))]
    for i in range(len(data)):
        #num = len(data) - 1
        final_lst[i] = np.zeros((len(data[i][:,0]), 5))#.astype(object) 

    for j in range(0,len(final_lst),1):
        for i in range(0, len(data[j][:,0]), 1):
            if i == 0: #Write the source line
                final_lst[j][i,0] = N(obs_z[j]*1000.0-1.0,9) #Source Depth in meters obs_z[j]*1000.0
                final_lst[j][i,1] = N((obs_x[j]*1000.0)+(bulk_shift*1000.0),9)
                final_lst[j][i,2] = 0
                final_lst[j][i,3] = 0
                final_lst[j][i,4] = 0
            else: #Write receiver lines
                final_lst[j][i,0] = N(0.005*1000.0,9) #Receiver Depth in meters
                final_lst[j][i,1] = N(data[j][i-1,0]*1000.0,9)
                final_lst[j][i,2] = round(data[j][i-1,1]*sr, 3) 
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

def setup_geom(f, setup_flag=1, sources=all):
    """
    -> Import toy2dac geometry file to read in the forward modeled data
    -> You have to use the geometry file you used to generate the data
    Params:
    -> f = full path to geometry file
    -> setup_flag = either 1 or 0. If 0, data are assumed to be constant fold. If 1, a variable fold is allowed,
    and the constant case is also allowed.
    -> sources = what source indicies to inverse fourier transform back to the time domain. Due to poor code design
    initially, going thorugh all sources would result in memory issues. Now, using only relevant frequencies, you
    should be able to inverse FT all sources. Specify the source index to just inverse FT one source, or specify all
    to inverse FT all sources."""

    ##############################
    ### IMPORT/EDIT ACQUI File ###
    ##############################

    if setup_flag != 0:
        #Full path to acqui file
        acqui = np.genfromtxt(f) #Save the tx.in file

        #Figure out acquisition info from acqui file
        shot_spacing = acqui[2,1]-acqui[1,1]#Maybe a better estimate for this in the future
        OBS = 0
        OBS_index = []
        shots = 0
        for i in range(0,len(acqui[:,1]),1):
            if acqui[i,4] == 0:
                OBS += 1
                OBS_index.append(i)
            else:
                shots += 1

        num_traces = []
        #Compute number of traces per OBS in geom. file
        for i in range(len(OBS_index)-1):
            num_traces.append(OBS_index[i+1]-OBS_index[i]-1)#-1 causing issues
            if i == len(OBS_index)-2:
                num_traces.append(len(acqui)-OBS_index[i+1]-1)
                
    #Print pertinent infromation
    if sources != all:
        print "EXTRACTING INFROMATION FOR A SINGLE SHOT, ", sources
    else:
        print "EXTRACTING INFORMATION FOR ALL SOURCES"
    print 

    if setup_flag != 0:
        print "SEISMIC DATA STATS..."
        print "Number of OBS: ", OBS
        print "Number of Shots: ", shots
        print "Shot Spacing: ", shot_spacing
        print "Average Number of Shots/OBS: ", shots/OBS
        print "Line Length: ", max(acqui[:,1]), "meters"
        print 
        #Specify the nchannels variable as the maximum fold
        nchannels = max(num_traces)
    else:
        print "DATA ASSUMED TO HAVE COMMON NUMBER OF CHANNELS PER SOURCE..."
    print
    
    if setup_flag != 0:
        return num_traces, nchannels
    else:
        "not", required
        
def restructure_toy2dac_data(data, nsources, nfreq, nrelfreq, nchannels, num_traces, sources=all):
    """
    -> This function restructures the toy2dac data_modeling file as a complex matrix of shape 
    [nsources, nfreq, maxfold]
    -> Any shot that is less than the maximum fold will be padded by zeros
    -> This will not affect the IFT taken in the vertical direction
    -> To successfully inverse FT the data, high frequencies are required, but forward modeling these frequencies 
    takes time. Therefore, only forward modeling relevant, low frequencies are required. The remainder of the frequencies
    will be zero. 
    -> A taper is later applied to smooth out this transition.
    Params:
    -> data = toy2dac datamodeling file
    -> nsources = number of shots
    -> nfreq = total number of discrete frequencies required for the inverse fourier transform (IFT)
    -> nrelfreq = number of relevant frequencies, i.e. discrete frequencies that were actually forward modeled in 
    toy2dac. 
    -> nchannels = maximum fold of the data, should be output form setup_geom previously
    -> num_traces = the fold for each shot, should be output from setup_geom previouisly
    -> sources = number of sources to restructure, all would restructure all sources, or enter the source index to 
    only restrucutre that source. 
    NEW UPDATES: DEC 05 19
    -> I've added functionality when only one source is modeled.
    """

    ########################
    ### RESTRUCTURE DATA ###
    ########################

    #FORM THE REAL MATRIX
    it=0
    print "RESTRUCTURING DATA_MODELING FILE..."
    #All sources
    if sources == all:
        RealMatrix= np.zeros((nsources, nfreq, nchannels))
        ComplexMatrix = np.zeros((nsources, nfreq, nchannels))
        #We don't have to read in all frequencies, only frequencies relevant to the data
        while it < relfreq:
            for s in range(0,nsources,1): 
                for c in range(0,2*num_traces[s],2):
                    RealMatrix[s,it,c/2] = data[c+(sum(num_traces[0:s]*2))+(it*sum(num_traces)*2)]
                    ComplexMatrix[s,it,c/2] = data[c+(sum(num_traces[0:s]*2))+(it*sum(num_traces)*2)+1]
            it+=1
    #One source
    elif sources == "one":
        RealMatrix= np.zeros((nfreq, nchannels))
        ComplexMatrix = np.zeros((nfreq, nchannels))
        while it < relfreq:
            for c in range(0,2*num_traces,2):
                RealMatrix[it,c/2] = data[c+((it*num_traces)*2)]
                ComplexMatrix[it,c/2] = data[c+(it*num_traces*2)+1]
            it+=1
    #A select sources out of many
    else:
        RealMatrix= np.zeros((nfreq, nchannels))
        ComplexMatrix = np.zeros((nfreq, nchannels))
        while it < relfreq:
            for c in range(0,2*num_traces[sources],2):
                RealMatrix[it,c/2] = data[c+(sum(num_traces[0:sources]*2))+(it*sum(num_traces)*2)]
                ComplexMatrix[it,c/2] = data[c+(sum(num_traces[0:sources]*2))+(it*sum(num_traces)*2)+1]
            it+=1
            
    #FORM THE COMPLEX MATRIX
    #All sources
    if sources == all:
        FD_Comp_SG = np.zeros((nsources, int(nfreq), nchannels), dtype=complex)
        for s in range(nsources):
            for it in range(int(relfreq)):
                for c in range(num_traces[s]):
                    FD_Comp_SG[s,it,c] = complex(RealMatrix[s,it,c],ComplexMatrix[s,it,c])
    elif sources == "one":
        FD_Comp_SG = np.zeros((int(nfreq), nchannels), dtype=complex)
        for it in range(int(relfreq)):
            for c in range(num_traces):
                FD_Comp_SG[it,c] = complex(RealMatrix[it,c],ComplexMatrix[it,c])
    else:
        FD_Comp_SG = np.zeros((int(nfreq), nchannels), dtype=complex)
        for it in range(int(relfreq)):
            for c in range(num_traces[sources]):
                FD_Comp_SG[it,c] = complex(RealMatrix[it,c],ComplexMatrix[it,c])
                
    print "RESTURCTURING COMPLETE"
    print
    
    return FD_Comp_SG

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

####################
### TOY2DAC DATA ###
####################

#Maximum TOY2DAC geometry file
f_geom = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/000_DATA/2_RealData/17_110StartMdl_THEEND/01_acqui_t0/acqui_100" 

####################
"""PREDICTED DATA"""
####################

#Open toy2dac data_modeling file
dirr = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/009_Thesis_Inversions/2_EMReal/0_Manu_Results_Final_043/8_fwdmdl/" #Path to root directory
model = "dm_fgfinal_mdl043"
filename = dirr + model #Name of binary file
with open(filename, 'rb') as f:
    data_pred = np.fromfile(f, dtype=np.float32)
    
#Plotting label
PredDataLabel = "Predicted Data"
    
###################
"""OBSERVED DATA"""
###################

#Open toy2dac data_modeling file
dirr = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/000_DATA/2_RealData/17_110StartMdl_THEEND/_FreqTest_V2/06_FreqProgressive0.25Hz_2Hz/" #Path to root directory
model = "DM_FG20_LC3_100"
filename = dirr + model #Name of binary file
with open(filename, 'rb') as f:
    data_obs = np.fromfile(f, dtype=np.float32)
    
#Plotting Label
ObsDataLabel = "Oberserved Data"
    
#########################
"""DATA FOR COMPARISON"""
#########################

#Open toy2dac data_modeling file
dirr = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/000_DATA/2_RealData/17_110StartMdl_THEEND/_0FinalFreqStartMdl/" #Path to root directory
model = "12_Feb0820_Start_LFC0.5_FP0.25Hz"
filename = dirr + model #Name of binary file
with open(filename, 'rb') as f:
    data_start = np.fromfile(f, dtype=np.float32)
    
#Plotting label
CompareDataLabel = "Starting Data"

#################
### VARIABLES ###
#################

####################
"""FOR ACQUI FILE"""
####################

#Full path to tx.in file
f_p = "/home/chris/Documents/1_data/EasternMed/Chris_Williams/first_breaks/pick9a.final.zerophase"
picks = np.genfromtxt(f_p) 

#OBS location file
f_ox = "/home/chris/Documents/1_data/EasternMed/Chris_Williams/locations/obs.pqrsdz" 
obs = np.genfromtxt(f_ox)
obs = np.sort(obs[:,0], axis=None)
for i in range(0,len(obs),1):
    obs[i] = (round(obs[i], 3))

#OBS depth file
f_oz = "/home/chris/Documents/1_data/EasternMed/Chris_Williams/locations/OBS_Bathy.dat" 
bathy = np.genfromtxt(f_oz)
obs_z = bathy[:,2]

#SEGY seismic data location and extension
segydir = '/home/chris/0_GlobeClaritas/0_Projects/Ref_EMed_Line1/COMMON/DATA/5_FinalFlow_Sept30th/3_FWIData_RMSNorm/'
name_start = ''
name_end = '_VERT_Final_FWI.sgy'
fileloc = [segydir, name_start, name_end]

#################
"""FOR FD DATA"""
#################

nsources = 16 #Number of sources 
nfreq = 22 #Number of discrete freqeuncy groups
relfreq = nfreq #Just do this so I don't need to change the function

##################
### MAKE ACQUI ###
##################

#Read the first-beak picks    
OBS, pick_numbers = read_fbpicks(picks, obs, [0,1,3,6], 1, [0,20.0], 0.004) #phase, unc, tt_window, SR

#Extract obsoffsets from the header and the gathers from sgy file
gathers, obsoffsets = extract_segy(fileloc, len(obs), 1821, 40.0/0.004+1) #maxfold, number of samples

del gathers

#Indexes for first break picks in segy data
min_idx, max_idx = indices_from_segy(OBS, obsoffsets)

#Clean up gathers
OBS, pick_num_data = data_cleanup(OBS, obsoffsets, [min_idx, max_idx])

del obsoffsets

#Set the maximum offset
OBS, indicesdel = max_offset(OBS, 100) #offset max

#Shift the data with the OBS model locations and a bulk shift
OBS = shift_data(OBS, obs, 0) #bulk shift

#Set the maximum and minimum values for model cropping
OBS, indicesdel_crop = crop_model(OBS, [30, 220]) #min/max

#Output final version of acqui file
acqui, del_idx = gen_toy2dac_acqi(OBS, obs, obs_z, 0, [30, 220], [0,1,2,11,13], 0.004) #Shift, min/max, del_obs, SR

###################
### FD DATA ORG ###
###################

#Setup geometry from toy2dac acqui file, setup flag=1 if non constant shots/OBS
num_traces, nchannels = setup_geom(f_geom, setup_flag=1, sources=all)

#Restructure toy2dac data_modeling file into matrix format
FD_Comp_SG = restructure_toy2dac_data(data_obs, nsources, nfreq, relfreq, nchannels, num_traces, sources=all)

FD_Obs_Traces = FD_Comp_SG.real

#################
### PRED DATA ###
#################

#Setup geometry from toy2dac acqui file, setup flag=1 if non constant shots/OBS
num_traces, nchannels = setup_geom(f_geom, setup_flag=1, sources=all)

#Restructure toy2dac data_modeling file into matrix format
FD_Comp_SG = restructure_toy2dac_data(data_pred, nsources, nfreq, relfreq, nchannels, num_traces, sources=all)

FD_Pred_Traces = FD_Comp_SG.real

##################
### START DATA ###
##################

#Setup geometry from toy2dac acqui file, setup flag=1 if non constant shots/OBS
num_traces, nchannels = setup_geom(f_geom, setup_flag=1, sources=all)

#Restructure toy2dac data_modeling file into matrix format
FD_Comp_SG = restructure_toy2dac_data(data_start, nsources, nfreq, relfreq, nchannels, num_traces, sources=all)

FD_Start_Traces = FD_Comp_SG.real

##################################
### PLOTTING - OBS & FREQUENCY ###
##################################

#Frequencies modeled
freq = [2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.25,4.5,4.75,5.0,5.25,5.5,5.75,6.0,6.25,6.5,6.75,7.0,7.25]

#Choose a frequency index
freq_idx = 12

#Choose a shot python idx
shot = 14

#Laplace Constant
TAU=0.5

#Scale the fontsize in the plots
fontscale=1.5
    
#Define Figure Size
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 24
fig_size[1] = 6
plt.rcParams["figure.figsize"] = fig_size    

#Setup subplots
#fig = plt.figure() 
#ax1 = fig.add_subplot(411)
#ax2 = fig.add_subplot(412)
#ax3 = fig.add_subplot(413)
#ax4 = fig.add_subplot(414)

#Setup subplots
fig = plt.figure() 
ax3 = fig.add_subplot(111)
#ax2 = fig.add_subplot(312)
#ax3 = fig.add_subplot(313)

#Figure Title
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "Shot "+str(shot+1)+" Frequency "+str(freq[freq_idx])

#Computed data-weights from the first breaks
fbs = [[] for _ in range(len(acqui))]
for i in range(len(acqui)):
    fbs[i].append(acqui[i][1::,2])
fb_weights = [[] for _ in range(len(acqui))]
for i in range(len(fbs)):
    fb_weights[i].append(np.exp(fbs[i][0][:]*TAU))

#########################################
"""First Plot: Data Weighting Function"""
#########################################

#Alter x_vals and dw_vals corresponding to DW function
x_vals = range(0,100001,5000)
#dw_vals = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00173, 0.00691, 0.01553, 0.02759, 0.04307, 0.06194, 0.08419, 0.10977, 0.13865, 0.1708, 0.20617, 0.24471, 0.28636, 0.33107, 0.37879, 0.42943, 0.48293, 0.53923, 0.59823, 0.65986, 0.72404, 0.79067, 0.85966, 0.93092, 1.00434, 1.07984, 1.1573, 1.23661, 1.31767, 1.40036, 1.48457, 1.57019, 1.65709, 1.74516, 1.83427, 1.9243, 2.01512, 2.10662, 2.19866, 2.29111, 2.38386, 2.47676, 2.5697, 2.66254, 2.75516, 2.84742, 2.9392, 3.03038, 3.12082, 3.21041, 3.29901, 3.38651, 3.47279, 3.55772, 3.64118, 3.72307, 3.80327, 3.88167, 3.95816, 4.03263, 4.10499, 4.17513, 4.24295, 4.30836, 4.37128, 4.4316, 4.48926, 4.54417, 4.59625, 4.64544, 4.69166, 4.73485, 4.77496, 4.81192, 4.84568, 4.8762, 4.90344, 4.92735, 4.94792, 4.96509, 4.97886, 4.98921, 4.99611, 4.99957, 4.9659, 4.69868, 4.1932, 3.50424, 2.70645, 1.88629, 1.13263, 0.52715, 0.13546, 0.0] #Change me
dw_vals = [0.0, 0.0, 0.11, 0.21, 0.32, 0.42, 0.53, 0.63, 0.74, 0.84, 0.95, 1.05, 1.16, 1.26, 1.37, 1.47, 1.58, 1.68, 1.79, 1.89, 2.0]

tmp_dw = interpolate.interp1d(x_vals, dw_vals, kind="linear")
dw_values = tmp_dw(abs(acqui[shot][1:len(acqui[shot]),1]-acqui[shot][0,1]))

#Make Plot
#ax1.plot((acqui[shot][1:len(acqui[shot]),1]-acqui[shot][0,1])/1000., dw_values, color="black")
#Labels
#ax1.set_ylabel("Amplitude", fontsize=24*fontscale)
#ax1.set_xlabel("Offset [km]",fontsize=24*fontscale)
#ax1.set_title("Data Weighting Function", fontsize=28*fontscale)
#ax1.tick_params(axis='both', which='major', labelsize=20*fontscale)
#ax1.tick_params(axis='both', which='minor', labelsize=20*fontscale)
#ax1.grid()

######################################
"""Second Plot: Un-Weighted FD Data"""
######################################

#Plots
#ax2.plot((acqui[shot][1:len(acqui[shot]),1]-acqui[shot][0,1])/1000.,
#         FD_Obs_Traces[shot][freq_idx][0:num_traces[shot]]*fb_weights[shot][0][:], 
#         lw=1, color="blue", label=ObsDataLabel)

#ax2.plot((acqui[shot][1:len(acqui[shot]),1]-acqui[shot][0,1])/1000.,
#         FD_Start_Traces[shot][freq_idx][0:num_traces[shot]]*fb_weights[shot][0][:], 
#         lw=1, ls="--", color="blue", label=CompareDataLabel) 

#ax2.plot((acqui[shot][1:len(acqui[shot]),1]-acqui[shot][0,1])/1000.,
#         FD_Pred_Traces[shot][freq_idx][0:num_traces[shot]]*fb_weights[shot][0][:], 
#         lw=1, ls="-", color="green", label=PredDataLabel)

#Labels
#ax2.set_ylabel("Amplitude", fontsize=24*fontscale)
#ax2.set_xlabel("Offset [km]",fontsize=24*fontscale)
#ax2.set_title("Un-Weighted FD Data", fontsize=28*fontscale)
#ax2.tick_params(axis='both', which='major', labelsize=20*fontscale)
#ax2.tick_params(axis='both', which='minor', labelsize=20*fontscale)
#ax2.legend(fontsize=14)
#ax2.set_ylim(-0.1,0.1)
#ax2.grid()

##################################
"""Third Plot: Weighted FD Data"""
##################################

ax3.plot((acqui[shot][1:len(acqui[shot]),1]-acqui[shot][0,1])/1000.,
         (FD_Obs_Traces[shot][freq_idx][0:num_traces[shot]])*dw_values*fb_weights[shot][0][:], 
         lw=1, ls="-", color="blue", label=ObsDataLabel)

ax3.plot((acqui[shot][1:len(acqui[shot]),1]-acqui[shot][0,1])/1000.,
         (FD_Start_Traces[shot][freq_idx][0:num_traces[shot]])*dw_values*fb_weights[shot][0][:], 
         lw=1, ls="--", color="red", label=CompareDataLabel) 

ax3.plot((acqui[shot][1:len(acqui[shot]),1]-acqui[shot][0,1])/1000.,
         (FD_Pred_Traces[shot][freq_idx][0:num_traces[shot]])*dw_values*fb_weights[shot][0][:], 
         lw=1, ls="-", color="green", label=PredDataLabel)

#Labels
ax3.set_ylabel("Amplitude", fontsize=24*fontscale)
ax3.set_xlabel("Offset [km]",fontsize=24*fontscale)
#ax3.set_title("Weighted FD Data", fontsize=28*fontscale)
ax3.tick_params(axis='both', which='major', labelsize=20*fontscale)
ax3.tick_params(axis='both', which='minor', labelsize=20*fontscale)
ax3.legend(fontsize=14*fontscale)
#ax3.set_ylim(-0.1,0.1)
ax3.grid()

#####################################
"""Fourth Plot: Weighted Residuals"""
#####################################

#Compute average residuals
#abs_startresid = abs(FD_Start_Traces[shot][freq_idx][0:num_traces[shot]]- \
#                     (FD_Obs_Traces[shot][freq_idx][0:num_traces[shot]]))*dw_values*fb_weights[shot][0][:]
#abs_mean_startresid = sum(abs_startresid)/len(acqui[shot][1:len(acqui[shot]),1]-acqui[shot][0,1])

#abs_predresid = abs(FD_Pred_Traces[shot][freq_idx][0:num_traces[shot]]- \
#                    (FD_Obs_Traces[shot][freq_idx][0:num_traces[shot]]))*dw_values*fb_weights[shot][0][:]
#abs_mean_predresid = sum(abs_predresid)/len(acqui[shot][1:len(acqui[shot]),1]-acqui[shot][0,1])

#mean_startresids = []
#mean_predresids = []
#for i in range(len(acqui[shot][1:len(acqui[shot]),1]-acqui[shot][0,1])):
#    mean_startresids.append(abs_mean_startresid)
#    mean_predresids.append(abs_mean_predresid)
    
#ax4.plot(acqui[shot][1:len(acqui[shot]),1]-acqui[shot][0,1],
#         abs(FD_Start_Traces[shot][freq_idx][0:num_traces[shot]]-(FD_Obs_Traces[shot][freq_idx][0:num_traces[shot]]))*dw_values, 
#         label=CompareDataLabel+" Residauls", lw=1, ls="--", color="blue") 

#ax4.plot(acqui[shot][1:len(acqui[shot]),1]-acqui[shot][0,1], mean_startresids, "b", label=CompareDataLabel)

#ax4.plot(acqui[shot][1:len(acqui[shot]),1]-acqui[shot][0,1],
#         abs(FD_Pred_Traces[shot][freq_idx][0:num_traces[shot]]-(FD_Obs_Traces[shot][freq_idx][0:num_traces[shot]]))*dw_values,  
#         label=PredDataLabel+" Residuals", lw=1, ls="--", color="green")

#ax4.plot(acqui[shot][1:len(acqui[shot]),1]-acqui[shot][0,1], mean_predresids, "g", label=PredDataLabel) #

#Labels
#ax4.set_ylabel("Amplitude", fontsize=24*fontscale)
#ax4.set_xlabel("Offset [km]",fontsize=24*fontscale)
#ax4.set_title("Weighted FD Data Residuals", fontsize=28*fontscale)
#ax4.tick_params(axis='both', which='major', labelsize=20*fontscale)
#ax4.tick_params(axis='both', which='minor', labelsize=20*fontscale)
#ax4.legend()
#ax4.grid()


fig.tight_layout()
plt.show()


#####################################
### PLOTTING - 2D OBS & FREQUENCY ###
#####################################

#Plot aspect ratio, <1 is stretched wide
pltaspect=1.5

#Fill a matrix with the residuals for the pred & compare data
Residuals_Start = np.zeros((len(acqui),len(freq)))
Residuals_Pred = np.zeros((len(acqui),len(freq)))

#Loop
for shot in range(len(acqui)):
    for freq_idx in range(len(freq)):
        
        #Alter x_vals and dw_vals corresponding to survey size
        tmp_dw = interpolate.interp1d(x_vals, dw_vals, kind="linear")
        dw_values = tmp_dw(abs(acqui[shot][1:len(acqui[shot]),1]-acqui[shot][0,1]))
        #Find starting (comparison) residuals
        abs_startresid = abs(FD_Start_Traces[shot][freq_idx][0:num_traces[shot]]- \
                             (FD_Obs_Traces[shot][freq_idx][0:num_traces[shot]]))*dw_values*fb_weights[shot][0][:]
        #Find predicted residuals
        abs_predresid = abs(FD_Pred_Traces[shot][freq_idx][0:num_traces[shot]]- \
                            (FD_Obs_Traces[shot][freq_idx][0:num_traces[shot]]))*dw_values*fb_weights[shot][0][:]
        #Fill residual matricies
        Residuals_Start[shot,freq_idx] = sum(abs_startresid)/len(acqui[shot][1:len(acqui[shot]),1]-acqui[shot][0,1])
        Residuals_Pred[shot,freq_idx] = sum(abs_predresid)/len(acqui[shot][1:len(acqui[shot]),1]-acqui[shot][0,1])

#Define Figure Size
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 16
fig_size[1] = 24
plt.rcParams["figure.figsize"] = fig_size    
        
#Setup subplots
fig = plt.figure() 
#ax5 = fig.add_subplot(311)
#ax6 = fig.add_subplot(312)
ax7 = fig.add_subplot(111)

#################################################
"""FIRST PLOT: STARTING (COPMARISON) RESIDUALS"""
#################################################

#pr_scale=max(Residuals_Pred.flatten())
#im1 = ax5.imshow(Residuals_Start, vmin=0, vmax=pr_scale, cmap="nipy_spectral_r", aspect=pltaspect)
#ax5.set_title(CompareDataLabel+" Absolute Residuals", fontsize=28)
#ax5.set_xlabel('Frequency [Hz]', fontsize=24)
#ax5.set_ylabel('OBS Number', fontsize=24)

#ax5.set_xticks(np.arange(0, len(freq), 5))
#ax5.set_xticklabels(freq[::5], fontsize=16)
#ax5.set_yticks(np.arange(0, len(acqui), 1))
#ax5.set_yticklabels(np.arange(1, len(acqui)+1, 1), fontsize=16)

#fig.colorbar(im1, ax=[ax5], location='right')
#cbar.ax.tick_params(labelsize=16.0)

######################################
"""SECOND PLOT: PREDICTED RESIDUALS"""
######################################

#im2 = ax6.imshow(Residuals_Pred, vmin=0, vmax=pr_scale, cmap="nipy_spectral_r", aspect=pltaspect)
#ax6.set_title(PredDataLabel+" Absolute Residuals", fontsize=28)
#ax6.set_xlabel('Frequency [Hz]', fontsize=24)
#ax6.set_ylabel('OBS Number', fontsize=24)

#ax6.set_xticks(np.arange(0, len(freq), 5))
#ax6.set_xticklabels(freq[::5], fontsize=16)
#ax6.set_yticks(np.arange(0, len(acqui), 1))
#ax6.set_yticklabels(np.arange(1, len(acqui)+1, 1), fontsize=16)

#fig.colorbar(im2, ax=[ax6], location='right')
#cbar.ax.tick_params(labelsize=16.0)

#############################################################
"""THIRD PLOT: PREDICTED - STARTING (COPMARISON) RESIDUALS"""
#############################################################

diff_scale=Residuals_Pred-Residuals_Start
im3 = ax7.imshow(Residuals_Pred-Residuals_Start, vmin=-max((diff_scale).flatten()),
                 vmax=max((diff_scale).flatten()), cmap="bwr", aspect=pltaspect)
#ax7.set_title(PredDataLabel+"-"+CompareDataLabel+" Residuals", fontsize=28)
ax7.set_xlabel('Frequency [Hz]', fontsize=24)
ax7.set_ylabel('OBS Number', fontsize=24)

ax7.set_xticks(np.arange(0, len(freq), 5))
ax7.set_xticklabels(freq[::5], fontsize=20)
ax7.set_yticks(np.arange(0, len(acqui), 1))
ax7.set_yticklabels(np.arange(1, len(acqui)+1, 1), fontsize=20)

cbar = fig.colorbar(im3, ax=[ax7], orientation='horizontal')
cbar.ax.tick_params(labelsize=18.0)
cbar.set_label("Amplitude", fontsize=20.0)

#######################################################
"""PRINT OUT THE STATS FROM THE COMPARISON OF MISFIT"""
#######################################################

#Print out: Number of cells with the smaller residual, mean residual, standard deviation of residuals
low_cells_compare = 0
low_cells_pred = 0
for shot in range(len(acqui)):
    for freq_idx in range(len(freq)):
        if Residuals_Pred[shot, freq_idx] > Residuals_Start[shot, freq_idx]:
            low_cells_compare += 1
        elif Residuals_Pred[shot, freq_idx] < Residuals_Start[shot, freq_idx]:
            low_cells_pred += 1
        else:
            print "Error, equal values"
            
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "~~~~~~~~~~ STATISTICS ~~~~~~~~~~"
print
print "Number of smaller avg. residuals:"
print "    Predicted Data: "+str(low_cells_pred)
print "    Starting/Predicted2 Data: "+str(low_cells_compare)
print "    Ratio: "+str(float(low_cells_pred)/float(low_cells_compare))
print "    Percentage: "+str((float(low_cells_pred)/(float(low_cells_compare)+float(low_cells_pred)))*100)
print 
print "Mean residual:"
print "    Predicted Data: "+str(np.mean(Residuals_Pred))
print
print "Standard Deviation:"
print "    Predicted Data: "+str(np.std(Residuals_Pred))
print
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

plt.show()

#####################################
### PLOTTING - 2D OBS & FREQUENCY ###
#####################################

"""WEIGHTED BY THE NUMBER OF TRACES IN EACH OBS"""

#Fill a matrix with the residuals for the pred & compare data
Residuals_Start_weight = np.zeros((len(acqui),len(freq)))
Residuals_Pred_weight = np.zeros((len(acqui),len(freq)))

#Build a weighting matrix based on number of picks per OBS num_picks/max(num_picks)
num_picks=np.zeros((len(acqui)))
weight=np.zeros((len(acqui)))
for i in range(len(acqui)):
    num_picks[i]=len(acqui[i])
for i in range(len(weight)):
    weight[i]=num_picks[i]/max(num_picks)

#Loop - fill residual matricies
for shot in range(len(acqui)):
    for freq_idx in range(len(freq)):
        
        #Alter x_vals and dw_vals corresponding to survey size
        tmp_dw = interpolate.interp1d(x_vals, dw_vals, kind="linear")
        dw_values = tmp_dw(abs(acqui[shot][1:len(acqui[shot]),1]-acqui[shot][0,1]))
        #Find starting (comparison) residuals
        abs_startresid = abs(FD_Start_Traces[shot][freq_idx][0:num_traces[shot]]- \
                             (FD_Obs_Traces[shot][freq_idx][0:num_traces[shot]]))*dw_values*fb_weights[shot][0][:]
        #Find predicted residuals
        abs_predresid = abs(FD_Pred_Traces[shot][freq_idx][0:num_traces[shot]]- \
                            (FD_Obs_Traces[shot][freq_idx][0:num_traces[shot]]))*dw_values*fb_weights[shot][0][:]
        #Fill residual matricies
        Residuals_Start_weight[shot,freq_idx] = sum(abs_startresid)/len(acqui[shot][1:len(acqui[shot]),1]-acqui[shot][0,1])*weight[shot]
        Residuals_Pred_weight[shot,freq_idx] = sum(abs_predresid)/len(acqui[shot][1:len(acqui[shot]),1]-acqui[shot][0,1])*weight[shot]

#Define Figure Size
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 16
fig_size[1] = 24
plt.rcParams["figure.figsize"] = fig_size    
        
#Setup subplots
fig = plt.figure() 
#ax5 = fig.add_subplot(311)
#ax6 = fig.add_subplot(312)
ax7 = fig.add_subplot(111)

#################################################
"""FIRST PLOT: STARTING (COPMARISON) RESIDUALS"""
#################################################

#im1 = ax5.imshow(Residuals_Start_weight, vmin=0, vmax=pr_scale, cmap="nipy_spectral_r", aspect=pltaspect)
#ax5.set_title(CompareDataLabel+" Residuals Weighted", fontsize=28)
#ax5.set_xlabel('Frequency [Hz]', fontsize=24)
#ax5.set_ylabel('OBS Number', fontsize=24)

#ax5.set_xticks(np.arange(0, len(freq), 5))
#ax5.set_xticklabels(freq[::5], fontsize=16)
#ax5.set_yticks(np.arange(0, len(acqui), 1))
#ax5.set_yticklabels(np.arange(1, len(acqui)+1, 1), fontsize=16)

#fig.colorbar(im1, ax=[ax5], location='right')
#cbar.ax.tick_params(labelsize=16.0)

######################################
"""SECOND PLOT: PREDICTED RESIDUALS"""
######################################

#im2 = ax6.imshow(Residuals_Pred_weight, vmin=0, vmax=pr_scale, cmap="nipy_spectral_r", aspect=pltaspect)
#ax6.set_title(PredDataLabel+" Residuals Weighted", fontsize=28)
#ax6.set_xlabel('Frequency [Hz]', fontsize=24)
#ax6.set_ylabel('OBS Number', fontsize=24)

#ax6.set_xticks(np.arange(0, len(freq), 5))
#ax6.set_xticklabels(freq[::5], fontsize=16)
#ax6.set_yticks(np.arange(0, len(acqui), 1))
#ax6.set_yticklabels(np.arange(1, len(acqui)+1, 1), fontsize=16)

#fig.colorbar(im2, ax=[ax6], location='right')
#cbar.ax.tick_params(labelsize=16.0)

#############################################################
"""THIRD PLOT: PREDICTED - STARTING (COPMARISON) RESIDUALS"""
#############################################################

im3 = ax7.imshow(Residuals_Pred_weight-Residuals_Start_weight, vmin=-0.02,#-max((diff_scale).flatten()),max((diff_scale).flatten())
                 vmax=0.02, cmap="bwr", aspect=pltaspect)
#ax7.set_title(PredDataLabel+"-"+CompareDataLabel+" Weighted", fontsize=28)
ax7.set_xlabel('Frequency [Hz]', fontsize=24)
ax7.set_ylabel('OBS Number', fontsize=24)

ax7.set_xticks(np.arange(0, len(freq), 4))
ax7.set_xticklabels(freq[::4], fontsize=20)
ax7.set_yticks(np.arange(0, len(acqui), 1))
ax7.set_yticklabels(np.arange(1, len(acqui)+1, 1), fontsize=20)

cbar = fig.colorbar(im3, ax=[ax7], orientation='horizontal')
cbar.ax.tick_params(labelsize=18.0)
cbar.set_label("Amplitude", fontsize=20.0)

#######################################################
"""PRINT OUT THE STATS FROM THE COMPARISON OF MISFIT"""
#######################################################

#Print out: Number of cells with the smaller residual, mean residual, standard deviation of residuals
low_cells_compare = 0
low_cells_pred = 0
for shot in range(len(acqui)):
    for freq_idx in range(len(freq)):
        if Residuals_Pred[shot, freq_idx] > Residuals_Start[shot, freq_idx]:
            low_cells_compare += 1
        elif Residuals_Pred[shot, freq_idx] < Residuals_Start[shot, freq_idx]:
            low_cells_pred += 1
        else:
            print "Error, equal values"
            
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "~~~~~~~~~~ STATISTICS ~~~~~~~~~~"
print
print "Number of smaller avg. residuals:"
print "    Predicted Data: "+str(low_cells_pred)
print "    Starting/Predicted2 Data: "+str(low_cells_compare)
print "    Ratio: "+str(float(low_cells_pred)/float(low_cells_compare))
print "    Percentage: "+str((float(low_cells_pred)/(float(low_cells_compare)+float(low_cells_pred)))*100)
print 
print "Mean residual:"
print "    Predicted Data: "+str(np.mean(Residuals_Pred))
print
print "Standard Deviation:"
print "    Predicted Data: "+str(np.std(Residuals_Pred))
print
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

plt.show()

import matplotlib.gridspec as gridspec

#Define colorbar range
plot_scale = 0.02

#Define Figure Size
size = [32,16]

#Re-define aspect ratio
pltaspect = 1

#Setup subplots
fig = plt.figure(figsize=(size[0], size[1]))
gs = gridspec.GridSpec(nrows=1, ncols=2,)

#############################################################
"""THIRD PLOT: PREDICTED - STARTING (COPMARISON) RESIDUALS"""
#############################################################

axa = fig.add_subplot(gs[0, 0])
im = axa.imshow(Residuals_Pred-Residuals_Start, vmin=-plot_scale,
                 vmax=plot_scale, cmap="bwr", aspect=pltaspect)
#axa.set_title(PredDataLabel+"-"+CompareDataLabel+" Residuals", fontsize=28)
axa.set_xlabel('Frequency [Hz]', fontsize=32*fontscale)
axa.set_ylabel('OBS Number', fontsize=32*fontscale)

axa.set_xticks(np.arange(0, len(freq), 2))
axa.set_xticklabels(freq[::2], fontsize=24*fontscale)
axa.set_yticks(np.arange(0, len(acqui), 1))
axa.set_yticklabels([4,5,6,7,8,9,10,11,13,15,16,17,18,19,20,21], fontsize=24*fontscale)

#############################################################
"""THIRD PLOT: PREDICTED - STARTING (COPMARISON) RESIDUALS"""
#############################################################

axb = fig.add_subplot(gs[0, 1])
im = axb.imshow(Residuals_Pred_weight-Residuals_Start_weight, vmin=-0.02,#-max((diff_scale).flatten()),max((diff_scale).flatten())
                 vmax=0.02, cmap="bwr", aspect=pltaspect)
#axb.set_title(PredDataLabel+"-"+CompareDataLabel+" Weighted", fontsize=28)
axb.set_xlabel('Frequency [Hz]', fontsize=32*fontscale)
#axb.set_ylabel('OBS Number', fontsize=24)

axb.set_xticks(np.arange(0, len(freq), 2))
axb.set_xticklabels(freq[::2], fontsize=24*fontscale)
axb.set_yticks(np.arange(0, len(acqui), 1))
axb.set_yticklabels([4,5,6,7,8,9,10,11,13,15,16,17,18,19,20,21], fontsize=24*fontscale)

#Colorbar Amplitude difference
ax_cbar = fig.add_axes([0.25, 0.05, 0.5, 0.04]) #[0.46, 0.355, 0.0175, 0.3] VERT #[0.053, -0.05, 0.434, 0.035] HORIZ W YAX
cbar = fig.colorbar(im, cax=ax_cbar, orientation='horizontal')
cbar.ax.tick_params(labelsize=24.0*fontscale)
cbar.set_label("Amplitude", fontsize=28.0*fontscale)

plt.tight_layout(w_pad=5,h_pad=0)

plt.show()