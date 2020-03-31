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
            
    #UNIQUE TO FORWARD MODELING: Delete the source positions in observed data
    for o in np.arange(len(shotloc)):
        for i in np.arange(len(OBS[o])-1):
            if OBS[o][i,3] == 0:
                OBS[o] = np.delete(OBS[o],i,0)

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

def synthetic_fwdmdl_picks(fi, mode="true"):
    """
    A function to read in synthetic forward modeled traveltimes, converted from TOMO2D, that are considered
    to be the true FB picks. 
    I figured how to read in a unstructured file here that'll be useful for future endavours...
    -> Mode: If "true", restructuring the true synthetic forward model traveltimes
    -> Mode: If "model", restructuring the forward model traveltimes from a starting or inverted model
    """

    #Read in an unstructured text file (APPLY THIS TO OTHER CODES DO I DON'T HAVE TO ADD THE DAMN ZEROS)
    out = []
    for line in fi:
        line = line.strip()
        columns = line.split()
        if columns[0] == "s":
            out_tmp = []
            if line != 1:
                out.append(out_tmp)
            if mode == "model":
                out_tmp.append([str(columns[0])])
        if columns[0] == "r":
            pos = float(columns[1])
            tt = float(columns[4])
            out_tmp.append([pos, tt])

    #Restructure to the desired OBS format
    #First some indexing stuff
    maxout, maxouttmp, totalout, indexout = 0, 0, 0, []
    for i in range(np.shape(out)[0]):
        if i == 0:
            indexout.append(maxouttmp)
        else:
            indexout.append(maxouttmp+indexout[i-1])
        maxouttmp = np.shape(out[i])[0]
        if maxouttmp > maxout:
            maxout = maxouttmp
        totalout+=maxouttmp
    #Save in proper format
    if mode == "true":
        OBS = np.zeros((np.shape(out)[0],totalout,4))
        for o in range(np.shape(out)[0]):
            for r in range(np.shape(out[o])[0]):
                OBS[o,r,0] = out[o][r][0]
                OBS[o,r,1] = out[o][r][1]
    if mode == "model":
        OBS = np.zeros((totalout,6))
        for o in range(np.shape(out)[0]):
            for r in range(np.shape(out[o])[0]):
                if out[o][r][0] == "s":
                    print "Skipping row, source ", o, " detected."+"\t",
                else:
                    OBS[indexout[o]+r,1] = out[o][r][0]
                    OBS[indexout[o]+r,4] = out[o][r][1]
    print
    
    return OBS

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

def read_fwdmdl_picks(forward, numobs):
    """
    A code to read forward modeled traveltimes from TOMO2D:
    NOTE: The forward model file output by TOMO2D will have 4 out of the 6 column entries where there is a source.
    Numpy matricies require full data entries. Therefore, I fudged these outputs placing two zero entries at each
    source line. This is completely managable for me as I only have ~20 sources. If you have more and wish not to
    do this, you'll have to re-write the code NOT using np.genfromtxt. 
    -> And now I've returned with 144 sources and not the time for this. For synthetic files, I'll design a code
    that adds the zeros at each source automatically.
    """
    
    indices_fwd = []
    for i in range(0, len(forward), 1): #Want to figure out where the breaks in the forward modeled file are
        if forward[i,4] == 0 and forward[i,5] == 0: #Recall we set these digits to 0 0
            indices_fwd.append(i)
            
    OBS_fw = [[] for _ in range(numobs)]
    for i in np.arange(numobs):
        if i == numobs-1:
            OBS_fw[i] = forward[indices_fwd[i]:len(forward)]
        else:
            OBS_fw[i] = forward[indices_fwd[i]:indices_fwd[i+1]]

    #Delete "s" and "r" values and source rows
    for o in np.arange(numobs):
        for i in np.arange(len(OBS_fw[o])-1): #MINUS 1???
            if OBS_fw[o][i,4] == 0 and OBS_fw[o][i,5] == 0: 
                OBS_fw[o] = np.delete(OBS_fw[o],i,0)
                
    return OBS_fw

def calc_cycleskip_risk(OBS, OBS_fw, freq):
    """
    Compute the cycle skipping risk based on picked and modeled traveltimes.
    -> Convert frequency to seconds/cycle, and half this.
    -> If the difference between the picked and modeled traveltimes are greater than 0.5 seconds/cycle for the
       desired frequency, then the cycle skipping criteria is not fufilled. 
    """

    #Cycle skip threshold
    wl_t = 1.0/(2.0*freq)

    #Differences between modeled and observed traveltimes
    OBS_diff = [[] for _ in range(len(OBS))]
    OBS_pos = [[] for _ in range(len(OBS))]
    for i in np.arange(len(OBS)):
        OBStmp = OBS[i][OBS[i][:,0].argsort()]
        OBS_fwtmp = OBS_fw[i][OBS_fw[i][:,1].argsort()]
        if len(OBStmp) > len(OBS_fwtmp):
            n = 0
            while OBStmp[n,0] != OBS_fwtmp[0,1]: #Find out if there are unmatched values prior too
                n += 1
            for j in np.arange(len(OBS_fwtmp)):
                assert OBStmp[j+n][0] == OBS_fwtmp[j][1]
                OBS_pos[i].append(OBS_fwtmp[j][1])
                OBS_diff[i].append(abs(OBStmp[j+n][1] - OBS_fwtmp[j][4]))   
        elif len(OBStmp) < len(OBS_fwtmp):
            n = 0
            while OBS_fwtmp[n,0] != OBStmp[0,1]: #Find out if there are unmatched values prior too
                n += 1
            for j in np.arange(len(OBStmp)):
                assert OBStmp[j][0] == OBS_fwtmp[j+n][1]
                OBS_pos[i].append(OBStmp[j][0])
                OBS_diff[i].append(abs(OBStmp[j][1] - OBS_fwtmp[j+n][4]))
        else:
            for j in np.arange(len(OBStmp)):
                assert OBStmp[j][0] == OBS_fwtmp[j][1]
                OBS_pos[i].append(OBStmp[j][0])
                OBS_diff[i].append(abs(OBStmp[j][1] - OBS_fwtmp[j][4]))

    #List of model positions were every shot was taken
    modelpositions = []
    
    for sublist in OBS_pos:
        for item in sublist:
            modelpositions.append(item)
    modelpositions = sorted(modelpositions)
    
    #NEW ADDITION, I HAD DUPLICATES, CHECK WITH REAL DATA!
    modelpositions = list(OrderedDict.fromkeys(modelpositions))

    #Matrix to hold cycle skipping values
    cycle_skip = np.zeros((len(OBS), len(modelpositions)))
    cycle_skip = cycle_skip
    tt_residual = np.zeros((len(OBS), len(modelpositions)))
    cycle_skip = cycle_skip 
    for i in np.arange(len(OBS)):
        n = 0
        for j in np.arange(len(OBS_diff[i])):
            if OBS_pos[i][j] == modelpositions[n]: #If the positions are equal...
                cycle_skip[i,n] = (wl_t - abs(OBS_diff[i][j]))
                n += 1
            else: #Increase modelposition until it matches OBS_pos (not a pick at every shot?)
                while OBS_pos[i][j] != modelpositions[n]:
                    n += 1
                cycle_skip[i,n] = (wl_t - abs(OBS_diff[i][j]))

        for j in np.arange(len((modelpositions))): #Find the max and min indicies for interpolation
            if modelpositions[j] == min(OBS_pos[i]):
                index_l = j
            if modelpositions[j] == max(OBS_pos[i]):
                index_h = j

        x = np.arange(len(modelpositions)) #Interpolate between first and last first break pick for each OBS
        idx = np.nonzero(cycle_skip[i]) #Only interpolate zero values
        interp = interp1d(x[idx], cycle_skip[i][idx],  fill_value="extrapolate")
        tt_residual[i,index_l:index_h] = interp(range(index_l, index_h))

        for j in np.arange(len((modelpositions))):
            if tt_residual[i,j] == 0:
                tt_residual[i,j] = np.nan #-wl_t*10.0

    return tt_residual, wl_t, modelpositions, OBS_diff

import matplotlib as mpl

def plot_cycleskip_risk(tt_residual, wl_t, loc, numobs, freq, savepath=None, savename=None,
                        fontscalar=0, size=[26,13]):
    """
    Define a unique colorbar and a unique plot to highlight the risk of cycle skipping to the data. 
    """
    
    ##################################################
    ### Make a colorbar for traveltime misfit plot ###
    ##################################################

    def make_colormap(seq):
        """Return a LinearSegmentedColormap
        seq: a sequence of floats and RGB-tuples. The floats should be increasing
        and in the interval (0,1).
        """
        seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
        cdict = {'red': [], 'green': [], 'blue': []}
        for i, item in enumerate(seq):
            if isinstance(item, float):
                r1, g1, b1 = seq[i - 1]
                r2, g2, b2 = seq[i + 1]
                cdict['red'].append([item, r1, r2])
                cdict['green'].append([item, g1, g2])
                cdict['blue'].append([item, b1, b2])
        return mcolors.LinearSegmentedColormap('CustomMap', cdict)

    #Edit this bit if you want to play with the colorbar
    c = mcolors.ColorConverter().to_rgb
    tt_misfit = make_colormap(
        [c('grey'), 0.25, c('darkseagreen'), 0.5, c('yellowgreen'), 0.75, c('forestgreen')])
        #[c('white'), 0.01, c('lightsteelblue'), 0.625, c('royalblue'), 0.75, c('blue'), 0.875, c('mediumblue')])
    #    tt_misfit = make_colormap([c('white'), 0.01, c('red'), 0.88, c('red'), c('orange'), 0.90, c('orange'), c('yellow'),
    #                           0.92, c('yellow'), c('green')])
 
    #Customize colorbar
    #bounds=[-wl_t,-wl_t*0.5,0,wl_t*0.5,wl_t]
    #norm = mcolors.BoundaryNorm(bounds, tt_misfit.N)

    ###################################
    ### Plot cycle skipping summary ###
    ###################################

    fig, ax = plt.subplots(figsize=(size[0],size[1]))

    cax = plt.imshow(tt_residual, cmap=tt_misfit, vmin = -wl_t, vmax = wl_t, aspect='auto')#,
                    #norm=norm)#-(wl_t*9)

    #Figure out what shot number corresponds to an OBS location
    for x in np.arange(len(loc)):
        n = 0
        while loc[x] > modelpositions[n]:
            n += 1
        #plt.axvline(n)
        line = mlines.Line2D([n, n], [(-0.5+1*x), (0.5+1*x)], color="black", linewidth=2*fontscalar)
        ax.add_line(line)

    #plt.title("FWI Inital Model Assesment, frequency: " + str(freq) + " Hz", fontweight="bold", 
    #          fontsize=20*fontscalar)
    plt.xlabel('X-Position [km]', fontsize = 18*fontscalar)
    plt.ylabel('OBS', fontsize = 18*fontscalar)
    
    #Customized colorbar
    cbar = plt.colorbar(cax, ticks=[(-wl_t*0.75),-(wl_t*0.25), (wl_t*0.25),(wl_t*0.75)], 
                        orientation='horizontal', pad=0.06)
    cbar.ax.set_xticklabels(['< -1*Lambda', '< -3/4*Lambda', '< -1/2*Lambda', '< -1/4*Lambda'])
    #cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=18*fontscalar)

    """Format Axis"""
    #For Starting FAT model 1
    #plt.xticks(np.arange(0,len(tt_residual[0,:]),len(tt_residual[0,:])/(len(np.arange(-30,230,20))-0.1)), 
    #           np.arange(-30,230,20), fontsize=16*fontscalar)
    #plt.yticks(np.arange(0,21,1), np.arange(1, numobs+1, 1), fontsize=16*fontscalar)

    #For Starting FAT model 2
    #plt.xticks(np.arange(0,len(tt_residual[0,:]),len(tt_residual[0,:])/(len(np.arange(25,225,25))-0.1)), 
    #           np.arange(25,225,25), fontsize=16*fontscalar)
    #plt.yticks(np.arange(0,16,1), [4,5,6,7,8,9,10,11,13,15,16,17,18,19,20,21], fontsize=16*fontscalar)
    
    #For Starting FAT model 3 - Appendices
    plt.xticks(np.arange(0,len(tt_residual[0,:]),len(tt_residual[0,:])/(len(np.arange(30,220,20))-0.1)), 
               np.arange(30,220,20), fontsize=16*fontscalar)
    plt.yticks(np.arange(0,16,1), [4,5,6,7,8,9,10,11,13,15,16,17,18,19,20,21], fontsize=16*fontscalar)
    
    #For real data final crop configuration (Starting FAT model three) - MANUSCRIPT
    #plt.xticks(np.arange(0,len(tt_residual[0,:]),len(tt_residual[0,:])/(len(np.arange(0,190,20)))), 
    #           np.arange(0,190,20), fontsize=16*fontscalar)
    #plt.yticks(np.arange(0,16,1), [4,5,6,7,8,9,10,11,13,15,16,17,18,19,20,21], fontsize=16*fontscalar)
    
    ax.yaxis.grid()
    
    if savepath != None and savename != None:
        plt.gcf().savefig(savepath+savename+".png", bbox_inches = 'tight', pad_inches = 0.1)
        
    plt.tight_layout()
        
    plt.show()

### RUN VERSION 1.1 ###

import sys
import numpy as np
from scipy.interpolate import interp1d
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.colors as mcolors
import matplotlib.lines as mlines

#############
### FILES ###
#############

#Mode geometry, problem specific source geometry -> mode=0, generic x,z source geometry file, mode=1
mode_geometry = 1

#Mode first breaks, first breaks from tx.in file -> mode=0, first breaks from tomo2d geometry file, mode=1
mode_fb = 1

#Mode forward modeled data, zeros added to sources -> mode=0, add zeros to sources, mode=1
mode_fwd = 1

#Full path to tx.in file (the TRUE dataset)
if mode_fb == 0:
    f = "/home/chris/Documents/1_data/EasternMed/Chris_Williams/first_breaks/pick9.final.zerophase" 
    #f = "/home/chris/Documents/1_data/EasternMed/Chris_Williams/first_breaks/pick9.2" 
    picks = np.genfromtxt(f)
if mode_fb ==1:
    fi = open("/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/00_TOMO2DCONVERTED/2_EmedSynthetic/1_NovStartMdl_FBPickZeroPh/Forward_TrueMdl_EMSynth", "r")
    OBS = synthetic_fwdmdl_picks(fi, mode="true")
    
#Full path to forward modeled traveltimes (the MODELED dataset)
if mode_fwd == 0:
    f="/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/110_pick9.final.zerophase_108run_Oct29/3_modeled_traveltimes/09_fwdtts"
    #f="/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/66_pick9.2_65run_Aug02/3_modeled_traveltimes/09_fwdtts"
    forward=np.genfromtxt(f, skip_header=1)
elif mode_fwd == 1:
    fi = open("/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/00_TOMO2DCONVERTED/2_EmedSynthetic/1_NovStartMdl_FBPickZeroPh/Forward_StartMdl_EMSynth", "r")
    forward=synthetic_fwdmdl_picks(fi, mode="model")

if mode_geometry == 0:
    #Build the source geometry file from avaliable dataset
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
if mode_geometry == 1:
    #Build the source geometry file from 2 column [x,z] ascii file
    f = "/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/00_TOMO2DCONVERTED/2_EmedSynthetic/0_FinalSourceGeoms/Synthetic_SourceGeom_Final_True"
    source_geom = np.genfromtxt(f)
    obs = source_geom[:,0]
    obs_z = source_geom[:,1]
    
#Path to save file, NONE if you don't want to save it. Both savepath and savename have to be defined.
savepath = None

#Name to save file, .png extension will be automatically used.
savename = None

#################
### VARIABLES ###
#################

#Define variables for read_fbpicks
phase = [0,1,3,6]
unc = 1
tt_window = [0,20.0]

#Define variables for model cropping
min_max = [-30,30]#[25,225]#[-30,260]#[30,220]

#Indicies for OBS to be removed, PYTHON INDEXING, -1!!
del_obs = []#[0, 1, 2, 11, 13]#For real data tomography, [] for synthetic examples 

#Frequency to compute cycle skipping risk plots
freq = 0.75

############################
### PLOTTING DEFINITIONS ###
############################

#Figure size, 3 times normal for high res image
size = [24,24]

#Multiply the fontsize by this number
fontscale = 2.6 #2.6 for manuscript, same for the histogram below, 2.75 for Appendix montage images

####################
### FIRST BREAKS ###
####################

#np.set_printoptions(threshold=sys.maxsize)

#Adjust first breaks for tx.in file
if mode_fb == 0:
    OBS, pick_numbers = read_fbpicks(picks, obs, phase, unc, tt_window)
    OBS = crop_model(OBS, min_max)
    OBS = np.delete(OBS, del_obs)

#Adjust first breaks synthetic forward modeled tomo2d file
if mode_fb == 1:
    print np.shape(OBS)
    OBS = crop_model(OBS, min_max)
    OBS = np.delete(OBS, del_obs, axis=0)
    print np.shape(OBS)
    
#Delete any omitted obs
obs = np.delete(obs, del_obs)

#Read in the forward modeled arrivals
OBS_fw = read_fwdmdl_picks(forward, len(obs))

###################################
### PRESENT CYCLE SKIPPING RISK ###
###################################

#Comput cycle skipiping risk for a select frequency
tt_residual, wl_t, modelpositions, obs_diff = calc_cycleskip_risk(OBS, OBS_fw, freq)

#Make a plot of the cycle skip risk
plot_cycleskip_risk(tt_residual, wl_t, obs, len(OBS), freq, savepath=savepath, savename=savename,
                    fontscalar=fontscale, size=size)

#Some QC stuff
bad_idx = [[] for _ in range(len(OBS))]
tot = 0
idxtot = 0
numbad = 0
for i in range(len(obs_diff)):
    for j in range(len(obs_diff[i])):
        if obs_diff[i][j] <= wl_t:
            bad_idx[i].append(j)
            tot += obs_diff[i][j]
            numbad += 1
    tot+=round(sum(obs_diff[i])*1000.0,1)
    idxtot+=len(obs_diff[i])
    
#Save the aveage
avg_misfit = round(tot/idxtot,1)

#Save the bad indicies, that didn't pass cycle skipping for the given frequency
#np.save("/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/110_pick9.final.zerophase_108run_Oct29/3_modeled_traveltimes/good_idx_2.5Hz",
#            bad_idx)

print "Total misfit, ms: ", tot
print "Number of picks: ", idxtot
print "Number of picks at risk of cycle-skipping: ", numbad
print "RMS Misfit, ms: ", avg_misfit

#Plot Histogram of traveltime misfits

#A buffer to move the text off the cycle skipping lines
txt_buff = 5
#How far up to plot the text
max_num_per_bin = 1400
#The maximum range to plot
max_xrange = 550
#The width of the histogram bins
bin_width = 10
#Fontscale
#fontscale=1.0

residuals = []
for i in range(len(obs_diff)):
    for j in range(len(obs_diff[i])):
        residuals.append(obs_diff[i][j]*1000.0)

fig, ax = plt.subplots()
        
n,bins,patches=ax.hist(residuals, bins=range(0,550,10), facecolor='blue', edgecolor='black', linewidth=1.0*fontscale)

fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 24
fig_size[1] = 24
plt.rcParams["figure.figsize"] = fig_size

ax.set_xlabel("Residuals [ms]",fontsize=18*fontscale)
ax.set_ylabel("Number Observed",fontsize=18*fontscale)
#ax.set_title("Distribution of Pick Residuals", fontsize=28*fontscale, fontweight="bold")

plt.xlim(0,max_xrange)
plt.xticks(fontsize=16*fontscale)
plt.yticks(fontsize=16*fontscale)

"""Plot Mean Residual"""

#Average
plt.axvline(x=avg_misfit, color="black")
plt.text(avg_misfit+txt_buff, max_num_per_bin, 'Average Residual: '+str(avg_misfit)+" [ms]", rotation=90, 
         fontsize=16*fontscale, fontweight="bold")

"""Cycle-skipping thresholds"""

#1.0 Hz (ms)
wl_t_10 = (1.0/(2.0*1.0))*1000.0
plt.axvline(x=wl_t_10, color="black", ls="--")
plt.text(wl_t_10+txt_buff, max_num_per_bin, r'$\bf{1.0 Hz}$'+' Cycle-Skipping',rotation=90, fontsize=16*fontscale)

#1.5 Hz (ms)
#wl_t_15 = (1.0/(2.0*1.5))*1000.0
#plt.axvline(x=wl_t_15, color="black", ls="--")
#plt.text(wl_t_15+txt_buff, max_num_per_bin, '1.5 Hz Cycle-Skipping',rotation=90, fontsize=16*fontscale)

#2.0 Hz (ms)
wl_t_20 = (1.0/(2.0*2.0))*1000.0
plt.axvline(x=wl_t_20, color="black", ls="--")
plt.text(wl_t_20+txt_buff, max_num_per_bin, r'$\bf{2.0 Hz}$'+' Cycle-Skipping',rotation=90, fontsize=16*fontscale)

#2.5 Hz (ms)
#wl_t_25 = (1.0/(2.0*2.5))*1000.0
#plt.axvline(x=wl_t_25, color="black", ls="--")
#plt.text(wl_t_25+txt_buff, max_num_per_bin, '2.5 Hz Cycle-Skipping',rotation=90, fontsize=16*fontscale)

#3.0 Hz (ms)
#wl_t_30 = (1.0/(2.0*3.0))*1000.0
#plt.axvline(x=wl_t_30, color="black", ls="--")
#plt.text(wl_t_30+txt_buff, max_num_per_bin, '3.0 Hz Cycle-Skipping',rotation=90, fontsize=16*fontscale)

#4.0 Hz (ms)
wl_t_40 = (1.0/(2.0*4.0))*1000.0
plt.axvline(x=wl_t_40, color="black", ls="--")
plt.text(wl_t_40+txt_buff, max_num_per_bin, r'$\bf{4.0 Hz}$'+' Cycle-Skipping',rotation=90, fontsize=16*fontscale)

#6.0 Hz (ms)
#wl_t_60 = (1.0/(2.0*6.0))*1000.0
#plt.axvline(x=wl_t_60, color="black", ls="--")
#plt.text(wl_t_60+txt_buff, max_num_per_bin, '6.0 Hz Cycle-Skipping',rotation=90, fontsize=16*fontscale)

plt.show()
