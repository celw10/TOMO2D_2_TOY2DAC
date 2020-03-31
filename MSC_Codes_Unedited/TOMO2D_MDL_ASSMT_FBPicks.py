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

def data_cleanup(data, obsoffset, idx, reverse=True):
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
    for i in range(len(data)):
        tmp.append(data[i][np.argsort(data[i][:,0])])

    #Only want unique values
    tmp0 = []
    for i in range(len(data)):
        tmp0.append(tmp[i][np.unique(tmp[i][:,0], return_index=True, axis=0)[1]])

    #Interpolate first breaks over min/max traces
    tmp1 = []
    for i in range(len(tmp0)):
        f = interpolate.interp1d(tmp0[i][:,0],tmp0[i][:,1])
        tmp1.append(np.vstack((obsoffset[i][idx[0][i]:idx[1][i]],
                              f(obsoffset[i][idx[0][i]:idx[1][i]]))).T)

    #Reverse shot order
    if reverse == True:
        data_out = []
        for i in range(len(tmp1)):
            it = len(tmp1) - 1 - i
            data_out.append(tmp1[it])
    else:
        data_out = tmp1
    
    return data_out

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

def gen_toy2dac_acqi(data, obs_x, obs_z, bulk_shift, min_max, del_obs):
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

#MAKE FOURIER GATHERS

def fourier_gathers(OBS, gathers, obsoffset, obs, amp_corr, 
                    del_idx, indicesdel_min, indicesdel, min_idx, 
                    nsamples, maxfold, buff, div, sr, off_min,
                    LFC=0):
    """
    |_________________________________________________________________________________
    | PURPOSE:
    |-> Porperly organize gathers to apply the fourier transform to
    |-> Mute above the first break pick plus a buffer
    |-> Only retain shot records in the cropped model
    |-> Optionally apply the time damping constant to the data
    |-> Optionally plot the gathers with FB picks 
    |-> Remove data within a minimum offset
    |_________________________________________________________________________________
    """
    
    if offset_min > 0:
        print "WARNING: Minimum offset has been restricted to: ", offset_min
    
    #Change the minimum indicies based on cropped model
    min_idx = min_idx[::-1]
    for i in range(len(min_idx)):
        min_idx[i] = min_idx[i] + len(indicesdel_min[i]) + indicesdel[i]
    
    #Flip data (MAKE THIS OPTIONAL, NOT ALL DATASETS WILL HAVE TO BE FLIPPED) and delete
    #FLIPPED OR NOT FLIPPED AT THIS POINT WAS MESSY
    gathers = gathers[::-1]
    gathers = np.delete(gathers,del_idx,0)
    OBS = np.delete(OBS,del_idx) 
    index = np.delete(min_idx,del_idx)
    shotgeom = np.delete(obs, del_idx)
    
    #Find the maximum number of first break picks
    max_pick_num = 0
    for i in range(len(OBS)):
        if len(OBS[i][:,0]) > max_pick_num:
            max_pick_num = len(OBS[i][:,0])
    
    #Make the FFT Gathers AND optionally apply the Laplace-Fourier damping approach to gathers
    gathers_fourier = np.zeros((len(OBS),int((nsamples-1)/div),max_pick_num))
    total=0
    for num in range(len(OBS)):
        for i in range(len(OBS[num][:,0])-1):
            #Zero out the gather if within the minimum offset
            if abs(OBS[num][i,0]-shotgeom[num]) < off_min:
                gathers_fourier[num][:,i] = 0
                total+=1
            else:
                for j in range(0,int((nsamples-1)/div),1): 
                    #Apply the buffer above first break pick, optionally apply time damping to the data
                    #Apply laplace damping from first break - Gorszczyk et al., 2017 formulation
                    if (j + buff) > int(OBS[num][i,1]) and LFC>0:
                        gathers_fourier[num][j,i] = \
                        gathers[num].T[j,i+int(index[num])]*e**(-(j*sr-OBS[num][i,1]*sr)/LFC)*amp_corr[num][i] #AMP corrections from optimization routine
                    #Apply laplace damping from t=0 - Shin and Cha 2009, Brossier et al., 2009 formulation
                    elif (j + buff) > int(OBS[num][i,1]) and LFC<0:
                        gathers_fourier[num][j,i] = \
                        gathers[num].T[j,i+int(index[num])]*e**(-(j*sr)*abs(LFC))*amp_corr[num][i] #AMP corrections from optimization routine
                    #Just apply amplitude correction
                    elif (j + buff) > int(OBS[num][i,1]) and LFC==0:
                        gathers_fourier[num][j,i] = \
                        gathers[num].T[j,i+index[num]]*amp_corr[num][i] #AMP corrections from optimization routine
                        
    print "Traces removed, min offsets: : ", total
    
    return gathers_fourier, OBS, max_pick_num

#Crop FB Picks

def crop_fb_picks(OBS, gathers, obsoffset, obs, amp_corr, 
                    del_idx, indicesdel_min, indicesdel, min_idx, 
                    nsamples, maxfold, buff, div, sr, off_min,
                    LFC=0):
    """
    |_________________________________________________________________________________
    | PURPOSE:
    |-> Similar function to fourier_gathers, but without sgy data read-in
    |_________________________________________________________________________________
    """
    
    #Change the minimum indicies based on cropped model
    min_idx = min_idx[::-1]
    for i in range(len(min_idx)):
        min_idx[i] = min_idx[i] + len(indicesdel_min[i]) + indicesdel[i]
    
    #Flip data (MAKE THIS OPTIONAL, NOT ALL DATASETS WILL HAVE TO BE FLIPPED) and delete
    #FLIPPED OR NOT FLIPPED AT THIS POINT WAS MESSY
    OBS = np.delete(OBS,del_idx) 
    index = np.delete(min_idx,del_idx)
    shotgeom = np.delete(obs, del_idx)
    
    #Find the maximum number of first break picks
    max_pick_num = 0
    for i in range(len(OBS)):
        if len(OBS[i][:,0]) > max_pick_num:
            max_pick_num = len(OBS[i][:,0])
            
    return OBS

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
        print "EXTRACTING INFROMATION FOR A SINGLE SHOT, ", source
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
    
    return num_traces, nchannels

"""
Chris Williams, September 2019
-> Funciton to plot a single, high quality shot gather. 
-> Choose between IFT shot gather or segy shot gather, both are formatted differently.
-> TO ADD (wishlist): first break picks
FUNCTION DEFINITIONS:
-> data: Dataset extracted from previous
-> shot: Shot # to plot
-> mode: either "segy", indicating that the data are initially formatted as segy, or "data_modeling" indicating that the data are output from forward modeling by TOY2DAC.
-> vrange: Colorbar range, leave as a single float to multiply the min, abs(min) values of the data by a defined constant, or set to [min,max] to manually define the min/max values of the colorbar.
-> sr: sample rate in seconds, defined in the code above.
-> maxplot_t: the maximum time desired in the plotted gather
-> tick_inc: the incriment in km to plot the x tick.
-> reduction_vel: reducting velocity to produced a reduced-time plot in m/s
-> re_maxplot_t: the maximum time to show if using a reduction velocity
-> cmap: colormap from matplotlib's colormap library
-> aspect: Define the aspect ratio: leave as "auto" to automatically define a suitable aspect ratio, or manually define
-> fontsize: a sclar when defined will multiple all the text size by this number
-> size: the size of the image
-> max_off: the maximum offset avaliable in the data
"""

def plot_single_sg(data, shot, acqui, mode, vrange, sr, maxplot_t, tick_inc, cbar_pos, 
                   reduction_vel=None, red_maxplot_t=None,
                   cmap="seismic_r", aspect="auto", fontsize=1, size=[32,32], 
                   max_off=100000, savepath=None, savename=None, del_neg=True,
                   fb_pick1=None, fb_pick2 = None, fb_labels=None):
    
    #Define Figure Size
    size = size
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = size[0]
    fig_size[1] = size[1]
    plt.rcParams["figure.figsize"] = fig_size
    
    """PLOTTING MODE"""
    #SEGY or from TOY2DAC InvFT data_modeling file
    if mode=="segy":
        print "Plotting mode selected: ", mode
        print "Plotting re-formatted SEGY data."
        print
        mode_title = "Observed Data"
    elif mode=="data_modeling":
        print "Plotting mode selected: ", mode
        print "Plotting re-formatted TOY2DAC data_modeling file."
        print
        mode_title = "Modeled FWI"
    else:
        "Please select a correct plotting mode using the parameter: mode"
    
    """DEFINE DATASET"""
    if mode == "segy":
        plot_data = data[shot][0:int(maxplot_t/sr),0:len(acqui[shot][1::,1])] 
    elif mode == "data_modeling":
        print "FIRST LENGTH"
        print int(len(data[shot][:,0])-(maxplot_t/sr)) 
        print "SECOND LENGTH"
        print len(data[shot][:,0])
        plot_data = data[shot][int(len(data[shot][:,0])-(maxplot_t/sr)):len(data[shot][:,0]), 
                               0:len(acqui[shot][1::,1])] 
        print "SHAPE"
        print np.shape(plot_data)
        
        plot_data = np.flip(plot_data, axis=0)
        
    
    """DEFINE COLORBAR RANGE"""
    if type(vrange)==float:
        elev_min, elev_max = min(plot_data.flatten())*vrange, abs(min(plot_data.flatten()))*vrange
    else:
        elev_min, elev_max = vrange[0], vrange[1]
    print "Colorbar range (min/max): (", elev_min, "/", elev_max, ")"
    
    """DEFINE PLOT ASPECT RATIO"""
    if aspect != "auto":
        aspect_ratio = aspect
    else:
        aspect_ratio = aspect
        
    print "SHAPE OF THE PLOTTED DATA"
    print np.shape(plot_data)
    print
    
    """CONVERT TO REDUCED TIME"""
    ### SEG-Y MODE ###
    if reduction_vel != None and mode == "segy":
        #Plotting specific definitions
        y_axis_label = "Time - Offset/"+str(reduction_vel/1000.0)+" [s]"
        #Extract reduced time nodes
        tnodes=[]
        print "Converting plot to reduced time using a reduction velocitiy of: ", reduction_vel
        for i in range(len(acqui[shot][1::,1])-1):
            x = acqui[shot][i+1,1] - acqui[shot][0,1]
            for j in range(int(maxplot_t/sr)):
                tnodes.append(int((j*sr-(abs(x)/reduction_vel))/sr))
        #Extract min/max fvalues
        max_tnode, min_tnode = max(tnodes), min(tnodes)
        max_redtime, min_redtime = max_tnode*sr, min_tnode*sr
        #Fill reduced-time shot gather
        red_plot_data = np.zeros((int(max_tnode-min_tnode),len(acqui[shot][1::,1])))
        for i in range(len(acqui[shot][1::,1])-1):
            for j in range(int(maxplot_t/sr)):
                red_plot_data[tnodes[int(i*maxplot_t/sr+j)]-min_tnode-1,i] = plot_data[j,i]
        #Tidy up time-reduced shot gather
        max_minredtime = maxplot_t - (max_off/reduction_vel)
        if max_minredtime < 0:
            print "Please use a larger reduction velocity, negatives encountered."
        else:
            data_yextent = abs(min_redtime)
            data_textent = 0
        figure = red_plot_data[int(data_yextent/sr):int((red_maxplot_t+data_yextent)/sr), 0:len(acqui[shot][:,0])-1]
    
    ### DATA_MODELING MODE ### NOT SURE I NEED TWO SEPARATE ANYMORE EITHER!
    elif reduction_vel != None and mode == "data_modeling":
        #Plotting specific definitions
        y_axis_label = "Time - Offset/"+str(reduction_vel/1000.0)+" [s]"
        #Extract reduced time nodes
        tnodes=[]
        print "Converting plot to reduced time using a reduction velocitiy of: ", reduction_vel
        for i in range(len(acqui[shot][1::,1])-1):
            x = acqui[shot][i+1,1] - acqui[shot][0,1]
            for j in range(int(maxplot_t/sr)):
                tnodes.append(int((j*sr-(abs(x)/reduction_vel))/sr))
        #Extract min/max fvalues
        max_tnode, min_tnode = max(tnodes), min(tnodes)
        max_redtime, min_redtime = max_tnode*sr, min_tnode*sr
        #Fill reduced-time shot gather
        red_plot_data = np.zeros((int(max_tnode-min_tnode),len(acqui[shot][1::,1])))
        print np.shape(red_plot_data)
        for i in range(len(acqui[shot][1::,1])-1):
            for j in range(int(maxplot_t/sr)):
                red_plot_data[tnodes[int(i*maxplot_t/sr+j)]-min_tnode-1,i] = plot_data[j,i]
        #Tidy up time-reduced shot gather
        max_minredtime = maxplot_t - (max_off/reduction_vel)
        if max_minredtime < 0:
            print "Please use a larger reduction velocity, negatives encountered."
        else:
            data_yextent = abs(min_redtime)
            data_textent = 0
        print "TERM 1"
        print int(data_yextent/sr)
        print "TERM 2"
        print int((red_maxplot_t+data_yextent)/sr)
        figure = red_plot_data[int(data_yextent/sr):int((red_maxplot_t+data_yextent)/sr), 0:len(acqui[shot][:,0])-1]
        #figure = red_plot_data[int(max_tnode-min_tnode)-int((red_maxplot_t+data_yextent)/sr):
        #                       int(max_tnode-min_tnode)-int(data_yextent/sr),
        #                       0:len(acqui[shot][:,0])-1]
        print "SHAPE OF THE ALTERED FIGURE"
        print np.shape(figure)
        print
        
    else:
        y_axis_label = 'Time (s) '
        figure = plot_data
    
    """First Break Picks to Reduced Time"""
    
    if fb_pick1 != None:
        #Extract reduced time nodes
        tnodes_fb1 = np.zeros((len(fb_pick1[shot][:,0]),2))
        #I need to fill out the x somehow
        print "Converting first break picks to reduced time using a reduction velocitiy of: ", reduction_vel
        for i in range(len(fb_pick1[shot][:,0])):
            x = (fb_pick1[shot][i,0]-min_max[0])*1000.0 - acqui[shot][0,1]
            tnodes_fb1[i,1] = ((fb_pick1[shot][i,1])*sr-(abs(x)/reduction_vel))
    
    if fb_pick2 != None:
        #Extract reduced time nodes
        tnodes_fb2 = np.zeros((len(fb_pick2[shot][:,0]),2))
        #I need to fill out the x somehow
        print "Converting first break picks to reduced time using a reduction velocitiy of: ", reduction_vel
        for i in range(len(fb_pick2[shot][:,0])):
            x = (fb_pick2[shot][i,0]-min_max[0])*1000.0 - acqui[shot][0,1]
            tnodes_fb2[i,1] = ((fb_pick2[shot][i,1])*sr-(abs(x)/reduction_vel))

    #We need to figure out the relative shift between first break picks
    #The first group of first break picks are the reference, of which the shot gather whas generated
    if fb_pick1 != None and fb_pick2 != None:
        i, j, flag = 0, 0, "missing"
        while i < len(fb_pick1[shot][:,0]):
            while j < len(fb_pick2[shot][:,0]):
                if fb_pick2[shot][j,0] == fb_pick1[shot][i,0]:
                    out_shift, flag = i-j, "found"
                    i, j = len(fb_pick1[shot][:,0]), len(fb_pick2[shot][:,0])
                else:
                    j+=1
            if flag == "found":
                print "Shift required relative to reference FB Picks: ", out_shift
            else:
                i+=1
        #X-nodes for pick2
        for i in range(out_shift, len(fb_pick2[shot][:,0])+out_shift,1):
            tnodes_fb2[i-out_shift,0] = i
        #Optionally delete negative node values, not extending the reduced time plot
        if del_neg == True:
            neg_idx = []
            for i in range(len(tnodes_fb2[:,0])):
                if tnodes_fb2[i,0] < 0:
                    neg_idx.append(i)
            print neg_idx
            tnodes_fb2 = np.delete(tnodes_fb2, neg_idx, axis=0)
    #Just a range of nodes for the reference FB pick set
    if fb_pick1 != None:
        for i in range(len(fb_pick1[shot][:,0])):
            tnodes_fb1[i,0] = i

    
    """DEFINE PLOT EXTENT""" # I DON"T NEED THIS DIFFERENCE FOR segy and data_modeling ANYMORE
    if reduction_vel != None and mode == "segy":
        extent=[0,len(acqui[shot][1::,1]),red_maxplot_t,0]
    elif reduction_vel != None and mode == "data_modeling":
        extent=[0,len(acqui[shot][1::,1]),red_maxplot_t,0]
    elif mode == "segy":
        extent=[0,len(acqui[shot][1::,1]),maxplot_t,0] 
    elif mode == "data_modeling":
        extent=[0,len(acqui[shot][1::,1]),maxplot_t,0]
    
    #Plot the shot gather
    fig = plt.figure() 
    ax = fig.add_subplot(111)
    im = ax.imshow(figure, aspect=aspect, vmin=elev_min, vmax=elev_max, cmap=cmap, extent=extent)
    if fb_pick1!=None:
        ax.plot(tnodes_fb1[:,0],tnodes_fb1[:,1], "b-", linewidth=3.0, label=fb_labels[0])
    if fb_pick2!=None:
        ax.plot(tnodes_fb2[:,0],tnodes_fb2[:,1], "r--", linewidth=3.0, label=fb_labels[1])
    
    """CONFIGURE X-AXIS TICKS"""
    #Define find nearest function for x-tick construction
    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx
    #Format x-ticks, there is no zero offset
    all_ticks = range(-max_off, max_off, tick_inc)
    min_tick_idx = find_nearest(all_ticks, int(min((acqui[shot][1::,1]-acqui[shot][0,1])/1000.)))
    if all_ticks[min_tick_idx] < int(min((acqui[shot][1::,1]-acqui[shot][0,1])/1000.)):
        min_tick = all_ticks[min_tick_idx+1]
    else:
        min_tick = all_ticks[min_tick_idx]
    #Find the ticks and best indicies for the x-axis
    x_tick_index = []
    x_tick_values = range(min_tick, int(max((acqui[shot][1::,1]-acqui[shot][0,1])/1000.)), tick_inc)
    for i in range(len(x_tick_values)):
        x_tick_index.append(find_nearest(((acqui[shot][1::,1]-acqui[shot][0,1])/1000.), x_tick_values[i]))
    plt.xticks(x_tick_index,x_tick_values)
    
    #Titles and plot labels 
    #plt.title(mode_title+', OBS '+str(shot+1),fontweight="bold", size=20.0*fontsize, x=0.10, y=0.95) 
    plt.ylabel(y_axis_label, fontsize = 18.0*fontsize) # Y label
    ax.set_xlabel('Offset [km]', fontsize = 18.0*fontsize) # X label
    ax.xaxis.set_label_position('top') 
    ax.xaxis.tick_top()
    plt.tick_params(axis='both', which='major', labelsize=16.0*fontsize)
    plt.tick_params(axis='both', which='minor', labelsize=16.0*fontsize)
    plt.grid()
    
    #fig.patches.extend([plt.Rectangle((0.63, 0.16),0.32, 0.036,
    #                              fill=True, color='k', alpha=0.5, zorder=1000,
    #                              transform=fig.transFigure, figure=fig)])

    #Colorbar
    ax_cbar = fig.add_axes([cbar_pos[0], cbar_pos[1], 0.22, 0.03]) #0.68, 0.263 OBS 12
    cbar = fig.colorbar(im, cax=ax_cbar, ax=ax, orientation='horizontal', format='%.0e')
    cbar.ax.tick_params(labelsize=12.0*fontsize)
    cbar.set_label("Amplitude", fontsize=14.0*fontsize)
    
    #Invert axis for toy2dac file
    #if mode == "data_modeling":
    #    ax.invert_yaxis()
    
    #Plot legend
    if fb_labels!=None:
        ax.legend(loc="lower left", fontsize=16*fontsize)
        
    if savepath != None and savename != None:
        plt.gcf().savefig(savepath+savename+".png", bbox_inches = 'tight', pad_inches = 0.1)
    
    plt.show()
    
###TEST WHAT I'VE DONE SO FAR###

import numpy as np
import segyio #This doesn't come standard with python
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.colors as mcolors
import struct

from scipy import interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sympy import N
from math import e

#############
### FLAGS ###
#############

#Number of first break picks to plot (maximum 2, minimum 1), 0 to only plot "final" first breaks
num_fb_picks = 1

#############
### FILES ###
#############

#Full path to tx.in file (Final FB Pick File)
f = "/home/chris/Documents/1_data/EasternMed/Chris_Williams/first_breaks/pick9.3a"
picks = np.genfromtxt(f) 

#Full path to FB pick file 1 (for plotting)
f = "/home/chris/Documents/1_data/EasternMed/Chris_Williams/first_breaks/pick8a.in"
picks1 = np.genfromtxt(f) 

#Full path to FB pick file 2 (for plotting)
f = "/home/chris/Documents/1_data/EasternMed/Chris_Williams/first_breaks/pick9a.in"
picks2 = np.genfromtxt(f) 

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
segydir = '/home/chris/0_GlobeClaritas/0_Projects/Ref_EMed_Line1/COMMON/DATA/5_FinalFlow_Sept30th/3_FWIData_RMSNorm/'
name_start = ''
name_end = '_VERT_Final_FWI.sgy'
fileloc = [segydir, name_start, name_end]

#Amplitude correction file generated externally (usually saved with the data_modeling file for the starting model)
f = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/6_ForwardModeling/006_Aug02_StartMdl66/AmpCorr_Sept072019.txt"
amp_corr = np.genfromtxt(f)

#################
### VARIABLES ###
#################

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
offset_max = 60
offset_min = 0
bulk_shift = 0 # Not really implimented
min_max = [30, 220]

#Indicies of OBS that have already been cropped out, blank if none: PYTHON INDEXING, OBS number - 1
del_obs = [8, 10]

#Define variables for fourier transform
buff = 200 #retain this many time samples above the first break
div = 2 #divide the total number of time samples by this number (saves memory)
LFC = 0.0 #the Laplace-Fourier Constant, set to zero if it isn't to be used, negative constant applys damping from t=0

#Saving reduced-time plot
savepath = None
savename = None

################
### GEOMETRY ###
################

#Read the first-beak picks    
OBS, pick_numbers = read_fbpicks(picks, obs, phase, unc, tt_window, sr)
if num_fb_picks == 1 or num_fb_picks == 2:
    #Optional FB picks/forward modeled data #1
    OBS_1, pick_numbers = read_fbpicks(picks1, obs, phase, unc, tt_window, sr)
if num_fb_picks == 2:
    #Optional FB picks/forward modeled data #2
    OBS_2, pick_numbers = read_fbpicks(picks2, obs, phase, unc, tt_window, sr)
    
#################
### SEGY DATA ###
#################

#Extract obsoffsets from the header and the gathers from sgy file
gathers, obsoffsets = extract_segy(fileloc, len(obs), maxfold, nsamples)
#Indexes for first break picks in segy data
min_idx, max_idx = indices_from_segy(OBS, obsoffsets)
if num_fb_picks == 1 or num_fb_picks == 2:
    #Optional FB picks/forward modeled data #1
    min_idx1, max_idx1 = indices_from_segy(OBS_1, obsoffsets)
if num_fb_picks == 2:
    #Optional FB picks/forward modeled data #2
    min_idx2, max_idx2 = indices_from_segy(OBS_2, obsoffsets)

################
### CLEAN UP ###
################

#Ensure uniqueness and interpolate gaps in dataset.
OBS = data_cleanup(OBS, obsoffsets, [min_idx, max_idx])
if num_fb_picks == 1 or num_fb_picks == 2:
    #Optional FB picks/forward modeled data #1
    OBS_1 = data_cleanup(OBS_1, obsoffsets, [min_idx1, max_idx1])
if num_fb_picks == 2:
    #Optional FB picks/forward modeled data #2
    OBS_2 = data_cleanup(OBS_2, obsoffsets, [min_idx2, max_idx2])

####################
### MODEL EXTENT ###
####################

#Set the maximum and minimum values for model cropping
OBS, indicesdel = max_offset(OBS, offset_max)
#Shift the data with the OBS model locations and a bulk shift
OBS = shift_data(OBS, obs, bulk_shift)
#Set the maximum and minimum values for model cropping
OBS, indicesdel_crop = crop_model(OBS, min_max)

#Optional FB picks/forward modeled data #1
if num_fb_picks == 1 or num_fb_picks == 2:
    #Set the maximum and minimum values for model cropping
    OBS_1, indicesdel_1 = max_offset(OBS_1, offset_max)
    #Shift the data with the OBS model locations and a bulk shift
    OBS_1 = shift_data(OBS_1, obs, bulk_shift)
    #Set the maximum and minimum values for model cropping
    OBS_1, indicesdel_crop_1 = crop_model(OBS_1, min_max)

#Optional FB picks/forward modeled data #2
if num_fb_picks == 2:
    #Set the maximum and minimum values for model cropping
    OBS_2, indicesdel_2 = max_offset(OBS_2, offset_max)
    #Shift the data with the OBS model locations and a bulk shift
    OBS_2 = shift_data(OBS_2, obs, bulk_shift)
    #Set the maximum and minimum values for model cropping
    OBS_2, indicesdel_crop_2 = crop_model(OBS_2, min_max)

#########################
### OUTPUT ACQUI FILE ###
#########################

#Output final version of acqui file
acqui, del_idx = gen_toy2dac_acqi(OBS, obs, obs_z, bulk_shift, min_max, del_obs)
if num_fb_picks == 1 or num_fb_picks == 2:
    #Optional FB picks/forward modeled data #1
    acqui_1, del_idx_1 = gen_toy2dac_acqi(OBS_1, obs, obs_z, bulk_shift, min_max, del_obs)
if num_fb_picks == 2:
    #Optional FB picks/forward modeled data #2
    acqui_2, del_idx_2 = gen_toy2dac_acqi(OBS_2, obs, obs_z, bulk_shift, min_max, del_obs)

######################
### OUTPUT GATHERS ###
######################

#Prep gathers for fourier transform
RD_SG, OBS, max_pick_num = fourier_gathers(OBS, gathers, obsoffsets, obs, amp_corr,
                                             del_idx, indicesdel_crop, indicesdel, min_idx,
                                             nsamples, maxfold, buff, div, sr, offset_min,
                                             LFC=LFC)
#Delete proper indicies for OBS 1
OBS_1 = np.delete(OBS_1, del_idx_1)

#Delete proper indicies for OBS 2
#OBS_2 = np.delete(OBS_2, del_idx_2)

#Save some memory 
del obsoffsets, gathers

#################################
### PRODUCE REDUCED TIME PLOT ###
#################################

#Plotting Definitions
shot = 12 #Python Indexing, index+1
vrange = [-5e-5, 5e-5] #Optional, dependent on amp corr file
cbar_pos = [0.68, 0.288]#[0.68, 0.35]##60km #First left and right, second up and down

#Labels
fb_labels=["Pick 8", "Pick 9"]

#Plot the shot gather
plot_single_sg(RD_SG, shot, acqui, "segy", vrange, sr, int(tmax/div), 10, cbar_pos,
                   reduction_vel=6500, red_maxplot_t=10,
                   cmap="gray", aspect=20, fontsize=2, size=[48,24], 
                   max_off=100000, savepath=savepath, savename=savename, 
                   fb_pick1=OBS, fb_pick2=OBS_1, fb_labels=fb_labels)