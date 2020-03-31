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
    taper_full = np.hanning(buff*2)
    taper_half = taper_full[0:buff]
    
    for num in range(len(OBS)):
        for i in range(len(OBS[num][:,0])-1):
            #Zero out the gather if within the minimum offset
            if abs(OBS[num][i,0]-shotgeom[num]) < off_min:
                gathers_fourier[num][:,i] = 0
                total+=1
            else:
                for j in range(0,int((nsamples-1)/div),1): 
                    #Apply laplace damping from first break - Gorszczyk et al., 2017 formulation
                    if (j + buff) > int(OBS[num][i,1]) and LFC>0:
                        #Taper the buffer zone with a hanning (no hard edges)
                        if j < int(OBS[num][i,1]):
                            gathers_fourier[num][j,i] = \
                            gathers[num].T[j,i+int(index[num])]*(taper_half[buff-(int(OBS[num][i,1])-j)])*amp_corr[num][i] 
                            #Amplitude corrections from AMPLITUDE_PROCESS
                        #Apply damping after the first break outside the buffer zone
                        else:
                            gathers_fourier[num][j,i] = \
                            gathers[num].T[j,i+int(index[num])]*e**(-(j*sr-OBS[num][i,1]*sr)/LFC)*amp_corr[num][i] 
                            #Amplitude corrections from AMPLITUDE_PROCESS
                            
                    #Apply laplace damping from t=0 - Shin and Cha 2009
                    elif (j + buff) > int(OBS[num][i,1]) and LFC<0:
                        #Taper the buffer zone with a hanning (no hard edges)
                        if j < int(OBS[num][i,1]):
                            gathers_fourier[num][j,i] = \
                            gathers[num].T[j,i+int(index[num])]*(taper_half[buff-(int(OBS[num][i,1])-j)])*amp_corr[num][i] 
                            #Amplitude corrections from AMPLITUDE_PROCESS
                        #Apply damping after the first break outside the buffer zone
                        else:
                            gathers_fourier[num][j,i] = \
                            gathers[num].T[j,i+int(index[num])]*e**(-(j*sr)*abs(LFC))*amp_corr[num][i] 
                            #AMP corrections from optimization routine
                            
                    #Just apply amplitude correction
                    elif (j + buff) > int(OBS[num][i,1]) and LFC==0:
                        #Taper the buffer zone with a hanning (no hard edges)
                        if j < int(OBS[num][i,1]):
                            gathers_fourier[num][j,i] = \
                            gathers[num].T[j,i+int(index[num])]*(taper_half[buff-(int(OBS[num][i,1])-j)])*amp_corr[num][i] 
                            #Amplitude corrections from AMPLITUDE_PROCESS
                        #First break outside the buffer zone
                        else:       
                            gathers_fourier[num][j,i] = \
                            gathers[num].T[j,i+index[num]]*amp_corr[num][i] 
                        #AMP corrections from optimization routine
        
    return gathers_fourier, OBS, max_pick_num

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

def taper_FD_data(FD_Comp_SG, taper_design, nsources, nfreq, nrelfreq, nchannels, sources=all):
    """
    -> Applies a taper to the frequency domain data to only retain the relevant frequencies
    -> Tapering the FD data prior to IFT results in less noise generated through the IFT
    Params:
    -> FD_Comp_SG = is the data output form the data_modeling restructuring completed previously in
    restructure_toy2dac_data
    -> taper_design = a list of three values containing the number of discrete frequency points defining the shape
    of a modified hanning window. taper_design[0] containts the width of the low taper. 
    taper_design[1] contains the width of the 100% passed portion of the modified hanning window, and taper_design[2]
    specifies the width of the high taper.
    -> nsources = number of shots
    -> nfreq = total number of discrete frequencies required for the inverse fourier transform (IFT)
    -> nrelfreq = number of relevant frequencies, i.e. discrete frequencies that were actually forward modeled in 
    toy2dac. 
    -> nchannels = maximum fold of the data, should be output form setup_geom previously
    -> sources = number of sources to restructure, all would restructure all sources, or enter the source index to 
    only restrucutre that source. 
    """

    ######################
    ### TAPER F-D DATA ###
    ######################

    #Low Window
    window_l = np.hanning(taper_design[0]*2)
    taper_l = window_l[0:taper_design[0]]
    #High Window
    window_h = np.hanning(taper_design[2]*2)
    taper_h = window_h[taper_design[2]:taper_design[2]*2]
    #Combine three sections
    taper = np.zeros((int(nfreq)))
    taper[0:taper_design[0]] = taper_l
    taper[taper_design[0]:taper_design[0]+taper_design[1]] = 1.0
    taper[taper_design[0]+taper_design[1]:taper_design[0]+taper_design[1]+taper_design[2]] = taper_h

    #Apply the taper
    print "APPLYING TAPER..."
    if sources == all:
        FD_Comp_SG_taper = np.zeros((nsources, int(nfreq), nchannels),dtype=complex)
        for s in range(nsources):
            for f in range(int(relfreq)):
                FD_Comp_SG_taper[s,f,:]=FD_Comp_SG[s,f,:]*taper[f]
    else:
        FD_Comp_SG_taper = np.zeros((int(nfreq), nchannels),dtype=complex)
        for f in range(int(relfreq)):
            FD_Comp_SG_taper[f,:]=FD_Comp_SG[f,:]*taper[f]
            
    print "TAPER APPLIED"
    print
    
    return FD_Comp_SG_taper

def inv_fft(FD_Comp_SG_taper, nsources, nfreq, nchannels, sources=all):
    """
    -> Apply the inverse fourier transform
    Params:
    -> FD_Comp_SG_taper = complex matrix containg the data_modeling file with the taper applied
    -> nsources = number of shots
    -> nfreq = total number of discrete frequencies required for the inverse fourier transform (IFT)
    -> nchannels = maximum fold of the data, should be output form setup_geom previously
    -> sources = number of sources to restructure, all would restructure all sources, or enter the source index to 
    only restrucutre that source. 
    README:
    -> Frequency spacing governs your maximum time: tmax=1/df
    -> Maximum frequency governs your time sampling: dt=1/fmax
    """
    
    ##################
    ### APPLY IFFT ###
    ##################

    #Source: https://docs.scipy.org/doc/numpy/reference/generated/numpy.fft.ifft.html
    print "BEGINNING INVERSE FOURIER TRANSFORM..."
    if sources == all:
        TD_SG = np.zeros((nsources, int(nfreq), nchannels))
        for s in range(nsources):
            TD_SG[s] = np.fft.ifft(FD_Comp_SG_taper[s],axis=-2)
    else:
        TD_SG = np.fft.ifft(FD_Comp_SG_taper,axis=-2)
    print "COMPLETED INVERSE FOURIER TRANSFORM"
    print
    
    return TD_SG.real

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

def plot_single_sg(data, shot, acqui, mode, vrange, sr, maxplot_t, tick_inc, cbar_pos, title_pos,
                   reduction_vel=None, red_maxplot_t=None,
                   cmap="seismic_r", aspect="auto", fontsize=1, size=[32,32], 
                   max_off=100000, savepath=None, savename=None):
    
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
    elif mode=="data_modeling":
        print "Plotting mode selected: ", mode
        print "Plotting re-formatted TOY2DAC data_modeling file."
        print
    elif mode=="single_data_modeling":
        print "Plotting mode selected: ", mode
        print "Plotting a single source TOY2DAC data_modeling file."
        print
    else:
        "Please select a correct plotting mode using the parameter: mode"
    
    """DEFINE DATASET"""
    if mode == "segy":
        plot_data = data[shot][0:int(maxplot_t/sr),0:len(acqui[shot][1::,1])] 
    elif mode == "data_modeling":
        plot_data = data[shot][int(len(data[shot][:,0])-(maxplot_t/sr)):len(data[shot][:,0]), 
                               0:len(acqui[shot][1::,1])] 
        plot_data = np.flip(plot_data, axis=0)
    elif mode == "single_data_modeling":
        acqui = acqui[0]
        plot_data = data[int(len(data[:,0])-(maxplot_t/sr)):len(data[:,0]), 
                               0:len(acqui[1::,1])]
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
    
    ### DATA_MODELING MODE, MULTIPLE SHOTS ### 
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
        figure = red_plot_data[int(data_yextent/sr):int((red_maxplot_t+data_yextent)/sr), 0:len(acqui[shot][:,0])-1]
        
    ### DATA MODELING MODE, SINGLE SHOT ###
    elif reduction_vel != None and mode == "single_data_modeling":
        #Plotting specific definitions
        y_axis_label = "Time - Offset/"+str(reduction_vel/1000.0)+" [s]"
        #Extract reduced time nodes
        tnodes=[]
        print "Converting plot to reduced time using a reduction velocitiy of: ", reduction_vel
        for i in range(len(acqui[1::,1])-1):
            x = acqui[i+1,1] - acqui[0,1]
            for j in range(int(maxplot_t/sr)):
                tnodes.append(int((j*sr-(abs(x)/reduction_vel))/sr))
        #Extract min/max fvalues
        max_tnode, min_tnode = max(tnodes), min(tnodes)
        max_redtime, min_redtime = max_tnode*sr, min_tnode*sr
        #Fill reduced-time shot gather
        red_plot_data = np.zeros((int(max_tnode-min_tnode),len(acqui[1::,1])))
        print np.shape(red_plot_data)
        for i in range(len(acqui[1::,1])-1):
            for j in range(int(maxplot_t/sr)):
                red_plot_data[tnodes[int(i*maxplot_t/sr+j)]-min_tnode-1,i] = plot_data[j,i]
        #Tidy up time-reduced shot gather
        max_minredtime = maxplot_t - (max_off/reduction_vel)
        if max_minredtime < 0:
            print "Please use a larger reduction velocity, negatives encountered."
        else:
            data_yextent = abs(min_redtime)
            data_textent = 0
        figure = red_plot_data[int(data_yextent/sr):int((red_maxplot_t+data_yextent)/sr), 0:len(acqui[:,0])-1]
    
    else:
        y_axis_label = 'Time (s) '
        figure = plot_data

    
    """DEFINE PLOT EXTENT""" 
    
    if reduction_vel != None and mode == "segy":
        extent=[0,len(acqui[shot][1::,1]),red_maxplot_t,0]
    elif reduction_vel != None and mode == "data_modeling":
        extent=[0,len(acqui[shot][1::,1]),red_maxplot_t,0]
    elif reduction_vel != None and mode == "single_data_modeling":
        extent=[0,len(acqui[1::,1]),red_maxplot_t,0]
    elif mode == "segy":
        extent=[0,len(acqui[shot][1::,1]),maxplot_t,0] 
    elif mode == "data_modeling":
        extent=[0,len(acqui[shot][1::,1]),maxplot_t,0] 
    elif mode == "single_data_modeling":
        extent=[0,len(acqui[1::,1]),maxplot_t,0] 
    
    """PLOT THE SHOT GATHER"""
    
    fig = plt.figure() 
    ax = fig.add_subplot(111)
    im = ax.imshow(figure, aspect=aspect, vmin=elev_min, vmax=elev_max, cmap=cmap, extent=extent)
    
    """CONFIGURE X-AXIS TICKS"""
    
    #Define find nearest function for x-tick construction
    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx
    #Define all as acqui to avoid re-definitions
    if mode == "segy" or mode == "data_modeling":
        acqui = acqui[shot]
    #Format x-ticks, there is no zero offset
    all_ticks = range(-max_off, max_off, tick_inc)
    min_tick_idx = find_nearest(all_ticks, int(min((acqui[1::,1]-acqui[0,1])/1000.)))
    if all_ticks[min_tick_idx] < int(min((acqui[1::,1]-acqui[0,1])/1000.)):
        min_tick = all_ticks[min_tick_idx+1]
    else:
        min_tick = all_ticks[min_tick_idx]
    #Find the ticks and best indicies for the x-axis
    x_tick_index = []
    x_tick_values = range(min_tick, int(max((acqui[1::,1]-acqui[0,1])/1000.)), tick_inc)
    for i in range(len(x_tick_values)):
        x_tick_index.append(find_nearest(((acqui[1::,1]-acqui[0,1])/1000.), x_tick_values[i]))
    plt.xticks(x_tick_index,x_tick_values)
    
    """PLOT LABELS AND TITLES"""
    
    #plt.title('Modeled Data',fontweight="bold", size=24.0*fontsize, x=0.10+title_pos[0], #+str(shot+1)
    #          y=0.95+title_pos[1]) 
    plt.ylabel(y_axis_label, fontsize = 18.0*fontsize) # Y label
    ax.set_xlabel('Offset [km]', fontsize = 18.0*fontsize) # X label
    ax.xaxis.set_label_position('top') 
    ax.xaxis.tick_top()
    plt.tick_params(axis='both', which='major', labelsize=16.0*fontsize)
    plt.tick_params(axis='both', which='minor', labelsize=16.0*fontsize)
    plt.grid()

    """COLORBAR"""
    
    ax_cbar = fig.add_axes([cbar_pos[0], cbar_pos[1], 0.28, 0.05]) #0.68, 0.263 OBS 12
    cbar = fig.colorbar(im, cax=ax_cbar, ax=ax, orientation='horizontal', format='%.0e')
    cbar.ax.tick_params(labelsize=14.0*fontsize)
    cbar.set_label("Amplitude", fontsize=16.0*fontsize)
        
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
from sympy import N
from math import e

#############
### FLAGS ###
#############

#Define what datasets to read, 0=real and modeled, 1=real only, 2=modeled only
dataset = 1

#If more than one shot gather leave as zero, if not set to 1
data_num = 0

#############
### FILES ###
#############

###############
"""REAL DATA"""
###############

#Full path to tx.in file
f = "/home/chris/Documents/1_data/EasternMed/Chris_Williams/first_breaks/pick9a.final.zerophase"
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
segydir = '/home/chris/0_GlobeClaritas/0_Projects/Ref_EMed_Line1/COMMON/DATA/5_FinalFlow_Sept30th/4_FWIData_FilterforImages/'
name_start = ''
name_end = '_VERT_Filter_FWI.sgy'
fileloc = [segydir, name_start, name_end]

#Amplitude correction file generated externally (usually saved with the data_modeling file for the starting model)
f = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/6_ForwardModeling/009_Oct29_StartMdl110/AmpCorr_QP050_sdNov12.txt"
amp_corr = np.genfromtxt(f)

##################
"""TOY2DAC DATA"""
##################

#Open toy2dac data_modeling file
dirr = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/008_2020RealDataTests/043_MiscTestsFGProg1_SmthCnstL1/8_fwdmdl/" #Path to root directory
model = "dm_15Hz_mdl043"
filename = dirr + model #Name of binary file
with open(filename, 'rb') as f:
    data = np.fromfile(f, dtype=np.float32)

#Maximum TOY2DAC geometry file
f = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/000_DATA/2_RealData/17_110StartMdl_THEEND/01_acqui_t0/acqui_100" 

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
bulk_shift = 0 # Not really implimented
min_max = [30, 220] #This is in the unit of the tx.in file

#Indicies of OBS that have already been cropped out, blank if none: PYTHON INDEXING, OBS number - 1
del_obs = [11,13]#[0,1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,19,20]#

#Define variables for fourier transform
buff = 100 #retain this many time samples above the first break
div = 1.5 #divide the total number of time samples by this number (saves memory)
LFC = 0.0 #the Laplace-Fourier Constant, set to zero if it isn't to be used, negative constant applys damping from t=0

#Plot number
plot = False #set as False if to plots are desired.

##################
"""TOY2DAC DATA"""
##################

#Data-related
nchannels = 853 #Number of channels, obsolete if setup_flag != 0 
nsources = 16 #Number of sources 
nfreq = 50000 #Number of discrete freqeuncy groups
relfreq = 3000 #Number of "relevant" frequency groups, i.e. sub 30 Hz
maxfreq = 250 #Maximum frequency in forward modeled data

#Algorithm-related
#Taper design in terms of # of discrete frequencies
taper_design = [400, 1000, 1600] #Save taper design into compact list, [low taper, 100% pass length, high taper]/5

################
### GEOMETRY ###
################

"""REAL DATA"""
#Read the first-beak picks    
OBS, pick_numbers = read_fbpicks(picks, obs, phase, unc, tt_window, sr)

if dataset == 0 or dataset == 2:
    """TOY2DAC DATA"""
    #Setup geometry from toy2dac acqui file, setup flag=1 if non constant shots/OBS
    if data_num == 0:
        num_traces, nchannels = setup_geom(f, setup_flag=1, sources=all)
    if data_num == 1:
        num_traces=nchannels
    
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

if dataset == 0 or  dataset == 2:
    #Restructure toy2dac data_modeling file into matrix format
    if data_num == 0:
        FD_Comp_SG = restructure_toy2dac_data(data, nsources, nfreq, relfreq, nchannels, num_traces, sources=all)
    elif data_num == 1:
        FD_Comp_SG = restructure_toy2dac_data(data, nsources, nfreq, relfreq, nchannels, num_traces, sources="one")        
    #Save some memory 
    del data

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

"""REAL DATA"""
#Output final version of acqui file
acqui, del_idx = gen_toy2dac_acqi(OBS, obs, obs_z, bulk_shift, min_max, del_obs)

#Use the del_obs values to remove everything other than the OBS of interest 
#if data_num == 1:
#    acqui = np.delete(acqui,del_obs)
    
if dataset == 0 or dataset == 1:
    
    ######################
    ### OUTPUT GATHERS ###
    ######################
    
    """REAL DATA"""
    #Prep gathers for fourier transform
    RD_SG, OBS, max_pick_num = fourier_gathers(OBS, gathers, obsoffsets, obs, amp_corr,
                                                 del_idx, indicesdel_crop, indicesdel, min_idx,
                                                 nsamples, maxfold, buff, div, sr, offset_min,
                                                 LFC=LFC)
    if dataset == 1:
        #Save some memory 
        del obsoffsets, gathers
    
if dataset == 0 or dataset == 2:
    
    #Save some memory 
    del obsoffsets, gathers

    """TOY2DAC DATA DATA"""

    #Filter frequency domain data
    if data_num == 0:
        #Taper the FD data
        FD_Comp_SG_taper = taper_FD_data(FD_Comp_SG, taper_design, nsources, nfreq, relfreq, nchannels, sources=all)
        del FD_Comp_SG
        #Inverse fourier transform to time domain
        FM_SG = inv_fft(FD_Comp_SG_taper, nsources, nfreq, nchannels, sources=all)
    else:
        #Taper the FD data
        FD_Comp_SG_taper = taper_FD_data(FD_Comp_SG, taper_design, nsources, nfreq, relfreq, nchannels, sources=0)
        del FD_Comp_SG
        #Inverse fourier transform to time domain
        FM_SG = inv_fft(FD_Comp_SG_taper, nsources, nfreq, nchannels, sources=0)

    #Save some memory
    del FD_Comp_SG_taper
    
#Definitions (not previously defined)
#SOMETHING WRONG WITH OBS 6? (python 5), not filling out the entirity of the OBS geometry 
#OBS 9 has some high amplitude far offset noise
#AFTER LUNCH ADD REDUCTION TIME AND COLORBAR, TEST WITH FD DATA

from mpl_toolkits.axes_grid1 import make_axes_locatable

"""
When I'm making these final figures, I'm going to keep the maximum offset to 60 km to keep it from getting too ugly.
"""

data = RD_SG#FM_SG #Whatever the dataset is defined as
mode = "segy" #"data_modeling", "segy", "single_data_modeling
shot = 15 #Python Indexing
vrange = [-5e-4, 5e-4] #The range of the colorbar
maxplot_t = 25 #The maximum time for the plot
reduction_vel = 6500 #Reduction velocity to use, set to None if undesired
red_maxplot_t = 10 #Maximum plot time for the reduced-velocity plot
aspect = 25 #The height:length ratio of the plot, higher stretches the plot out
tick_inc = 10 #Incriment (in km) for the x-axis ticks
cbar_pos = [0.62, 0.12475] # OBS 10 [0.68, 0.334] asp 30,[0.675, 0.275] asp 25, [0.675, 0.321] asp 20
title_pos = [0.03,-0.08] # OBS 10 [0.005,-0.02] OBS 4 [0.02,-0.02] OBS 13[0.01,-0.02], [0.022,-0.02]

#To automatically save the figure
savepath = None
savename = None

#Plot the shot gather
plot_single_sg(data, shot, acqui, mode, vrange, sr, maxplot_t, tick_inc, cbar_pos, title_pos,
                   reduction_vel=reduction_vel, red_maxplot_t=red_maxplot_t,
                   cmap="gray", aspect="auto", fontsize=2.5, size=[48,24], 
                   max_off=100000, savepath=savepath, savename=savename)