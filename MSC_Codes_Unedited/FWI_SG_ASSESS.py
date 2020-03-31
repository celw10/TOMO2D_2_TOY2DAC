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

def read_fwdmdl_picks(forward, numobs, obs, del_obs_fwd, sr):
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
                
    #Re-arrange in proper format
    for o in np.arange(numobs):
        OBS_fw[o] = np.delete(OBS_fw[o],(0,2,3,5),axis=1)
        OBS_fw[o][:,1] = OBS_fw[o][:,1]/sr
        
    #print "THE THING THAT SHOULD BE DELETED: ", obs
    #Crop obs if requried
    if len(forward) != len (obs):
        obs = np.delete(obs,del_obs_fwd)
    #print "THE THING THAT SHOULD BE DELETED: ", obs
    
    #Shift to offset
    for o in np.arange(numobs):
        OBS_fw[o][:,0] = OBS_fw[o][:,0] - obs[o]
        
    #Reverse shot order
    OBS = []
    for i in range(len(OBS_fw)):
        it = len(OBS_fw) - 1 - i
        OBS.append(OBS_fw[it])

    return OBS

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

def data_cleanup(OBS, obsoffset, idx, reverse=True, del_obs_rev=None):
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
    
    #Delete any obsoffset values if using a cropped forward modeled first break pick set
    if del_obs_rev != None and len(OBS) != len(obsoffset):
        obsoffset = np.delete(obsoffset,del_obs_rev)
    
    #Delete any negative numbers (for some reason there are negative values in some FB pick files)
    tmp_del = [[] for _ in range(len(OBS))]
    for i in range(len(OBS)):
        for j in range(len(OBS[i][:,0])):
            if OBS[i][j,1] < 0: 
                tmp_del[i].append(j)
        OBS[i] = np.delete(OBS[i],(tmp_del[i]),axis=0)
    
    #Sort data by offset
    tmp = []
    for i in range(len(OBS)):
        tmp.append(OBS[i][np.argsort(OBS[i][:,0])])
        
    #Only want unique values
    tmp0 = []
    for i in range(len(OBS)):
        tmp0.append(tmp[i][np.unique(tmp[i][:,0], return_index=True, axis=0)[1]])
        
    #Round values to three decimal places
    for i in range(len(OBS)):
        for j in range(len(tmp0[i])):
            tmp0[i][j,0] = round(tmp0[i][j,0],3)
        
    #Interpolate first breaks over min/max traces
    tmp1 = []
    for i in range(len(tmp0)):
        #print "OBS, Min Index, Max Index: ", i, idx[0][i], idx[1][i]
        #print "OBS OFFSETS, Length of OBS Offsetes: ", obsoffset[i], len(obsoffset[i])
        #print "OBS offsets at min/max indicies: ", obsoffset[i][idx[0][i]], obsoffset[i][idx[1][i]]
        #print "Minimum and Maximum observed values: ", min(tmp0[i][:,0]), max(tmp0[i][:,0])
        #print
        #This is a fudge to fix errors in FB pick files
        #if obsoffset[i][idx[0][i]] > min(tmp0[i][:,0]):
        #    interp_min = min(tmp0[i][:,0])
        #else:
        #    interp_min = obsoffset[i][idx[0][i]]
        #if obsoffset[i][idx[1][i]] >
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

### Now we need to read in the data from SGY files

def indices_from_segy(OBS, offset, del_obs_rev=None):
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
    
    #Delete indixed OBS offsets if we are using forward modeled picks with a cropped OBS set
    if del_obs_rev != None and len(OBS) != len(offset):
        offset = np.delete(offset,del_obs_rev)
        
    #Store the maximum and minimum offsets for first break picks and from headers
    for i in np.arange(len(OBS)):
        minoffsets_ps[i] = round(min(OBS[i][:,0]),3)
        maxoffsets_ps[i] = round(max(OBS[i][:,0]),3)
        minoffsets_hdr[i] = round(min(offset[i]),3)
        maxoffsets_hdr[i] = round(max(offset[i]),3)
        
    #print minoffsets_ps, len(minoffsets_ps)
    #print maxoffsets_ps, len(maxoffsets_ps)
    #print
    
    #Find nearest
    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx

    #Save the maximum and minimum indicies of first break picks w.r.t. segy header offsets
    minimum_idx, maximum_idx, a, b = [], [], 0, 0
    for i in np.arange(len(offset)):
        #Minimum offsets
        j=0
        while j < len(offset[i]):
            if minoffsets_ps[i] == round(offset[i][j],3): #added round
                #print i,
                minimum_idx.append(j)
                #print minimum_idx, len(minimum_idx)-i
                j = len(offset[i])
            elif j == len(offset[i])-1:
                #If we haven't found anything get the offset value thats the closest
                print "Having trouble finding an matching minimum offset index for OBS: ", i+1
                minimum_idx.append(find_nearest(offset[i],minoffsets_ps[i]))
                j = len(offset[i])
            else:
                j += 1
        #Maximum offsets
        j=0
        while j < len(offset[i]):
            if maxoffsets_ps[i] == round(offset[i][j],3):
                #print i,
                maximum_idx.append(j)
                #print maximum_idx, len(maximum_idx)-i
                j = len(offset[i])
            elif j == len(offset[i])-1:
                #If we haven't found anything get the offset value thats the closest
                print "Having trouble finding an matching Maximum offset index for OBS: ", i+1
                maximum_idx.append(find_nearest(offset[i],maxoffsets_ps[i]))
                j = len(offset[i])
            else:
                j += 1
    
    # Length assertions
    assert len(minimum_idx) == len(offset)
    assert len(maximum_idx) == len(offset)
    
    #print
    #print minimum_idx
    #            
    #print minoffsets_ps
    #for i in range(len(offset)):
    #    print offset[i][minimum_idx[i]],
    #print
    #print maxoffsets_ps
    #for i in range(len(offset)):
    #    print offset[i][maximum_idx[i]],
    #    print
                
    return minimum_idx, maximum_idx

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

def shift_data(OBS, shotloc, bulk_shift, del_obs=None):
    """Shift entire dataset to have a origin of zero"""
    
    if bulk_shift != 0:
        print "WARNING: Bulk shift of ", bulk_shift, " applied to the acquisition file"
        print
        
    #Delete the appropiate shot locations if applicable, i.e. using a forward modeled dataset that has already been cropped
    if del_obs != None and len(OBS) != len(shotloc):
        shotloc = np.delete(shotloc,del_obs)
    
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
    """

    ########################
    ### RESTRUCTURE DATA ###
    ########################

    #Restructure data_modeling file, only reading frequencies we desire
    it=0
    print "RESTRUCTURING DATA_MODELING FILE..."
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
    else:
        RealMatrix= np.zeros((nfreq, nchannels))
        ComplexMatrix = np.zeros((nfreq, nchannels))
        while it < relfreq:
            for c in range(0,2*num_traces[source],2):
                RealMatrix[it,c/2] = data[c+(sum(num_traces[0:source]*2))+(it*sum(num_traces)*2)]
                ComplexMatrix[it,c/2] = data[c+(sum(num_traces[0:source]*2))+(it*sum(num_traces)*2)+1]
            it+=1
    #Form the complex matrix
    if sources == all:
        FD_Comp_SG = np.zeros((nsources, int(nfreq), nchannels), dtype=complex)
        for s in range(nsources):
            for it in range(int(relfreq)):
                for c in range(num_traces[s]):
                    FD_Comp_SG[s,it,c] = complex(RealMatrix[s,it,c],ComplexMatrix[s,it,c])
    else:
        FD_Comp_SG = np.zeros((int(nfreq), nchannels), dtype=complex)
        for it in range(int(relfreq)):
            for c in range(num_traces[source]):
                FD_Comp_SG[it,c] = complex(RealMatrix[it,c],ComplexMatrix[it,c])
    print "RESTURCTURING COMPLETE"
    print
    
    return FD_Comp_SG

#MAKE FOURIER GATHERS

def fourier_gathers(OBS, gathers, obsoffset, obs, amp_corr, 
                    del_idx, indicesdel_min, indicesdel, min_idx, 
                    nsamples, maxfold, buff, div, sr, off_min,
                    LFC=0, plot=2):
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
                            
                    #Apply laplace damping from t=0 - Shin and Cha 2009, Brossier et al., 2009 formulation
                    elif (j + buff) > int(OBS[num][i,1]) and LFC<0:
                        #Taper the buffer zone with a hanning (no hard edges)
                        #This is different for t=0 damping, first save the point, then taper
                        if j < int(OBS[num][i,1]):
                            gathers_fourier[num][j,i] = \
                            gathers[num].T[j,i+int(index[num])]*e**(-(j*sr)*abs(LFC))*amp_corr[num][i]
                            gathers_fourier[num][j,i] = gathers_fourier[num][j,i]*(taper_half[buff-(int(OBS[num][i,1])-j)])
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
                        
    print "Traces removed, min offsets: : ", total
        
    return gathers_fourier, OBS, max_pick_num

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
    
    return TD_SG

"""
Chris Williams, October 2019
-> plot_single_sg divided into two bits due to the nature of this problem
-> This bit, convert everything, observed & modeled datasets, and FB picks to reduced time
FUNCTION DEFINITIONS:

"""

def convert_red_time(data_observed, data_modeled, shot, min_max, acqui, sr, 
                     maxplot_t, reduction_vel, red_maxplot_t, del_neg=True,
                     max_off=100000, fb_pick1=False, fb_pick2=False):
    
    """CONVERT TO REDUCED TIME"""
    ### OBSERVED DATASET ###
    plot_data_obs = data_observed[0:int(maxplot_t/sr),0:len(acqui[shot][1::,1])] 
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
            red_plot_data[tnodes[int(i*maxplot_t/sr+j)]-min_tnode-1,i] = plot_data_obs[j,i]
    #Tidy up time-reduced shot gather
    max_minredtime = maxplot_t - (max_off/reduction_vel)
    if max_minredtime < 0:
        print "Please use a larger reduction velocity, negatives encountered."
    else:
        data_yextent = abs(min_redtime)
        data_textent = 0
    data_obs = red_plot_data[int(data_yextent/sr):int((red_maxplot_t+data_yextent)/sr), 0:len(acqui[shot][:,0])-1]
    print
    
    if type(data_modeled) != bool:
        ### PREDICTED DATASET ###
        plot_data_mdl = data_modeled[int(len(data_modeled[:,0])-(maxplot_t/sr)):len(data_modeled[:,0]), 
                               0:len(acqui[shot][1::,1])] 
        plot_data_mdl = np.flip(plot_data_mdl, axis=0)
        #Plotting specific definitions
        y_axis_label = "Time - Offset/"+str(reduction_vel/1000.0)+" [s]"
        #Extract reduced time nodes
        tnodes=[]
        print "Converting plot to reduced time using a reduction velocitiy of: ", reduction_vel
        for i in range(len(acqui[shot][1::,1])):
            x = acqui[shot][i+1,1] - acqui[shot][0,1]
            for j in range(int(maxplot_t/sr)):
                tnodes.append(int((j*sr-(abs(x)/reduction_vel))/sr))
        #Extract min/max fvalues
        max_tnode, min_tnode = max(tnodes), min(tnodes)
        max_redtime, min_redtime = max_tnode*sr, min_tnode*sr
        #Fill reduced-time shot gather
        red_plot_data = np.zeros((int(max_tnode-min_tnode),len(acqui[shot][1::,1])))
        print np.shape(red_plot_data)
        for i in range(len(acqui[shot][1::,1])):
            for j in range(int(maxplot_t/sr)):
                red_plot_data[tnodes[int(i*maxplot_t/sr+j)]-min_tnode-1,i] = plot_data_mdl[j,i]
        #Tidy up time-reduced shot gather
        max_minredtime = maxplot_t - (max_off/reduction_vel)
        if max_minredtime < 0:
            print "Please use a larger reduction velocity, negatives encountered."
        else:
            data_yextent = abs(min_redtime)
            data_textent = 0
        data_mdl = red_plot_data[int(data_yextent/sr):int((red_maxplot_t+data_yextent)/sr), 0:len(acqui[shot][:,0])-1]
        print
    else:
        data_mdl = data_modeled
    
    """First Break Picks to Reduced Time"""
    
    if type(fb_pick1) != bool:
        #Extract reduced time nodes
        tnodes_fb1 = np.zeros((len(fb_pick1[shot][:,0]),2))
        #I need to fill out the x somehow
        print "Converting first break picks to reduced time using a reduction velocitiy of: ", reduction_vel
        for i in range(len(fb_pick1[shot][:,0])):
            x = (fb_pick1[shot][i,0]-min_max[0])*1000.0 - acqui[shot][0,1]
            tnodes_fb1[i,1] = ((fb_pick1[shot][i,1])*sr-(abs(x)/reduction_vel))
    
    if type(fb_pick2) != bool:
        #Extract reduced time nodes
        tnodes_fb2 = np.zeros((len(fb_pick2[shot][:,0]),2))
        #I need to fill out the x somehow
        print "Converting first break picks to reduced time using a reduction velocitiy of: ", reduction_vel
        for i in range(len(fb_pick2[shot][:,0])):
            x = (fb_pick2[shot][i,0]-min_max[0])*1000.0 - acqui[shot][0,1]
            tnodes_fb2[i,1] = ((fb_pick2[shot][i,1])*sr-(abs(x)/reduction_vel))

    #We need to figure out the relative shift between first break picks
    #The first group of first break picks are the reference, of which the shot gather whas generated
    if type(fb_pick1) != bool and type(fb_pick2) != bool:
        i, j, flag = 0, 0, "missing"
        while i < len(fb_pick1[shot][:,0]):
            while j < len(fb_pick2[shot][:,0]):
                if round(fb_pick2[shot][j,0],3) == round(fb_pick1[shot][i,0],3):
                    out_shift, flag = i-j, "found"
                    i, j = len(fb_pick1[shot][:,0]), len(fb_pick2[shot][:,0])
                else:
                    j+=1
            if flag == "found":
                print "Shift required relative to reference FB Picks: ", out_shift
            else:
                i+=1
                j=0
                out_shift=0
        print "OUT SHIFT: ", out_shift
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
    if type(fb_pick1) != bool:
        for i in range(len(fb_pick1[shot][:,0])):
            tnodes_fb1[i,0] = i
    
    print
    print "Shape Observed Data: ", np.shape(data_obs)
    if type(data_mdl) != bool:
        print "Shape Modeled Data: ", np.shape(data_mdl)
    if type(fb_pick1) != bool:
        print "Shape FB Picks, 1: ", np.shape(tnodes_fb1)
    if type(fb_pick2) != bool:
        print "Shape FB Picks, 2: ", np.shape(tnodes_fb2)
    
    if type(fb_pick1) != bool and type(fb_pick2) == bool:
        return data_obs, data_mdl, tnodes_fb1
    elif type(fb_pick1) != bool and type(fb_pick2) != bool:
        return data_obs, data_mdl, tnodes_fb1, tnodes_fb2
    else:
        return data_obs, data_mdl
    
"""
Chris Williams, September 2019
-> plot_single_sg divided into two bits due to the nature of this problem
-> This bit, plot everything up, going to overlay a transparent modeled dataset over an observed dataset.
FUNCTION DEFINITIONS:
"""


def plot_SG_overlay(data_observed, data_modeled, sr, cbar_pos, title_pos, maxplot_t, reduction_vel, red_maxplot_t, 
                    vrange=None, vrange_diff=None,
                    crop=None, alpha=0.5, tick_inc=10,  aspect=20, fontsize=1, size=[32,32], 
                    max_off=100000, savepath=None, savename=None,
                    fb_pick1=bool, fb_pick2=bool, fb_labels=None, residuals=False):

    """DEFINE FIGURE SIZE"""
    size = size
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = size[0]
    fig_size[1] = size[1]
    plt.rcParams["figure.figsize"] = fig_size
    
    """DEFINE COLORBAR RANGE"""
    if vrange==None:
        elev_min, elev_max = min(data_observed.flatten())/100.0, -min(data_observed.flatten())/100.0
    else:
        elev_min, elev_max = vrange[0], vrange[1]
    print "Colorbar range (min/max): (", elev_min, "/", elev_max, ")"
    
    """CROP: WORK ON THIS IN THE FUTURE"""
    
    """COMPUTE RESIDUALS"""
    if residuals == True:
        #Modeled data subtracted from observed data
        diff = data_observed - data_modeled
        #Colorbar range for difference plot, define limits based on min/max difference or user specified
        if vrange_diff == None:
            elev_min_diff, elev_max_diff = min(diff.flatten())/100.0, -min(diff.flatten())/100.0
            mid_val=0
        else:
            elev_min_diff, elev_max_diff = vrange_diff[0], vrange_diff[1]
            mid_val=0
    
    """DEFINE PLOT EXTENT"""
    extent=[0,len(acqui[shot][1::,1]),red_maxplot_t,0]
    
    """MAKE PLOTS"""
    fig = plt.figure() 
    if residuals==True:
        ax = fig.add_subplot(211)
        dif = fig.add_subplot(212)
    else:
        ax = fig.add_subplot(111)
    #Shot Gathers
    im_obs = ax.imshow(data_observed, aspect=aspect, vmin=elev_min, vmax=elev_max, cmap="gray_r", extent=extent)
    if type(data_modeled) != bool:
        im_mdl = ax.imshow(data_modeled, aspect=aspect, vmin=vrange_diff[0], vmax=vrange_diff[1], cmap="seismic_r", #TMP COLOR BAR CHANGE vmin=elev_min, vmax=elev_max,
                           extent=extent, alpha=alpha)
    if residuals == True:
        #Residual Plot
        im_dif = dif.imshow(diff, aspect=aspect, vmin=elev_min_diff, vmax=elev_max_diff, cmap='seismic', 
                             extent=extent, norm=MidpointNormalize(midpoint=mid_val))
    #FB Picks
    if type(fb_pick1)!=bool:
        ax.plot(fb_pick1[:,0],fb_pick1[:,1], "b-", linewidth=2.0)#, label=fb_labels[0])
    if type(fb_pick2)!=bool:
        ax.plot(fb_pick2[:,0],fb_pick2[:,1], "r--", linewidth=3.0)#, label=fb_labels[1])
    
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
    ax.set_xticks(x_tick_index)
    ax.set_xticklabels(x_tick_values)
    
    #Titles and plot labels (First Plot)
    ax.set_title('OBS '+str(shot+1),fontweight="bold", #
                 size=20.0*fontsize, x=0.14+title_pos[0], y=0.95+title_pos[1]) 
    y_axis_label = "Time - Offset/"+str(reduction_vel/1000.0)+" [s]"
    ax.set_ylabel(y_axis_label, fontsize = 18.0*fontsize)#plt.ylabel(y_axis_label, fontsize = 18.0*fontsize) # Y label
    ax.set_xlabel('Offset [km]', fontsize = 18.0*fontsize) # X label
    ax.xaxis.set_label_position('top') 
    ax.xaxis.tick_top()
    ax.tick_params(axis='both', which='major', labelsize=16.0*fontsize)
    ax.tick_params(axis='both', which='minor', labelsize=16.0*fontsize)
    ax.grid()
    
    if residuals == True:
        #Titles and plot labels (Second Plot)
        dif.set_title('Residual Wavefield, OBS '+str(shot+1),fontweight="bold", 
                     size=20.0*fontsize, x=0.11+title_pos[0], y=0.95+title_pos[1]) 
        y_axis_label = "Time - Offset/"+str(reduction_vel/1000.0)+" [s]"
        dif.set_ylabel(y_axis_label, fontsize = 18.0*fontsize)#plt.ylabel(y_axis_label, fontsize = 18.0*fontsize) # Y label
        dif.set_xlabel('Offset [km]', fontsize = 18.0*fontsize) # X label
        dif.xaxis.set_label_position('top') 
        dif.xaxis.tick_top()
        dif.tick_params(axis='both', which='major', labelsize=16.0*fontsize)
        dif.tick_params(axis='both', which='minor', labelsize=16.0*fontsize)
        dif.grid()
    
    #Colorbar
    ax_cbar1 = fig.add_axes([cbar_pos[0], cbar_pos[1], 0.22, 0.03]) #0.68, 0.263 OBS 12
    cbar1 = fig.colorbar(im_obs, cax=ax_cbar1, ax=ax, orientation='horizontal', format='%.0e')
    cbar1.ax.tick_params(labelsize=10.0*fontsize)
    cbar1.set_label("Amplitude Observed", fontsize=14.0*fontsize)
    
    if type(data_modeled) != bool:
        ax_cbar2 = fig.add_axes([cbar_pos[0], cbar_pos[1]-0.05, 0.22, 0.03]) #0.68, 0.263 OBS 12
        cbar2 = fig.colorbar(im_mdl, cax=ax_cbar2, ax=ax, orientation='horizontal', format='%.0e')
        cbar2.ax.tick_params(labelsize=10.0*fontsize)
        cbar2.set_label("Amplitude Modeled", fontsize=14.0*fontsize)

    if residuals==True:
        ax_cbar3 = fig.add_axes([cbar_pos[0], 0.177, 0.22, 0.03]) #Second value to vertically position
        cbar3 = fig.colorbar(im_dif, cax=ax_cbar3, ax=dif, orientation='horizontal', format='%.0e')
        cbar3.ax.tick_params(labelsize=10.0*fontsize)
        cbar3.set_label("Amplitude Difference", fontsize=14.0*fontsize)
    
    #Plot legend
    if fb_labels!=None:
        ax.legend(loc="upper left", fontsize=14*fontsize)
        
    if savepath != None and savename != None:
        plt.gcf().savefig(savepath+savename+".png", bbox_inches = 'tight', pad_inches = 0.1)
    
    plt.show()
    
import matplotlib.colors as colors
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
    
 #############################################################################################

#The data generation step is a memory intensive process
#This code seems to free up memory nicely, the source is provided below
#https://stackoverflow.com/questions/32167386/force-garbage-collection-in-python-to-free-memory

def mem():
    print('Memory usage         : % 2.2f MB' % 
          round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0,1))

def memoryhog():
    #Bring in files from disc
    gathers = np.load("gathers/gathers.npy", allow_pickle=True)

    print
    print "Generating Shot Gathers..."
    print 
    
    ######################
    ### OUTPUT GATHERS ###
    ######################

    """REAL DATA"""
    #Prep gathers for fourier transform
    RD_SG, OBS_fourier, max_pick_num = fourier_gathers(OBS, gathers, obsoffsets, obs, amp_corr,
                                                       del_idx, indicesdel_crop, indicesdel, min_idx,
                                                       nsamples, maxfold, buff, div, sr, offset_min,
                                                       LFC=LFC)

    #if num_fb_picks == 2: 
        #Delete proper indicies for OBS 1, not required if you're already using an OBS with cropped picks!
    #    OBS_1 = np.delete(OBS_1, del_idx_1)

    #Save optional FB picks
    if num_fb_picks >= 1:
        fb_pick1 = OBS_fourier
    else:
        fb_pick1 = False

    if num_fb_picks == 2:
        fb_pick2 = OBS_1
    else:
        fb_pick2 = False

    #Only plotting one SG to save memory
    Observed_Dataset = RD_SG[shot]

    #Save some memory 
    del gathers, RD_SG
    gc.collect()

    if type_data == 0:
        """TOY2DAC DATA DATA"""

        #Filter frequency domain data
        FD_Comp_SG_taper = taper_FD_data(FD_Comp_SG, taper_design, nsources, nfreq, relfreq, nchannels, sources=all)

        #Inverse fourier transform to time domain
        FM_SG = inv_fft(FD_Comp_SG_taper, nsources, nfreq, nchannels, sources=all)

        #Only plotting one SG to save memory
        Modeled_Dataset = FM_SG[shot]

        #Save some memory
        del FD_Comp_SG_taper, FM_SG
        gc.collect()

    ##############################################
    ### OUTPUT REDUCED TIME MATRICIES/FB PICKS ###
    ##############################################

    if type_data == 1:
        Modeled_Dataset = False

    if num_fb_picks == 0:
        dobs, dmdl = convert_red_time(Observed_Dataset, Modeled_Dataset, shot, min_max, acqui, sr, 
                                                int(tmax/div), 6500, 10, del_neg=True, 
                                                max_off=100000, fb_pick1=fb_pick1, fb_pick2=fb_pick2)
        np.save("data/observed", dobs)
        np.save("data/modeled", dmdl)
    if num_fb_picks == 1:
        dobs, dmdl, fb1 = convert_red_time(Observed_Dataset, Modeled_Dataset, shot, min_max, acqui, sr, 
                                                int(tmax/div), 6500, 10, del_neg=True, 
                                                max_off=100000, fb_pick1=fb_pick1, fb_pick2=fb_pick2)
        np.save("data/observed", dobs)
        np.save("data/modeled", dmdl)
        np.save("data/picks1", fb1)
    if num_fb_picks == 2:
        dobs, dmdl, fb1, fb2 = convert_red_time(Observed_Dataset, Modeled_Dataset, shot, min_max, acqui, sr, 
                                                int(tmax/div), 6500, 10, del_neg=True, 
                                                max_off=100000, fb_pick1=fb_pick1, fb_pick2=fb_pick2)
        np.save("data/observed", dobs)
        np.save("data/modeled", dmdl)
        np.save("data/picks1", fb1)
        np.save("data/picks2", fb2)
        
###TEST WHAT I'VE DONE SO FAR###

### 2.0 Update: December 5th 2019, I've tried to make this a little more memory efficent

import numpy as np
import segyio #This doesn't come standard with python
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.colors as mcolors
import struct

from scipy import interpolate
from sympy import N
from math import e
import multiprocessing as mp
import struct
import resource
import gc

#############
### FLAGS ###
#############

#Number of first break picks to be overlain on top of the gathers, 0=none, 1=final, 2=final and optional FB#1
num_fb_picks = 2

#Type of optional first break pick, picked=0 or forward modeled=1
type_fb_pick = 1

#Include toy2dac forward modeled data=0, or time domain data only=1
type_data = 1

#############
### FILES ###
#############

###############
"""REAL DATA"""
###############

#Full path to tx.in file
f = "/home/chris/Documents/1_data/EasternMed/Chris_Williams/first_breaks/pick9a.final.zerophase"
picks = np.genfromtxt(f) 

#Full path to FB pick file 1 - Model first break picks
f = "/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/110_pick9.final.zerophase_108run_Oct29/3_modeled_traveltimes/09_fwdtts"
picks1 = np.genfromtxt(f,skip_header=1)

#Labels
fb_labels=["Picked", "Modeled"]

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
#segydir = '/home/chris/0_GlobeClaritas/0_Projects/Ref_EMed_Line1/COMMON/DATA/5_FinalFlow_Sept30th/3_FWIData_RMSNorm/'
segydir = '/home/chris/0_GlobeClaritas/0_Projects/Ref_EMed_Line1/COMMON/DATA/5_FinalFlow_Sept30th/3_FWIData_RMSNorm/'
name_start = ''
#name_end = '_VERT_Final_FWI.sgy'
name_end = '_VERT_Final_FWI.sgy'
fileloc = [segydir, name_start, name_end]

#Amplitude correction file generated externally (usually saved with the data_modeling file for the starting model)
f = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/6_ForwardModeling/009_Oct29_StartMdl110/AmpCorr_QP050_sdNov12.txt"
amp_corr = np.genfromtxt(f)

##################
"""TOY2DAC DATA"""
##################

#Open toy2dac data_modeling file
dirr = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/008_2020RealDataTests/038_MiscTests_FGProg1_lBFGS/8_fwdmdl/" #Path to root directory
model = "dm_fg15Hz_mdl043"
filename = dirr + model #Name of binary file
with open(filename, 'rb') as f:
    data = np.fromfile(f, dtype=np.float32)

#TOY2DAC Geometry file
f = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/000_DATA/2_RealData/17_110StartMdl_THEEND/01_acqui_tfb/15_30_45_60_100/acqui_100" 

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
offset_max = 100 #Has to be 100 reading in mdl data
offset_min = 0
bulk_shift = 0 # Not really implimented
min_max = [30, 220]

#Indicies of OBS that have already been cropped out, blank if none: PYTHON INDEXING, OBS number - 1
del_obs = [0, 1, 2, 11, 13]
#Indicate what OBS are being deleted if using forward modeled picks (Python indexing!)
del_obs_fwd = [0, 1, 2, 11, 13]
del_obs_rev = [7, 9, 18, 19, 20]

#Define variables for fourier transform
buff = 100 #retain this many time samples above the first break
div = 1.5 #divide the total number of time samples by this number (saves memory)
LFC = 0.0 #the Laplace-Fourier Constant, set to zero if it isn't to be used, negative constant applys damping from t=0

#Plotting Definitions
shot = 12 #Python Indexing, index+1

#Saving reduced-time plot
savepath = None
savename = None

##################
"""TOY2DAC DATA"""
##################

#Data-related
nchannels = 660 #Number of channels, obsolete if setup_flag != 0 
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
if num_fb_picks == 2 and type_fb_pick == 0:
    #Optional FB picks/forward modeled data #1
    OBS_1, pick_numbers = read_fbpicks(picks1, obs, phase, unc, tt_window, sr)
elif num_fb_picks == 2 and type_fb_pick == 1:
    OBS_1 = read_fwdmdl_picks(picks1, 16, obs, del_obs_fwd, sr)# Enter appropiate number of OBS

if type_data == 0: 
    """TOY2DAC DATA"""
    #Setup geometry from toy2dac acqui file
    num_traces, nchannels = setup_geom(f, setup_flag=1, sources=all)

#################
### SEGY DATA ###
#################

"""REAL DATA"""
#Extract obsoffsets from the header and the gathers from sgy file
gathers, obsoffsets = extract_segy(fileloc, len(obs), maxfold, nsamples)

##### TMP
#vrange = [-5e-4,5e-4]
#import matplotlib.pyplot as plt
#plt.imshow(gathers[0], vmin=vrange[0], vmax=vrange[1], cmap="gray")
#plt.show()
    
#plt.imsh
#########

#Save to disc
np.save("gathers/gathers", gathers)

#Save some memory extract_segy
del gathers
gc.collect()

#Indexes for first break picks in segy data
min_idx, max_idx = indices_from_segy(OBS, obsoffsets)
if num_fb_picks == 2:
    #Optional FB picks/forward modeled data #1
    min_idx1, max_idx1 = indices_from_segy(OBS_1, obsoffsets, del_obs_rev=del_obs_rev)
    
################
### CLEAN UP ###
################

"""REAL DATA"""
#Ensure uniqueness and interpolate gaps in dataset.
OBS, pick_num_data = data_cleanup(OBS, obsoffsets, [min_idx, max_idx])
if num_fb_picks == 2:
    #Optional FB picks/forward modeled data #1
    OBS_1 = data_cleanup(OBS_1, obsoffsets, [min_idx1, max_idx1], del_obs_rev=del_obs_rev)

if type_data == 0:    
    """TOY2DAC DATA"""
    #Restructure toy2dac data_modeling file into matrix format
    FD_Comp_SG = restructure_toy2dac_data(data, nsources, nfreq, relfreq, nchannels, num_traces, sources=all)
    #Save some memory 
    del data
    gc.collect()

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
#~~~ Optional FB picks/forward modeled data #1
if num_fb_picks == 2:
    OBS_1 = OBS_1[0]
    #Set the maximum and minimum values for model cropping
    OBS_1, indicesdel_1 = max_offset(OBS_1, offset_max)
    #Shift the data with the OBS model locations and a bulk shift
    OBS_1 = shift_data(OBS_1, obs, bulk_shift, del_obs=del_obs_fwd)
    #Set the maximum and minimum values for model cropping
    OBS_1, indicesdel_crop_1 = crop_model(OBS_1, min_max)

#########################
### OUTPUT ACQUI FILE ###
#########################

"""REAL DATA"""
#Output final version of acqui file
acqui, del_idx = gen_toy2dac_acqi(OBS, obs, obs_z, bulk_shift, min_max, del_obs)
if num_fb_picks == 2:
    #Optional FB picks/forward modeled data #1
    acqui_1, del_idx_1 = gen_toy2dac_acqi(OBS_1, obs, obs_z, bulk_shift, min_max, del_obs)

##################################
### OUTPUT GATHERS & RED. TIME ###
##################################

"""
Memory intensive step, defined as a function above "memoryhog"
"""

#Run Memory Intensive Code
mem()
proc = mp.Process(target=memoryhog)
proc.start()
proc.join()
mem()

#Everything for plotting is saved to disc under "data"

##################################
### PLOT OVERLAIN SHOT GATHERS ###
##################################

"""CLOSE: Two things left to fix. Colorbar positioning, and model cropping!"""

#If running in Jupyter notebook, you might need to re-open the notebook as such:
#jupyter notebook --NotebookApp.iopub_data_rate_limit=10000000000

dobs = np.load("data/observed.npy", allow_pickle=True)
#dmdl = np.load("data/modeled.npy", allow_pickle=True)
fb1 = np.load("data/picks1.npy", allow_pickle=True)
fb2 = np.load("data/picks2.npy", allow_pickle=True)

#Definitions
vrange = [-5e-4,5e-4] #Optional, dependent on amp corr file
vrange_diff = [-5e-4, 5e-4]
cbar_pos = [0.68, 0.389]#10=0.68, 0.3666
title_pos = [0.78,-0.08]#15

plot_SG_overlay(dobs, False, sr, cbar_pos, title_pos, int(tmax/div), 6500, 10, vrange, vrange_diff=vrange_diff,
                    alpha=0.25, tick_inc=10,  aspect=30, fontsize=2.0, size=[32,32], 
                    max_off=100000, savepath=None, savename=None,
                    fb_pick1=fb1, fb_pick2=fb2, fb_labels=["First Break Picks"], residuals=False)