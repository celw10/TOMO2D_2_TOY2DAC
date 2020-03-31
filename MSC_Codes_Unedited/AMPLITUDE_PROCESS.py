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

def fourier_gathers(OBS, gathers, obsoffset, obs, 
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
    #Record the noise profile above FB's to scale down noisy traces
    noise_profile = np.zeros((len(OBS), max_pick_num, 2))
    total=0
    for num in range(len(OBS)):
        for i in range(len(OBS[num][:,0])):
            #Zero out the gather if within the minimum offset
            if abs(OBS[num][i,0]-shotgeom[num]) < off_min:
                gathers_fourier[num][:,i] = 0
                total+=1
            else:
                for j in range(0,int((nsamples-1)/div),1): 
                    #Apply the buffer above first break pick, optionally apply time damping to the data
                    if (j + buff) > int(OBS[num][i,1]) and LFC!=0:
                        gathers_fourier[num][j,i] = \
                        gathers[num].T[j,i+int(index[num])]*e**(-(j*sr-OBS[num][i,1]*sr)/LFC)
                    elif (j + buff) > int(OBS[num][i,1]) and LFC==0:
                        gathers_fourier[num][j,i] = \
                        gathers[num].T[j,i+index[num]]
                    #Save the average amplitude values above the first break plus buffer, highlight noisy traces
                    #elif (j + 100) < int(OBS[num][i,1]):
                    else:
                        #First column for the sum of the amplitudes
                        noise_profile[num][i,0] += (gathers[num].T[j,i+int(index[num])])**2
                        #Second column for the number of time samples
                        noise_profile[num][i,1] += 1
                        
    print "Traces removed, min offsets: : ", total
    
    #Gathers that go into TOY2DAC have to be flipped
    plt_data = gathers_fourier[plot]
    gathers_fourier = np.flip(gathers_fourier,1) 

    #Plot Time-domain Gathers with FB picks
    if plot!=False:
        plt.rcParams["figure.figsize"] = (16,12)
        plt.imshow(plt_data, aspect=max_pick_num/(nsamples/div), vmin=-max(plt_data.flatten()/100),
                   vmax=max(plt_data.flatten()/100), cmap="seismic")
        cbar=plt.colorbar()
        plt.plot(range(len(OBS[plot][:,0])),OBS[plot][:,1])
        plt.title("Time-domain data, shot: "+str(plot))
        plt.xlabel('Channel')
        plt.ylabel('Time sample')
        plt.show()
        
    return gathers_fourier, OBS, max_pick_num, noise_profile

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

def compute_NRMS(noise_profile, OBS, minoffset, max_pick_num):
    """
    Compute the RMS amplitude above the first break, or the RMS noise for each trace.
    We are looking for anomalously noisy traces, and we wish to scale these traces down.
    We only want to scale down traces outside of the minimum offset.
    -> There are reverberations from the application of the filter above the first break at near offset
    -> Although the underlying traces themselves are not noisy
    """
    
    #Shift the geometry from model to offset domain
    OBS_offset = [[] for _ in range(len(OBS))]
    del_obs_loc = np.delete(obs,del_idx,0)
    for i in range(len(del_obs_loc)):
        OBS_offset[i] = OBS[i]-del_obs_loc[i]
        
    #RMS matricies
    RMS_preFB_amp = np.zeros((len(OBS), max_pick_num))
    RMS_preFB_amp_interp = np.zeros((len(OBS), max_pick_num))
    
    #Compute the noise RMS 
    for s in range(len(OBS)):
        #Remove pesky nans
        vals = np.isnan(RMS_preFB_amp[s])
        RMS_preFB_amp[s][vals] = 0
        for i in range(len(OBS[s])):
            #We need to keep a value at maximum and minimum offsets for interpolation
            minoffset_tmp = min(OBS_offset[s][:,0])
            maxoffset_tmp = max(OBS_offset[s][:,0])
            if OBS_offset[s][i,0] == minoffset_tmp:
                RMS_preFB_amp[s,i] = sqrt(noise_profile[s][i,0]/noise_profile[s][i,1])
            elif OBS_offset[s][i,0] == maxoffset_tmp:
                RMS_preFB_amp[s,i] = sqrt(noise_profile[s][i,0]/noise_profile[s][i,1])
            #We only want to scale down amplitudes of far offsets
            elif abs(OBS_offset[s][i,0]) > minoffset:
                RMS_preFB_amp[s,i] = sqrt(noise_profile[s][i,0]/noise_profile[s][i,1])
            #Set short offsets RMS amplitude values to zero for now, interpolate zeros later
            else:
                RMS_preFB_amp[s,i] = 0
        #Interpolate zeros        
        y = np.array(RMS_preFB_amp[s,0:len(OBS[s])])
        xnew = np.arange(len(OBS[s]))
        idx0 = np.where(y==float(0))
        xold = np.delete(xnew,idx0)
        yold = np.delete(y, idx0)
        f = interpolate.interp1d(xold,yold)
        RMS_preFB_amp_interp[s,0:len(OBS[s])] = f(xnew)
        
    return RMS_preFB_amp_interp, OBS_offset

def compute_SD(RMS_preFB_amp_interp, OBS_offset, max_pick_num):
    """
    Compute the standard deviations from the interpolated NRMS profile calculated previously.
    Note: The SD values may be quite small dependent on how you've scaled the data.
    -> I've multiplied the entire dataset for each OBS so the RMS amplitude is 1.0 for example
    -> I now have very small SD's
    """
    
    #Standard deviation matrix
    Strd_Dev = np.zeros((len(OBS_offset), max_pick_num))
    
    #Compute the standard deviation for each point
    for s in range(len(OBS_offset)):
        for i in range(len(OBS_offset[s])):
            a = (RMS_preFB_amp_interp[s,i]-sum(RMS_preFB_amp_interp[s])/len(OBS_offset[s]))**2
            b = (1/float(len(OBS_offset[s])-1))
            #Exaggerate the standard deviation
            Strd_Dev[s,i] = sqrt(a*b)
            
    return Strd_Dev

"""
Chris Williams, July 19th 2019.
We are going to fit our noise profile with a linear function.
The linear function that best fits the noise profile is computed for by a linear inversion.
We calculate and incorporate standard deviations into this computation.
The linear function is to be considered the background noise profile. 
We wish to lower noisy traces to this background noise level. 
"""

import sympy as sy

def GN_2term_fit(f, Bi, data_obs, S, std_dev, Plot=True, Verbose=True,
                 x = sy.Symbol('x'), B0 = sy.Symbol('B0'), B1 = sy.Symbol('B1')):

    #Initial Guess - put into a matrix, setup a matrix for weighted result if errorbar is flagged
    b0, b1 = Bi[0], Bi[1]

    B = np.matrix([[b0],[b1]])
    
    #Serup the Jacobian and residual matricies, first is for the control, second is weighted
    rows, cols = len(S), 2

    Jfw, rw, r = np.zeros((rows,cols)), np.zeros((rows,1)), np.zeros((rows,1)) 
    
    #Define the function, and setup the partial derivatives to be filled into the Jacobian matrix
    def model(f, b0, b1, xi, x = sy.Symbol('x'), 
                     B0 = sy.Symbol('B0'), B1 = sy.Symbol('B1')):
        return f.subs(x,xi).subs(B0, b0).subs(B1, b1)

    def partialDerB0(f, b0, b1, xi, x = sy.Symbol('x'), 
                     B0 = sy.Symbol('B0'), B1 = sy.Symbol('B1')):
        return - f.diff(B0,1).subs(x,xi).subs(B0, b0).subs(B1, b1)

    def partialDerB1(f, b0, b1, xi, x = sy.Symbol('x'), 
                     B0 = sy.Symbol('B0'), B1 = sy.Symbol('B1')):
        return - f.diff(B1,1).subs(x,xi).subs(B0, b0).subs(B1, b1)

    def residual(f, xi, observed, b0, b1, x = sy.Symbol('x'), 
                     B0 = sy.Symbol('B0'), B1 = sy.Symbol('B1')):
        return (observed - f.subs(x, xi).subs(B0, b0).subs(B1, b1))
    
    def LSQMisfit(f, xi, observed, b0, b1, x = sy.Symbol('x'), 
                     B0 = sy.Symbol('B0'), B1 = sy.Symbol('B1')):
        return sqrt((observed - f.subs(x, xi).subs(B0, b0).subs(B1, b1))**2)
    
    #Fill the Jacobian and residual matricies, solve for the updated fit
    TotalLSQMisfit = 0
    for j in range(0,rows,1):
        #Fill the residual matrix and the LSQ Misfit matrix
        rw[j,0] = (1.0/std_dev[j])*residual(f, S[j], data_obs[j], B[0], B[1])
        Jfw[j,0] = (1.0/std_dev[j])*partialDerB0(f, b0, b1, S[j])
        Jfw[j,1] = (1.0/std_dev[j])*partialDerB1(f, b0, b1, S[j])

        #Compute Misfit
        TotalLSQMisfit += LSQMisfit(f, S[j], data_obs[j], B[0], B[1])
    if Verbose==True:
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("Starting B values: " + '\n' + str(B) + '\n' + ", Starting LSQ misfit: " 
              + str(sy.N(TotalLSQMisfit,5)))
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

    #Solve for the updated points
    Jftw =  Jfw.T
    B = B - np.dot(np.dot(np.linalg.inv(np.dot(Jftw,Jfw)),Jftw),rw)

    #Compute updated misfit and un-weighted residual
    TotalLSQMisfit = 0
    for j in range(0,rows,1):
        TotalLSQMisfit += LSQMisfit(f, S[j], data_obs[j], B[0], B[1])
        r[j,0] = residual(f, S[j], data_obs[j], B[0], B[1])
    if Verbose==True:
        print("Updated B values: " + '\n' + str(B) + '\n' + ", Updated LSQ misfit: " 
              + str(sy.N(TotalLSQMisfit,5)))
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    
    return r, B

#Define function to get the NRMS trace scaling values

def NRMS_trace_scale(RMS_preFB_amp_interp, Strd_Dev, SD_Target, OBS_offset, maxoffset, scale,
                     plotval=None, fontscale=None):
    """
    Compute the value by which defined high amplitude noisy traces are to be multiplied by to shift them to a 
    background noise level.
    The SD_Target is a defined standard deviation level, above which the traces is considered anomalously noisy.
    A the residual of between the line of best fit and the points RMS is multiplied by a scaler [scale] prior to 
    shifting.
    Ensure you're QC'ing what the standard deviations are seeing by looking at the shot records.
    plotval is a QC tool to plot the results. Simply set it to the OBS you wish to plot. 
    """

    #Inputs for NRMS Scaling
    NRMS_scaling = [[] for _ in range(len(OBS_offset))]
    
    #Definitions for 2 term linear fit
    #Variables
    x, B0, B1 = sy.Symbol('x'), sy.Symbol('B0'), sy.Symbol('B1')
    #linear function
    f = B0*x + B1
    #Initial guess
    Bi = [1.0, 1.0]
    
    altered_trace_num = []
    altered_trace_num_offsets = []
    for val in range(len(OBS_offset)):
        #Run Gauss Newton interpolation solver
        if val == plotval:
            residuals, B = GN_2term_fit(f, Bi, RMS_preFB_amp_interp[val,0:len(OBS[val])], 
                                       np.arange(len(OBS[val])), Strd_Dev[val,0:len(OBS[val])],
                                       Plot=plotval, Verbose=False)
        else:
            residuals, N = GN_2term_fit(f, Bi, RMS_preFB_amp_interp[val,0:len(OBS[val])], 
                                       np.arange(len(OBS[val])), Strd_Dev[val,0:len(OBS[val])],
                                       Plot=plotval, Verbose=False)

        #Figure out constants to scale only the noisy traces by
        for i in range(len(residuals)):
            #If Standard deviation is above the desired threshold, and were within our maximum offset limit
            #Aside: Maximum offset limit as alised previous shot is dipping, doesn't affect underlying traces
            if Strd_Dev[val,i] > SD_Target[val] and abs(OBS_offset[val][i,0]) < maxoffset:
                #Append the trace number for plotting
                if val == plotval:
                    altered_trace_num.append(i)
                    altered_trace_num_offsets.append(OBS_offset[val][i,0])

                #How much we're over the limit
                func = RMS_preFB_amp_interp[val][i]-residuals[i]*scale
                NRMS_scaling[val].append((func[0]/RMS_preFB_amp_interp[val][i])) 
            else:
                NRMS_scaling[val].append(1.0)

        if val == plotval:
            fig_size = plt.rcParams["figure.figsize"]
            fig_size[0] = 16
            fig_size[1] = 8
            plt.rcParams["figure.figsize"] = fig_size
            
            #Find out where the adjusted traces are
            Altered_SD = []
            Altered_Amp = []
            print np.shape(Strd_Dev)
            print np.shape(NRMS_scaling)
            print np.shape(RMS_preFB_amp_interp)
            for i in range(len(OBS_offset[val])):
                if i in altered_trace_num:
                    Altered_SD.append(Strd_Dev[val,i]*NRMS_scaling[val][i])
                    Altered_Amp.append(RMS_preFB_amp_interp[val,i]*NRMS_scaling[val][i])

            #Compute noise floor
            xvals = np.arange(0,len(OBS_offset[val]))
            yvals = xvals*B.item(0)+B.item(1)
            
            #Standard Deviation Plot
            plt.plot(OBS_offset[val][:,0], Strd_Dev[val,0:len(OBS_offset[val])], "ro", label="Original SDs")
            plt.plot(OBS_offset[val][:,0], Strd_Dev[val,0:len(OBS_offset[val])]*NRMS_scaling[val], 
                     "bo", label="Final SDs")
            plt.plot(altered_trace_num_offsets, Altered_SD, "go", label="Updated SDs")
            plt.plot(OBS_offset[val][:,0], np.full((len(OBS_offset[val])),SD_Target[val]), "k-", label="SD Target")
            plt.title("Standard Deviations for NRMS Amplitudes", fontsize=20.0*fontscale)
            plt.xlabel("Offset [km]", fontsize=18.0*fontscale)
            plt.ylabel("Standard Deviation", fontsize=18.0*fontscale)
            mpl.rcParams['xtick.labelsize'] = 16.0*fontscale 
            mpl.rcParams['ytick.labelsize'] = 16.0*fontscale 
            plt.grid()
            plt.legend(fontsize=18*fontscale)
            plt.show()

            #Amplitude plot
            plt.plot(OBS_offset[val][:,0], RMS_preFB_amp_interp[val,0:len(OBS_offset[val])], 
                     "ro", label="Starting RMS Amplitudes")
            plt.plot(OBS_offset[val][:,0], RMS_preFB_amp_interp[val,0:len(OBS_offset[val])]*NRMS_scaling[val], 
                     "bo", label="Final RMS Amplitudes")
            plt.plot(altered_trace_num_offsets, Altered_Amp, "go", label="Updated RMS Amplitudes")
            plt.plot(OBS_offset[val][:,0], np.arange(0,len(OBS_offset[val]))*B.item(0)+B.item(1),
                     "k--", label="Noise Floor")
            plt.title("NRMS Amplitudes", fontsize=20.0*fontscale)
            plt.xlabel("Offset [km]", fontsize=18.0*fontscale)
            plt.ylabel("NRMS Amplitude", fontsize=18.0*fontscale)
            mpl.rcParams['xtick.labelsize'] = 16.0*fontscale 
            mpl.rcParams['ytick.labelsize'] = 16.0*fontscale 
            plt.grid()
            plt.legend(fontsize=18*fontscale)
            plt.show()
            
            plt.plot(OBS_offset[val][:,0], NRMS_scaling[val], "k")
            plt.title("Trace Amplitude Scaling Values", fontsize=20.0*fontscale)
            plt.xlabel("Offset [km]", fontsize=18.0*fontscale)
            plt.ylabel("Amplitude Scaling (%)", fontsize=18.0*fontscale)
            mpl.rcParams['xtick.labelsize'] = 16.0*fontscale 
            mpl.rcParams['ytick.labelsize'] = 16.0*fontscale 
            plt.grid()
            plt.legend(fontsize=18*fontscale)
            plt.show()
            
    return NRMS_scaling

#NOW I GOT TO CLEAN UP THIS GAUSS NEWTON FUNCTION FOR FOUR TERMS TO MAKE IT FOR ONE TERM, C
#THEN I LOOP OVER EVERYTHING MINIMIZING THE MISFIT BETWEEN MEAN AMPLITUDES 0-20 SECONDS FOR MODELED AND OBS DATA
#THIS WILL BE THE SECOND AMPLITDUE SHIFT

#WORKING

def GN_1term_min(amp_real, amp_mdl, rows, c0=1, verbose=True, scalar=None):
    """
    f is the function of which we are trying to minimize.
    c0 is an initial guess at a single constant (in this case) which is going to minimize our function f.
    The observed data here are the observed differences between the two amplitude profiles we
    are tyring to minimze.
    x is just an array containing intiger numbers representing the number of points
    The number of points (rows) are the total number of difference points.
    """

    #Initial Guess - put into a matrix
    C = np.matrix([c0])
    
    #Setup the Jacobian and residual matricies
    Jf = np.zeros((rows,1)) 
    r = np.zeros((rows,1))
    
    #Single iteration linear inversion
    Misfit = 0
    for j in range(rows):
        #Calculate the residual
        r[j,0] = 0 - (amp_real[j]*c0 - amp_mdl[j])
        #Compute the Jacobian (The derivative of the above with respect to c0, without data_obs)
        Jf[j,0] = -amp_real[j]
        #Compute initial misfit
        Misfit += sqrt((0 - amp_real[j]*c0 - amp_mdl[j])**2)

    Jft =  Jf.T

    if verbose==True:
        print "Starting LSQ misfit: "+str(sy.N(Misfit,5)),

    #Update the value of C
    C = C - np.dot(np.dot(np.linalg.inv(np.dot(Jft,Jf)),Jft),r)      
    
    #Re-compute misfit to ensure it decreased
    Misfit = 0
    for j in range(rows):
        Misfit += sqrt((0 - (amp_real[j]*C - amp_mdl[j]))**2)

    if verbose==True:
        print ", Final LSQ misfit: "+str(sy.N(Misfit,5)),
        if scalar == None:
            print ", Constant to best match amplitudes: "+str(sy.N(C,5))
        else:
            print ", Constant to best match amplitudes: "+str(sy.N(C*scalar,5))
        print
        
    if scalar == None:
        return float(C)
    else:
        return float(C)*scalar
    
### RUN JOB ###

import numpy as np
import segyio #This doesn't come standard with python
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.colors as mcolors
import struct
import gc

from scipy import interpolate
from sympy import N
from math import e, sqrt

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
segydir = '/home/chris/0_GlobeClaritas/0_Projects/Ref_EMed_Line1/COMMON/DATA/5_FinalFlow_Sept30th/3_FWIData_RMSNorm/'
name_start = ''
name_end = '_VERT_Final_FWI.sgy'
fileloc = [segydir, name_start, name_end]

##################
"""TOY2DAC DATA"""
##################

#Open toy2dac data_modeling file
dirr = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/6_ForwardModeling/009_Oct29_StartMdl110/" #Path to root directory
model = "data_modeling_qp50_16Hz"
filename = dirr + model #Name of binary file
with open(filename, 'rb') as f:
    data = np.fromfile(f, dtype=np.float32)

#Maximum TOY2DAC geometry file
f = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/000_DATA/2_RealData/17_110StartMdl_THEEND/1_acqui_t0/acqui_100" 

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
maxfold = 1821 #maximum fold, maximum number of traces

#Define variables for model cropping
offset_max = 100
offset_min = 0
bulk_shift = 0 # Not really implimented
min_max = [30, 220]

#Indicies of OBS that have already been cropped out, blank if none: PYTHON INDEXING, OBS number - 1
del_obs = [11, 13]

#Define variables for fourier transform
buff = 100 #retain this many time samples above the first break
div = 2.0 #divide the total number of time samples by this number (saves memory)
LFC = 0 #the Laplace-Fourier Constant, set to zero if it isn't to be used, negative constant applys damping from t=0

#Plot number
plot = False #set as False if to plots are desired.

##################
"""TOY2DAC DATA"""
##################

#Data-related
nchannels = 660 #Number of channels, obsolete if setup_flag != 0 
nsources = 16 #Number of sources 
nfreq = 50000 #Number of discrete freqeuncy groups
relfreq = 3200 #Number of "relevant" frequency groups, i.e. sub 30 Hz
maxfreq = 250 #Maximum frequency in forward modeled data

#Algorithm-related
#Taper design in terms of # of discrete frequencies
taper_design = [400, 1000, 1800] #Save taper design into compact list, [low taper, 100% pass length, high taper]

################
### GEOMETRY ###
################

"""REAL DATA"""

#Read the first-beak picks    
OBS, pick_numbers = read_fbpicks(picks, obs, phase, unc, tt_window, sr)

"""TOY2DAC DATA"""

#Setup geometry from toy2dac acqui file
num_traces, nchannels = setup_geom(f, setup_flag=1, sources=all)

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

#########################
### OUTPUT ACQUI FILE ###
#########################

"""REAL DATA"""

#Output final version of acqui file
acqui, del_idx = gen_toy2dac_acqi(OBS, obs, obs_z, bulk_shift, min_max, del_obs)

######################
### OUTPUT GATHERS ###
######################

"""REAL DATA"""

RD_SG, OBS, max_pick_num, noise_profile = fourier_gathers(OBS, gathers, obsoffsets, obs, 
                                             del_idx, indicesdel_crop, indicesdel, min_idx, 
                                             nsamples, maxfold, buff, div, sr, offset_min,
                                             LFC=LFC, plot=plot)

#Save some memory
del gathers
gc.collect()

#Compute RMS amplitudes for the real data
avg_vals_real = [[] for _ in range(len(RD_SG))]
for i in range(len(RD_SG)):
    obs_data = RD_SG[i]
    for idx in range(len(OBS[i])):
        avg_vals_real[i].append(sqrt(sum(np.square(obs_data[:,idx]))/len(obs_data[:,idx])))

#Save some memory, wont  be able to plot if you delete RD_SG
del obsoffsets,# RD_SG
gc.collect()

"""TOY2DAC DATA DATA"""

#Filter frequency domain data
FD_Comp_SG_taper = taper_FD_data(FD_Comp_SG, taper_design, nsources, nfreq, relfreq, nchannels, sources=all)

#Save some memory
del FD_Comp_SG
gc.collect()

#Inverse fourier transform to time domain
FM_SG = inv_fft(FD_Comp_SG_taper, nsources, nfreq, nchannels, sources=all)

#Compute RMS amplitudes for the TOY2DAC data
avg_vals_model = [[] for _ in range(len(FM_SG.real))]
for i in range(len(FM_SG.real)):
    #Crop to the OBS data
    fm_data = FM_SG.real[i][int(len(FM_SG.real[i][:,0])-(tmax/div/sr)):len(FM_SG.real[i][:,0]), 0:num_traces[i]] 
    fm_data = np.flip(fm_data, axis=0)
    for idx in range(num_traces[i]):
        avg_vals_model[i].append(sqrt(sum(np.square(fm_data[:,idx]))/len(fm_data[:,idx])))
        
#Save some memory, wont be able to plot if youre deleting FM_SG
del FD_Comp_SG_taper,# FM_SG
gc.collect()

#Specify the minimum and maximum offsets for NRMS scaling

from math import sqrt
import matplotlib as mpl
maxoffset, minoffset = 90, 8

#All standard deviations above this value will have their corresponding traces attenuated
# -> Recording SD targets for 30, 220 crop model for each (18) OBS, first 3 omitted
#SD_Target = [0.01, 0.008, 0.01, 0.01, 0.01, 0.008, 0.01, 0.01, 0.01, 0.01, 0.01, 0.002, 0.002, 0.002, 0.002,
#             0.002, 0.002, 0.01] #Old
#SD_Target = [0.01, 0.008, 0.01, 0.01, 0.01, 0.008, 0.01, 0.01, 0.01, 0.01, 0.01, 0.002, 0.002, 0.002, 0.002,
#             0.002, 0.002, 0.01] #Updated
#A little more agressive this time...
SD_Target = [0.0025, 0.006, 0.005, 0.002, 0.002, 0.01, 0.01, 0.01, 
             0.002, 0.002, 0.003, 0.0025, 0.006, 0.005, 0.002, 0.01,] #Updated Oct 30th, Start mdl 110, QP 50

#Bult amplitude shifts of noisy traces are scaled by this percentage
scale = 0.8

###################################
### COMPUTE NRMS AND SD OF NRMS ###
###################################

#NRMS and geometry in offset domain files
RMS_preFB_amp_interp, OBS_offset = compute_NRMS(noise_profile, OBS, minoffset, max_pick_num)

#Standard deviations for each NRMS point
Strd_Dev = compute_SD(RMS_preFB_amp_interp, OBS_offset, max_pick_num)

#########################################
### COMPUTE NRMS TRACE SCALING VALUES ###
#########################################
#TMP

#Compute the NRMS scaling values
NRMS_scaling = NRMS_trace_scale(RMS_preFB_amp_interp, Strd_Dev, SD_Target, OBS_offset, maxoffset, scale,
                     plotval=11, fontscale=1)

"""GENERATE AMPLITUDE CORRECTION FILE"""

#CAUTION: There is one more trace in the real data than the modeled data, I assume this is due to an OBS?
#or am I removing a trace in the modeled data that I shouldn't be?

#Apply the NRMS amplitude correction file
avg_vals_NRMS_adj = [[] for _ in range(len(avg_vals_real))]
#Amplitude correction file
amp_corr_vals = np.zeros((len(num_traces), max(num_traces)))
    
for i in range(len(avg_vals_real)):
    for j in range(len(avg_vals_real[i])):
        avg_vals_NRMS_adj[i].append(avg_vals_real[i][j]*NRMS_scaling[i][j])

#Multiply the output of the amplitude matching by a scalar, this is because the RMS observed data are inherently noiseier
scalar = 1.5
    
print "***START AMPLITUDE OPTIMIZATION***"
print
amp_cnst = []
for it in range(len(num_traces)):
    rows = num_traces[it]
    amp_cnst.append(GN_1term_min(avg_vals_NRMS_adj[it], avg_vals_model[it], rows, verbose=True, scalar=scalar))
    for i in range(num_traces[it]):
        amp_corr_vals[it][i] = NRMS_scaling[it][i]*(amp_cnst[it])

#Save the amplitude constant file with the starting model used to generate these values
#np.savetxt(dirr+"AmpCorr_QP050_sdNov12.txt", amp_corr_vals)

### PLOTTER, DEFINE LATER
"""Need to uncomment RD_SG and FM_SG above, or you'll only be able to plot OBS 16"""

plot_shot = 15
obs_data = RD_SG[plot_shot]
model_data = FM_SG[plot_shot][int(len(FM_SG.real[plot_shot][:,0])-(tmax/div/sr)):len(FM_SG.real[plot_shot][:,0]), 
                              0:num_traces[plot_shot]]
model_data = np.flip(model_data, axis=0)
RD_AMP_CNST = 0.0017043

plt.rcParams["figure.figsize"] = (16,12)
plt.imshow(obs_data[0:50000,0:len(OBS[plot_shot])]*RD_AMP_CNST, aspect=0.1, vmin=-5e-3,
           vmax=5e-3, cmap="hsv")
cbar=plt.colorbar()
plt.plot(range(len(OBS[plot][:,0])),OBS[plot][:,1])
plt.title("Time-domain data, shot: "+str(0))
plt.xlabel('Channel')
plt.ylabel('Time sample')
plt.gca().invert_yaxis()
plt.show()


plt.rcParams["figure.figsize"] = (16,12)
plt.imshow(model_data, aspect=0.1,  vmin=-5e-3,
           vmax=5e-3, cmap="hsv")
plt.title("Frequency domain data, shot: "+str(0))
plt.xlabel('Channel')
plt.ylabel('Time (s)')
cbar = plt.colorbar()

plt.show()