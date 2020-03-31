#RUN ME IN JUPYTER NOTEBOOK, SEPARATE OUT THE FUNCTIONS FROM THE MAIN CODE AT THE END

"""readgenericmodel version 1.0 - Chris Williams
Model must be input as a 3 column text file structured as follows:
x-coordinate z-coordinate velocity
"""

def readgenericmodel(model, size=[260,35], sigma=[0,0], bathy="none", watervel=1500, s=0, wvcu = "off",
                     interpolation="nearest", scaling="no", upscale="no", max_vel="no", min_vel="no"):
    
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

    #Read in entire file
    f = genfromtxt(model)
    f = abs(f)

    #Separate text file
    xcorr = f[:,0]
    zcorr = f[:,1]
    xcorr = np.unique(xcorr)
    zcorr = np.unique(zcorr)
    vels = f[:,2]*1000.0
    #Reshape the velocity mesh
    smesh = np.reshape(vels, (int(max(zcorr))+1, int(max(xcorr))+1), order="F")
    smesh = smesh[0:size[1]+1,0:size[0]+1]
    nx = np.shape(smesh)[1]
    nz = np.shape(smesh)[0]
    #Spacing -TMPPPP!!!!
    x_space = float(size[0])/float(nx-1)
    z_space = float(size[1])/float(nz-1)
    #Coordinates
    nx_max = size[0]
    nx_min = 0
    nz_max = size[1]
    nz_min = 0

    #Output
    print "Original x- and z-nodes: ", nx, ", ", nz
    print "Original x- and z-spacing: ", x_space, ", ", z_space
    print "Original x-coordinate minimum and maximum: ", nx_min, ", ", nx_max
    print "Original z-coordinate minimum and maximum: ", nz_min, ", ", nz_max    
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

    #Remove any model exaggeration
    E = x_space/z_space
    if E != 1.0:
    #if E < 1.0:
        print "Removing a vertical exaggeration of ", E, " from the mesh"
        xarange = np.arange(0,nx,1)
        zarange = np.arange(0,nz,1)
        #Rectangular Bivariate Spline interpolation, kx- ky- degrees to the bivariate spline, s-smoothing
        interp_spline = interpolate.RectBivariateSpline(zarange, xarange, smesh, kx=3, ky=3, s=s)
        zarange_update = np.arange(0, ((nz-1)/E)+1, 1)
        smesh = interp_spline(zarange_update, xarange)
        #Re-calculate model dimensions
        nz, nx = np.shape(smesh)
        x_space = float(nx_max)/float(nx-1)
        z_space = float(nz_max)/float(nz-1)
        print "Number of x-nodes: ", nx
        print "Number of z-nodes: ", nz
        print "Current x- and z- dimension mesh spacing: ", x_space, " and ", z_space
        print "Model exaggeration successfully corrected for at time: {0}s".format(time.time()-tt)
    else:
        print "TOMO2D Mesh is has no exaggeration - proceeding..."
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

    #Check for x-grid structure, raise a flag if mesh is unstructured
    flag_x = ""
    flag_z = ""
    for i in np.arange(nx-1):
        value = abs(xcorr[i+1] - xcorr[i])
        if value != x_space:
            flag_x = "xdimension"

    #Check for z-grid structure, raise a flag if mesh is unstructured
    for k in np.arange(0,int(E*nz)-int(E),int(E)):
        value = abs(zcorr[k+int(E)] - zcorr[k])
        if value != z_space:
            flag_z = "zdimension"

    #TOY2DAC Requires a structured mesh with constant x- and z- grid spacings
    if flag_x == "xdimension" or flag_z == "zdimension":
        print "WARNING: Mesh flagged as unstructured"
        print "Generate a sturctured mesh using " + interpolation + " interpolation."
        #List of every x- an z- coordinate
        points = np.zeros((nz*nx,2))
        for i in np.arange(nx):
            #Account for Vertical Exaggeration
            for k in np.arange(0,int(E*nz),int(E)):
                idx = nz*(i)+(k/int(E))
                points[idx,1] = xcorr[i] - truenx_min
                points[idx,0] = zcorr[k] - truenz_min  

        #Generate the structured mesh based on nodes
        grid_z, grid_x = np.mgrid[nz_min:nz_max+z_space:z_space,nx_min:nx_max+x_space:x_space]

        #Option for interpolation
        if interpolation == 'nearest':
            smesh = interpolate.griddata(points, (smesh.T).flatten(), (grid_z, grid_x), method='nearest')
        if interpolation == 'linear':
            smesh = interpolate.griddata(points, (smesh.T).flatten(), (grid_z, grid_x), method='linear')
        if interpolation == 'cubic':
            smesh = interpolate.griddata(points, (smesh.T).flatten(), (grid_z, grid_x), method='cubic')
        if interpolation != 'nearest' and interpolation != 'linear' and interpolation != 'cubic':
            raise ValueError("Interpolation type not specified, continuing but assuming nearest interpolation")
        x_space = float(nx_max)/float(nx-1)
        z_space = float(nz_max)/float(nz-1)
        print "Number of x-nodes: ", nx
        print "Number of z-nodes: ", nz
        print "Current x- and z- dimension mesh spacing: ", x_space, " and ", z_space 
        print "Structured mesh generation complete at time: {0}s".format(time.time()-tt)
    else: 
        print "TOMO2D Mesh is already structured - proceeding..."
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

    #Upscale model
    if upscale != "no":
        print "Upscaling model by a factor of ", upscale
        xarange = np.arange(0,nx,1)
        zarange = np.arange(0,nz,1)
        xarange_update = np.arange(0, (nx-1.0)+(1.0/upscale), 1.0/upscale)
        zarange_update = np.arange(0, (nz-1.0)+(1.0/upscale), 1.0/upscale)
        #Rectangular Bivariate Spline interpolation, kx- ky- degrees to the bivariate spline, s-smoothing
        interp_spline = interpolate.RectBivariateSpline(zarange, xarange, smesh)
        smesh = interp_spline(zarange_update, xarange_update)
        #Extract Output
        nz, nx = np.shape(smesh)
        x_space = float(nx_max)/float(nx-1)
        z_space = float(nz_max)/float(nz-1)
        print "Mesh x-nodes have been increased to: ", nx
        print "Mesh z-nodes have been increased to: ", nz
        print "Current x- and z- dimension mesh spacing: ", x_space, " and ", z_space 
        print "Model upscaling successfull completed at time: {0}s".format(time.time()-tt)
    else:
        print "Model upscaling flag turned off - proceeding..."
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

    #Optional Model Smoothing
    if sigma != [0,0]:
        print "Begin model smoothing..."
        print "Smoothing in the x- and z-dimension defined as: ", sigma[0], " and ", sigma[1], " respectfully."
        smoothing = [sigma[1], sigma[0]]
        smesh = sp.ndimage.filters.gaussian_filter(smesh, smoothing, mode='reflect')
        print "Model smoothing complete at time: {0}s".format(time.time()-tt)
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    
    #Water Column Velocity Clean-Up
    if sigma != [0,0] or s !=0 or wvcu == "on":
        if bathy == "none":
            print "Beginning water column velocity clean-up..."
            print "FAILURE: A bathymetry file is required for this operation."
        else:
            #Reset water velocities
            smesh_h = np.zeros((nz,nx))
            for b in np.arange(len(bathy)):
                #Number of nodes prior to bathymetry
                numwaternodes = min(range(nz), key=lambda i: abs((i*z_space*1000.0)-bathy[b]))
                waternodes = np.arange(numwaternodes)
                for z in np.arange(numwaternodes):
                    if z in waternodes:
                        #Set the velocity of water, default 1500m/s
                        smesh[z,b] = watervel

            print "Water column velocities tidy-up completed at time: {0}s".format(time.time()-tt)
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

    
    #Check x- and z-spacing
    if x_space == z_space:
        spacing = x_space
        print "Successfully outputted a structured mesh with constant x- and z- spacing of: ", spacing
    else:
        #Warning if x- and z-spacings not exact for more than 3 decimal places (meter scale)
        if round(x_space,3) == round(z_space,3):
            spacing = (x_space + z_space)/2.0
            print "WARNING: x- and z- spacing not exact but equal within 3 decimal places, meter scale."
            print "For dense meshes, this meter-scale error can manifest itself into large kilometer scale errors."
            print "This is managable if you are working with synthetic data."
            print "If working with real data, you may need to reconsider the TOMO2D mesh, or try a different scaling value."
            print "Output a structed mesh with an average x- and z- spacing of: ", spacing
        else:
            raise ValueError("Unsuccessful in outputing a structure mesh within 3 decimal places, x-spacing: ", 
                             x_space, ", z-spacing: ", z_space, 
                             ". Recomend trying a different scaling value. ")
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

    #Impose maximum and minimum velocity conditions
    if max_vel != "no" and min_vel != "no":
        for i in np.arange(nx):
            for k in np.arange(nz):
                if smesh[k,i] > max_vel:
                    smesh[k,i] = max_vel
                if smesh[k,i] < min_vel:
                    smesh[k,i] = min_vel
        print "Maximum and minimum velocity imposed as : ", max_vel, " and ", min_vel, " respectfully"
    if max_vel != "no" and min_vel == "no":
        for i in np.arange(nx):
            for k in np.arange(nz):
                if smesh[k,i] > max_vel:
                    smesh[k,i] = max_vel
        print "Maximum velocity imposed as : ", max_vel
    if min_vel != "no" and max_vel == "no":
        for i in np.arange(nx):
            for k in np.arange(nz):
                if smesh[k,i] < min_vel:
                    smesh[k,i] = min_vel
        print "Minimum velocity imposed as : ", min_vel
    if min_vel == "no" and max_vel == "no":
        print "No maximum or minimum velocity conditions imposed - proceeding..."

    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "ALL JOBS COMPLETED AT TIME: {0}s".format(time.time()-tt)
    
    #Output scaled to meteres for TOY2DAC - doing so prior results in a memory error
    return smesh, [nx, nz], [nx_min*1000.0, nx_max*1000.0], [nz_min*1000.0, nz_max*1000.0], spacing*1000.0

''' Generate Bathymetry 
    Plots OBS locations '''

def generatebathymetry(bath, OBS, shift, nx_max, spacing, plot='true', bathy_1="no"):
    
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "GENERATING BATHYMETRIC PROFILE: "
    if plot =="true":
        print "Will plot bathymetric profile and OBS locations additionally..."
    
    #Correct for OBS locations
    obs = OBS[:,1:3]
    obs[:,0] = obs[:,0] + shift

    #True locations of OBS's
    obs_trueloc = OBS[:,1:3]
    obs_trueloc[:,0] = obs_trueloc[:,0]*1000.0*spacing
    obs_trueloc[:,1] = obs_trueloc[:,1]*1000.0

    #Setup for interpolation
    x_coordinates = np.arange(0,nx_max+(spacing/1000.0),spacing/1000.0)
    x_coordinates_m = x_coordinates*1000.0

    #Bathymetric profile
    bath[:,0] = bath[:,0] + shift
    bath = np.insert(bath, 0, (0,bath[0,1]), axis=0)
    bath = np.insert(bath, len(bath), (nx_max, bath[len(bath)-1,1]), axis=0)

    #Interpolate to match mesh spacing
    tmp_bath = interpolate.interp1d(bath[:,0],bath[:,1])
    bathy_true = tmp_bath(x_coordinates)
    bathy_true = bathy_true*1000.0
    
    if plot == "true":
        #Plot
        fig_size = plt.rcParams["figure.figsize"]
        fig_size[0] = 24
        fig_size[1] = 12
        plt.rcParams["figure.figsize"] = fig_size
        
        plt.plot(x_coordinates_m, bathy_true, label = "True bathymetric profile")
        
        if bathy_1 != "no":
            plt.plot(x_coordinates_m, bathy_1, label = "Bathymetric Profile 1")
        
        plt.plot(obs_trueloc[:,0]/spacing,obs_trueloc[:,1],'bo', label="OBS locations")
        plt.xlabel('x position (m)', fontsize = 24)
        plt.ylabel('Scaled bathymetry (m)', fontsize = 24)
        plt.xticks(fontsize=20) 
        plt.yticks(fontsize=20)
        plt.title("Bathymetry", fontsize = 28)
        plt.gca().invert_yaxis()
        plt.legend(fontsize=20)
        plt.show()

    print "Profile generated at: {0}s".format(time.time()-tt)
    
        
    return bathy_true

''' Accessory functions '''

def sharpenmodel(model, bathy, spacing, alpha=30.0, max_vel=8500.0, min_vel=1450.0):
    
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "GENERATING SHARPENED MODEL: "
    
    #Sharpen a model by approximating the laplacian
    blurred_l = ndimage.gaussian_filter(model, 3)
    filter_blurred_l = ndimage.gaussian_filter(blurred_l, 1)
    alpha = alpha
    sharpened = blurred_l + alpha * (blurred_l - filter_blurred_l)

    #Min/max velocities/correct for waterbottom
    for i in np.arange(len(sharpened[0,:])):
        for j in np.arange(len(sharpened[:,0])):
            depth = spacing*j
            current_z = bathy[i]
            if sharpened[j,i] > max_vel:
                sharpened[j,i] = max_vel
            if sharpened[j,i] < min_vel:
                sharpened[j,i] = min_vel
            if depth < current_z:
                sharpened[j,i] = min_vel
                
    print "Model generated at: {0}s".format(time.time()-tt)
                
    return sharpened

def gardnerrelation(model, bathy, spacing, max_rho=3400.0, min_rho=1000.0):
    
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "GENERATING DENSITY MODEL: "
    
    #Generate Density Model using Gardner's Relation, unit m/s
    rho = np.zeros((len(model[:,0]),len(model[0,:])))
    for j in np.arange(len(model[:,0])):
        for i in np.arange(len(model[0,:])):
            rho[j,i] = 310*(model[j,i])**0.25
            
    #Min/max densities/correct for waterbottom
    for i in np.arange(len(rho[0,:])):
        for j in np.arange(len(rho[:,0])):
            depth = spacing*j
            current_z = bathy[i]
            if rho[j,i] > max_rho:
                rho[j,i] = max_rho
            if rho[j,i] < min_rho:
                rho[j,i] = min_rho
            if depth < current_z:
                rho[j,i] = min_rho
    
    print "Model generated at: {0}s".format(time.time()-tt)
    
    return rho
                
def qgeneration(model, bathy, spacing, max_q=10000.0, min_q=200.0):
    
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "GENERATING DENSITY MODEL: "
    
    #Simple attenuation model, 10,000.0 in water, 200.0 in rock
    q = np.zeros((len(model[:,0]),len(model[0,:])))
            
    #Min/max attenuation/correct for waterbottom
    for i in np.arange(len(q[0,:])):
        for j in np.arange(len(q[:,0])):
            depth = spacing*j
            current_z = bathy[i]
            if depth < current_z:
                q[j,i] = max_q
            if depth >= current_z:
                q[j,i] = min_q

    print "Model generated at: {0}s".format(time.time()-tt)
                
    return q

##################################
### COLORMAP FROM GMT COLORMAP ###
##################################

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

"""We are woking in meters from now on, make sure everything is in the correct unit."""

def cropmodel(model, space, cropx, cropz, bathy="none"):
    
    print "Cropping the model to x-min: ", cropx[0], " x-max: ", cropx[1]
    print "Croppint the model to z-min: ", cropz[0], " z-max: ", cropz[1]
    
    #Min max distance
    minx, maxx = cropx[0], cropx[1]
    minz, maxz = cropz[0], cropz[1]
    
    #Nodes in x and z
    nx = np.shape(model)[1]
    nz = np.shape(model)[0]
    
    #Find out crop indicies in nodes, properly round down for min and up for max
    cropx0, cropx1, cropz0, cropz1 = \
    int(float(minx)/float(space)), int(math.ceil(float(maxx)/float(space)))+1, \
    int(float(minz)/float(space)), int(math.ceil(float(maxz)/float(space)))+1
    
    #Cropped nodes
    nodes_crop = [cropx1-cropx0, cropz1-cropz0]
    
    #Crop the model
    out = model[cropz0:cropz1, cropx0:cropx1]
    
    print "New number of nodes in the x- and z- dimension: ", nodes_crop[0], " and ", nodes_crop[1]
    
    print "README: "
    print cropx0, cropx1
    
    if bathy != "none":
        print "Producing a cropped bathymetric profile..."
        print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        #Crop the bathymetric profile
        outbathy = bathy[cropx0:cropx1]
        return out, outbathy, nodes_crop
    
    else:
        print "Assuming a cropped bathymetric profile has already been output..."
        print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        return out, nodes_crop
    
#THREE TERM GAUSS NEWTON BEST FIT FUNCITON
#CHRIS WILLIAMS - 2020

import sympy as sy
from math import sqrt
from numpy.linalg import inv

#Define the Gauss Newton Interpolation function

def GN_3Term_fit(f, Bi, data_obs, x = sy.Symbol('x'), B0 = sy.Symbol('B0'), B1 = sy.Symbol('B1'), 
                 B2 = sy.Symbol('B2'), B3 = sy.Symbol('B3')):

    #Lengths
    S, num_pts = np.arange(len(data_obs)), len(data_obs)
    
    #Initial Guess - put into a matrix
    b0, b1, b2, b3 = Bi[0], Bi[1], Bi[2], Bi[3]
    B = np.matrix([[b0],[b1],[b2],[b3]])

    #Serup the Jacobian and residual matricies, first is for the control, second is weighted
    rows, cols = num_pts, 4
    Jf, r = np.zeros((rows,cols)), np.zeros((rows,1)) 
    
    #Define the function, and setup the partial derivatives to be filled into the Jacobian matrix
    def model(f, b0, b1, b2, b3, xi, 
              x = sy.Symbol('x'), B0 = sy.Symbol('B0'), 
              B1 = sy.Symbol('B1'), B2 = sy.Symbol('B2'),
              B3 = sy.Symbol('B3')):
        return f.subs(x,xi).subs(B0, b0).subs(B1, b1).subs(B2, b2).subs(B3, b3)

    def partialDerB0(f, b0, b1, b2, b3, xi, 
              x = sy.Symbol('x'), B0 = sy.Symbol('B0'), 
              B1 = sy.Symbol('B1'), B2 = sy.Symbol('B2'),
              B3 = sy.Symbol('B3')):
        return - f.diff(B0,1).subs(x,xi).subs(B0, b0).subs(B1, b1).subs(B2, b2).subs(B3, b3)

    def partialDerB1(f, b0, b1, b2, b3, xi, 
              x = sy.Symbol('x'), B0 = sy.Symbol('B0'), 
              B1 = sy.Symbol('B1'), B2 = sy.Symbol('B2'),
              B3 = sy.Symbol('B3')):
        return - f.diff(B1,1).subs(x,xi).subs(B0, b0).subs(B1, b1).subs(B2, b2).subs(B3, b3)

    def partialDerB2(f, b0, b1, b2, b3, xi, 
              x = sy.Symbol('x'), B0 = sy.Symbol('B0'), 
              B1 = sy.Symbol('B1'), B2 = sy.Symbol('B2'),
              B3 = sy.Symbol('B3')):
        return - f.diff(B2,1).subs(x,xi).subs(B0, b0).subs(B1, b1).subs(B2, b2).subs(B3, b3)

    def partialDerB3(f, b0, b1, b2, b3, xi, 
              x = sy.Symbol('x'), B0 = sy.Symbol('B0'), 
              B1 = sy.Symbol('B1'), B2 = sy.Symbol('B2'),
              B3 = sy.Symbol('B3')):
        return - f.diff(B3,1).subs(x,xi).subs(B0, b0).subs(B1, b1).subs(B2, b2).subs(B3, b3)

    def residual(f, xi, observed, b0, b1, b2, b3):
        return (observed - f.subs(x, xi).subs(B0, b0).subs(B1, b1).subs(B2, b2).subs(B3, b3))
    
    def LSQMisfit(f, xi, observed, b0, b1, b2, b3):
        return sqrt((observed - f.subs(x, xi).subs(B0, b0).subs(B1, b1).subs(B2, b2).subs(B3, b3))**2)
    
    #Fill the Jacobian and residual matricies, solve for the updated fit
    sumOfResid = 0
    TotalLSQMisfit = 0
    for j in range(0,rows,1):
        r[j,0] = residual(f, S[j], data_obs[j], B[0], B[1], B[2], B[3])
        Jf[j,0] = partialDerB0(f, b0, b1, b2, b3, S[j])
        Jf[j,1] = partialDerB1(f, b0, b1, b2, b3, S[j])
        Jf[j,2] = partialDerB2(f, b0, b1, b2, b3, S[j])
        Jf[j,3] = partialDerB3(f, b0, b1, b2, b3, S[j])  
        sumOfResid += (r[j,0] * r[j,0])
        TotalLSQMisfit += LSQMisfit(f, S[j], data_obs[j], B[0], B[1], B[2], B[3])

    #Solve for the updated points
    TotalLSQMisfit_Prior = TotalLSQMisfit
    Jft =  Jf.T
    B = B - np.dot(np.dot(inv(np.dot(Jft,Jf)),Jft),r)

    #Updated Misfit
    TotalLSQMisfit = 0
    for j in range(0,rows,1):
        TotalLSQMisfit += LSQMisfit(f, S[j], data_obs[j], B[0], B[1], B[2], B[3])
    TotalLSQMisfit_Post=TotalLSQMisfit
        
    assert TotalLSQMisfit_Prior > TotalLSQMisfit_Post

    yval = np.zeros((len(S),1))
    for k in range(0, len(S), 1):
        yval[k,0] = model(f, sy.N(B[0],5), sy.N(B[1], 5), sy.N(B[2], 5), sy.N(B[3],5), S[k])
        
    return yval

"""readtomo2dmodel version 1.0 - Chris Williams
~~~ README ~~~
- This code is meant to read in velocity files from TOMO2D (seismic tomography code) and reformat
    the models to be used in TOY2DAC (full waveform inversion - FWI code)
FLAGS
- "model" must be the full path to the TOMO2D velocity model
- nx and nz must be the number of mesh nodes in the x- and z-dimensions respectively on the TOMO2D model
- "interpolation" is used if the TOMO2D mesh is unstructured as TOY2DAC requires a structured mesh.
    It has three possible options: "nearest", "linear", or "cubic", "nearest" is default.
- "scaling" is optional and may be used if a different model size is desired. The TOMO2D 
    mesh will be automatically scaled in the x- and z-dimension based on the "scaling" value.
    If "scaling" is desired, set scaling=value (as a percentage), if not, set scaling="no". 
- "upscale" is another optional flag and is used to generate a coarser grid. This is useful as FWI requires
    a finer node spacing than tomography. "upscale" must be set equal to the FLOAT that you wish to upscale 
    the model equally in the x- and z-dimension, i.e. upscale=FLOAT. If you do not wish to upscale the model, 
    set upscale="no".
-"max_vel" and "min_vel" set the maximum and minimum velocities when regridding the TOMO2D mesh to TOY2DAC format.
    If this is not desired set max_vel="no" and min_vel="no"
FURTHER NOTES
- Upscaling flag is only meant to increase the mesh density, this will generally be the case when 
    moving from TOMO2D to TOY2DAC as FWI requires a denser mesh than tomography
-"origin" sets the minimum value for the x- and z-dimensions. The origin will always be set to zero. This would
    be a problem for land acquisitions...
-When using "scaling" try to use a "nice" number that will give easy to work with values. Critical that x- and
    z-spacing is equal. 

"""

def readtomo2dmodel(model, nx, nz, detrend=False, sigma=[0,0], bathy="none", watervel=1500, interpolation="nearest", 
                    scaling="no", upscale="no", max_vel="no", min_vel="no", plt_mdl=False):

    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

    """ PHASE 1: EXTRACT INFORMATION FROM TOMO2D FILE"""
    
    #Skip first four rows as they contain non-velocity related information
    vels = genfromtxt(model,skip_header=4)
    #Read in the other four rows, nodes, x-coordinates, bathymetry, z-coordinates
    nodes = genfromtxt(model, skip_footer = nx+3)
    xcorr = genfromtxt(model, skip_header = 1, skip_footer = nx+2)
    bathy = genfromtxt(model, skip_header = 2, skip_footer = nx+1)
    zcorr = genfromtxt(model, skip_header = 3, skip_footer = nx)
    #Establish maximimum and minimum dimensions of the model - convert to meters
    nx_max = xcorr[len(xcorr)-1]
    nx_min = xcorr[0]
    nz_max = zcorr[len(zcorr)-1] #NOT WORKING WITH A POORLY BEHAIVED TOMO2D MODEL
    nz_min = zcorr[0]
    truenx_max = (nx_max - nx_min)
    truenz_max = (nz_max - nz_min)
    truenx_min = nx_min
    truenz_min = nz_min
    x_space = float(truenx_max)/float(nx-1)
    z_space = float(truenz_max)/float(nz-1)
    
    print "Original x- and z-nodes: ", nx, ", ", nz
    print "Original x- and z-spacing: ", x_space, ", ", z_space
    print "Original x-coordinate minimum and maximum: ", nx_min, ", ", nx_max
    print "Original z-coordinate minimum and maximum: ", nz_min, ", ", nz_max    
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    
    #Convert velocities to m/s
    smesh = np.zeros((nz,nx))
    for b in np.arange(nx):
        for z in np.arange(nz):
            smesh[z,b] = vels[b,z] * 1000.0
    
    #Plot Origional TOMO2D Model
    if plt_mdl != False:
        print "Velocity hanging from bathymtery not implimented at this point."
        if x_space != z_space:
            print "WARNING: x and z spacing not equal, z-axis label will be incorrect for this plot."
        plot_onemdl(smesh, [nx,nz], x_space, vrange=[1500,8000], unit="[m]", fontscale=3.0, size=[72,36], 
                    cmap=plt_mdl[0])
    
    #Shift model to zero eleminating negative values for simplicity
    nz_min = 0
    nx_min = 0
    print
    print "Minimum x-dimension coordinate has been changed from: ", truenx_min, " to: ", nx_min
    print "Minimum z-dimension coordinate has been changed from: ", truenz_min, " to: ", nz_min
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    
    """ PHASE 2: Remove Any Model Exaggeration """
    
    #Remove any model exaggeration
    E = x_space/z_space
    if E != 1.0:
    #if E < 1.0:
        print "Removing a vertical exaggeration of ", E, " from the mesh"
        xarange = np.arange(0,nx,1)
        zarange = np.arange(0,nz,1)
        #Rectangular Bivariate Spline interpolation, kx- ky- degrees to the bivariate spline, s-smoothing
        interp_spline = interpolate.RectBivariateSpline(zarange, xarange, smesh, kx=3, ky=3, s=0)
        zarange_update = np.arange(0, nz, E)
        smesh = interp_spline(zarange_update, xarange)
        #Re-calculate model dimensions
        nz, nx = np.shape(smesh)
        x_space = float(truenx_max)/float(nx-1)
        z_space = float(truenz_max)/float(nz-1)
        print "Number of x-nodes: ", nx
        print "Number of z-nodes: ", nz
        print "Current x- and z- dimension mesh spacing: ", x_space, " and ", z_space
        print "Model exaggeration successfully corrected for at time: {0}s".format(time.time()-tt)
    else:
        print "TOMO2D Mesh is has no exaggeration - proceeding..."
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    
    """ Phase 3: Ensure Grid Mesh is Structured """
    
    #Check for x-grid structure, raise a flag if mesh is unstructured
    flag_x = ""
    flag_z = ""
    for i in np.arange(nx-1):
        value = abs(xcorr[i+1] - xcorr[i])
        if value != x_space:
            flag_x = "xdimension"
            
    #Check for z-grid structure, raise a flag if mesh is unstructured
    for k in np.arange(0,int(E*nz)-int(E),int(E)):
        value = abs(zcorr[k+int(E)] - zcorr[k])
        if value != z_space:
            flag_z = "zdimension"
    
    #I need to interpolate the previous z coordinates to the z coordinates without vertical exaggeration
    tmp_zcorr = interpolate.interp1d(range(len(zcorr)),zcorr)
    zcorr = tmp_zcorr(np.arange(0,len(zcorr),float(len(zcorr))/float(nz)))
    
    #TOY2DAC Requires a structured mesh with constant x- and z- grid spacings
    if flag_x == "xdimension" or flag_z == "zdimension":
        print "WARNING: Mesh flagged as unstructured"
        print "Generate a sturctured mesh using " + interpolation + " interpolation."
        #List of every x- an z- coordinate
        points = np.zeros((nz*nx,2))
        for i in np.arange(nx):
            #Account for Vertical Exaggeration
            for k in np.arange(nz):
                idx = (nz)*(i)+(k/int(E)) #MINUS 1 on nz???
                points[idx,1] = xcorr[i] - truenx_min
                points[idx,0] = zcorr[k] - truenz_min  
        
        #Generate the structured mesh based on nodes
        grid_z, grid_x = np.mgrid[nz_min:truenz_max+z_space:z_space,nx_min:truenx_max+x_space:x_space]
        
        #Option for interpolation
        if interpolation == 'nearest':
            smesh = interpolate.griddata(points, (smesh.T).flatten(), (grid_z, grid_x), method='nearest')
        if interpolation == 'linear':
            smesh = interpolate.griddata(points, (smesh.T).flatten(), (grid_z, grid_x), method='linear')
        if interpolation == 'cubic':
            smesh = interpolate.griddata(points, (smesh.T).flatten(), (grid_z, grid_x), method='cubic')
        if interpolation != 'nearest' and interpolation != 'linear' and interpolation != 'cubic':
            raise ValueError("Interpolation type not specified, continuing but assuming nearest interpolation")
        x_space = float(truenx_max)/float(nx-1)
        z_space = float(truenz_max)/float(nz-1)
        print "Number of x-nodes: ", nx
        print "Number of z-nodes: ", nz
        print "Current x- and z- dimension mesh spacing: ", x_space, " and ", z_space 
        print "Structured mesh generation complete at time: {0}s".format(time.time()-tt)
    else: 
        print "TOMO2D Mesh is already structured - proceeding..."
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

    """Phase 3.5 (optional): Detrend velocity mesh"""

    if detrend==True:
        #Terms for 3-term GN fit
        x = sy.Symbol('x')
        B0 = sy.Symbol('B0')
        B1 = sy.Symbol('B1')
        B2 = sy.Symbol('B2')
        B3 = sy.Symbol('B3')
        f = B0*x**3 + B1*x**2 + B2*x + B3 #3 term polynomial
        Bi = [1.0, 1.0, 1.0, 1.0] #Initial guess

        print "Detrending column(s)..."
        #loop through all columns

        for x in range(len(smesh[0,:])):
            print x,
            
            #Obtain the detrended velocities
            v_prof_detrended = GN_3Term_fit(f, Bi, smesh[:,x])

            #Construct that given column of the detrended velocity model
            smesh[:,x] = v_prof_detrended.T

        print "MODEL DETRENDED AT TIME: {0}s".format(time.time()-tt)
        print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

    """ Phase 4: Upscale Mesh """
    # Lets try this griddata thing instead of RectBivariateSpline
    if upscale != "no":
        grid_x, grid_z = np.meshgrid(np.arange(0,(nx-1.0)+(1.0/upscale),1.0/upscale), 
                                     np.arange(0,(nz-1.0)+(1.0/upscale),1.0/upscale))
        points = np.zeros((nx*nz,2))
        values = np.zeros((nx*nz,1))
        for i in range(nx):
            for j in range(nz):
                points[i*nz+j,0] = i
                points[i*nz+j,1] = j
                values[i*nz+j,0] = smesh[j,i]
        grid_z = interpolate.griddata(points, values, (grid_x, grid_z), method='cubic') #NEAREST
        smesh = grid_z.T[0].T    
        #Extract Output
        nz, nx = np.shape(smesh)
        x_space = float(truenx_max)/float(nx-1)
        z_space = float(truenz_max)/float(nz-1)
        print "Mesh x-nodes have been increased to: ", nx
        print "Mesh z-nodes have been increased to: ", nz
        print "Current x- and z- dimension mesh spacing: ", x_space, " and ", z_space 
        print "Model upscaling successfull completed at time: {0}s".format(time.time()-tt)        
    else:
        print "Model upscaling flag turned off - proceeding..."
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    
    """ Phase 5: Optional Smoothing and QC """
    
    #Optional Model Smoothing
    if sigma != [0,0]:
        print "Begin model smoothing..."
        print "Smoothing in the x- and z-dimension defined as: ", sigma[0], " and ", sigma[1], " respectfully."
        smoothing = [sigma[1], sigma[0]]
        smesh = sp.ndimage.filters.gaussian_filter(smesh, smoothing, mode='reflect')
        print "Model smoothing complete at time: {0}s".format(time.time()-tt)
        #We don't need to tidy up for bathymetry as hanging bathymetry is implimented later
    else:
        print "Model smoothing flag turned off - proceeding..."
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        
    #Check x- and z-spacing
    if x_space == z_space:
        spacing = x_space
        print "Successfully outputted a structured mesh with constant x- and z- spacing of: ", spacing
    else:
        #Warning if x- and z-spacings not exact for more than 3 decimal places (meter scale)
        if round(x_space,3) == round(z_space,3):
            spacing = (x_space + z_space)/2.0
            print "WARNING: x- and z- spacing not exact but equal within 3 decimal places, meter scale."
            print "For dense meshes, this meter-scale error can manifest itself into large kilometer scale errors."
            print "This is managable if you are working with synthetic data."
            print "If working with real data, you may need to reconsider the TOMO2D mesh, or try a different scaling value."
            print "Output a structed mesh with an average x- and z- spacing of: ", spacing
        else:
            raise ValueError("Unsuccessful in outputing a structure mesh within 3 decimal places, x-spacing: ", 
                             x_space, ", z-spacing: ", z_space, 
                             ". Recomend trying a different scaling value. ")
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    
    """ Phase 6: Hanging Bathymetery Implimentation """
    
    #Plot plot velocity model prior to hanging bathym implimentation - QC
    if plt_mdl != False:
        print "Plotting the inverse of slowness squared following model upscaling"
        plot_onemdl(smesh, [nx,nz], x_space, vrange=[1500,8000], unit="[m]", fontscale=3.0, size=[72,36], 
                    cmap=plt_mdl[0])
    
    #Account for hanging bathymetry last - take advantage of a higher resolution model
    xcorr_update = np.arange(0, truenx_max+spacing, spacing) #REMOVED SPACING, NOW I ADDED SPACING BACK???
    zcorr_update = np.arange(0, truenz_max+spacing, spacing)

    tmp_bathy = interpolate.interp1d(xcorr-truenx_min,bathy)
    bathy_update = tmp_bathy(xcorr_update)
    
    #Velocity model hanging from bathymetry
    smesh_h = np.zeros((nz,nx))
    for b in np.arange(len(bathy_update)):
        #Number of nodes prior to hanging velocity mesh - assume mesh is unstrutured
        numwaternodes = min(range(len(zcorr_update)), key=lambda i: abs(zcorr_update[i]-bathy_update[b]))
        waternodes = np.arange(numwaternodes)
        for z in np.arange(len(zcorr_update)):
            if z in waternodes:
                #Set the velocity of water
                smesh_h[z,b] = nodes[2]*1000.0
            else:
                #TOMO2D operates in km/s, TOY2DAC operates in m/s
                smesh_h[z,b] = smesh[z-numwaternodes,b]
                
    """ PHASE 7: Velocity Bounds """
    
    #Impose maximum and minimum velocity conditions
    if max_vel != "no" and min_vel != "no":
        for i in np.arange(nx):
            for k in np.arange(nz):
                if smesh_h[k,i] > max_vel:
                    smesh_h[k,i] = max_vel
                if smesh_h[k,i] < min_vel:
                    smesh_h[k,i] = min_vel
        print "Maximum and minimum velocity imposed as : ", max_vel, " and ", min_vel, " respectfully"
    if max_vel != "no" and min_vel == "no":
        for i in np.arange(nx):
            for k in np.arange(nz):
                if smesh_h[k,i] > max_vel:
                    smesh_h[k,i] = max_vel
        print "Maximum velocity imposed as : ", max_vel
    if min_vel != "no" and max_vel == "no":
        for i in np.arange(nx):
            for k in np.arange(nz):
                if smesh_h[k,i] < min_vel:
                    smesh_h[k,i] = min_vel
        print "Minimum velocity imposed as : ", min_vel
    if min_vel == "no" and max_vel == "no":
        print "No maximum or minimum velocity conditions imposed - proceeding..."
    
    #Final Model QC, Now with Bathymetry implimented
    if plt_mdl != False:
        print "Plotting the inverse of slowness squared following model upscaling"
        plot_onemdl(smesh_h, [nx,nz], x_space, vrange=[1500,8000], unit="[m]", fontscale=3.0, size=[72,36], 
                    cmap=plt_mdl[0])
    
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "ALL JOBS COMPLETED AT TIME: {0}s".format(time.time()-tt)

    #Output scaled to meteres for TOY2DAC - doing so prior results in a memory error
    return smesh_h, bathy_update*1000.0, [nx, nz], [nx_min, nx_max*1000.0], [nx_min*1000.0, nz_max*1000.0], spacing*1000.0

#########################
### PLOT SINGLE MODEL ###
#########################

def plot_onemdl(smesh_1, nodes, spacing, savepath=None, savename=None, vrange="true", 
                          unit="m", fontscale=1, shift=[0,0], letterlbl=None, size=[72,36], cmap='jet',
                         bathy=None, src_loc=False):
    
    """ 
    This code is developed to gereate a high quality image of the starting model and true model for
    a synthetic inversion.
    -> Modified from the original as we now ha
    ve TOMO2D data, 
    """

    #########################
    ### READ BINARY FILES ###
    #########################

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
    nz_min = nodes[1] * spacing
    nz_max = 0 * spacing

    #Define Figure Size
    size = [72,36]
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = size[0]
    fig_size[1] = size[1]
    plt.rcParams["figure.figsize"] = fig_size

    #Colorbar range, define limits based on starting or true model, or user specified
    if vrange == "true":
        elev_min, elev_max = min(smesh_1.flatten()), max(smesh_1.flatten())
        print "Auto Colorbar, min/max: ", elev_min, "/", elev_max
    else:
        elev_min, elev_max = vrange[0], vrange[1]

    #Setup subplots
    fig = plt.figure() 
    ax1 = fig.add_subplot(111)

    #Plot
    if type(bathy) != bool:
        ax1.plot(np.arange(len(bathy_mesh))*spacing, bathy_mesh, color='black', linewidth=6)
    if type(src_loc) != bool:
        ax1.scatter(src_loc[:,0], src_loc[:,1], s=100*fontscale, color="black", marker="v")
    im = ax1.imshow(smesh_1, extent=[nx_min,nx_max,nz_min,nz_max], vmin=elev_min, vmax=elev_max, cmap=cmap,
                   aspect=2.0)
    ax1.set_title('Attenuation Model',fontweight="bold", size=20.0*fontscale) # Title
    ax1.set_ylabel('Z-Position '+unit, fontsize = 18.0*fontscale) # Y label
    ax1.set_xlabel('X-Position '+unit, fontsize = 18.0*fontscale) # X label
    ax1.tick_params(axis='both', which='major', labelsize=16.0*fontscale)
    ax1.tick_params(axis='both', which='minor', labelsize=16.0*fontscale)

    #Colorbar
    #ax_cbar = fig.add_axes([0.81+shift[1], 0.3, 0.03, 0.40])
    #cbar = fig.colorbar(im, cax=ax_cbar, orientation='vertical')
    #cbar.ax.tick_params(labelsize=16.0*fontscale)
    #cbar.set_label("Attenuation", fontsize=18.0*fontscale)
    
    if letterlbl != None:
        #Add text
        fig.text(0.086+shift[0], 0.86, letterlbl[0], fontweight="bold", fontsize=18*fontscale)
        
    if savepath != None and savename != None:
        plt.gcf().savefig(savepath+savename+".png", bbox_inches = 'tight', pad_inches = 0.1)

    plt.show()

"""~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BREAK ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""
    
##################################
### TOMO2D STARTING MODEL REAL ###
##################################

import time
tt = time.time()

""" Model Generation """

#TOMO2D Starting Model
model_starting = "/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/110_pick9.final.zerophase_108run_Oct29/0_inverted_models/09_smesh"

#Node numbers for TOMO2D Model
nx = 951
nz = 151

#Node numbers for cropped models
cropx = [12922.0, 182692.0]
cropz = [0, 30000.0]

""" Plotting Files """

# Make a colormap from GMT file
cpt = gmtColormap('/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/rgb.kw.Nov09.cpt')#rgb.celw.Jan26.cpt')#
cpt_convert = colors.LinearSegmentedColormap('cpt', cpt)

#Bathymetry file
bathy = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/000_DATA/2_RealData/17_110StartMdl_THEEND/0_models/model_50m/fbathy"

#Import x,z OBS locations
dirr = "/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/00_TOMO2DCONVERTED/2_EmedSynthetic/0_FinalSourceGeoms/"
f = "Synthetic_SourceGeom_Final_True"
src_loc = np.genfromtxt(dirr+f)

#Save as plot model parameters, turning the optional model plot on
plt_mdl = [cpt_convert, bathy, src_loc]

######################
### STARTING MODEL ###
######################

#TOMO2D Model; better results if the model is structured
vp_s, bathy_tomo, nodes_s, xdim_s, zdim_s, spacing_s = readtomo2dmodel(model_starting, nx, nz, detrend=False,
                                                  interpolation="nearest", scaling = "no", upscale=4.0, 
                                                  max_vel="no", min_vel="no", plt_mdl=False)

print np.shape(bathy_tomo)

#Crop the generic model
vp_s, croppedbathy, nodes_crop = cropmodel(vp_s, spacing_s, cropx, cropz, bathy=bathy_tomo)

plt.plot(croppedbathy)
plt.show()

plot_onemdl(vp_s, nodes_crop, spacing_s, vrange=[0,9000], unit="[km]", fontscale=3.0, size=[72,36], 
            cmap=cpt_convert, bathy=croppedbathy, src_loc=src_loc)

#####################
### DENSITY MODEL ###
#####################

rho = gardnerrelation(vp_s, croppedbathy, spacing_s, max_rho=5000.0, min_rho=1000.0)

plot_onemdl(rho, nodes_crop, spacing_s, vrange=[1000,4000], unit="[km]", fontscale=3.0, size=[72,36], 
            cmap="jet", bathy=bathy, src_loc=src_loc)

#########################
### ATTENUATION MODEL ###
#########################
q = qgeneration(vp_s, croppedbathy, spacing_s, max_q=10000.0, min_q=50.0)

plot_onemdl(q, nodes_crop, spacing_s, vrange=[50,10000], unit="[km]", fontscale=3.0, size=[72,36], 
            cmap="jet", bathy=bathy, src_loc=src_loc)

"""~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BREAK ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""

####################################
### DETRENDED FWI STARTING MODEL ###
####################################

import time
tt = time.time()

#TOMO2D Starting Model
model_starting = "/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/110_pick9.final.zerophase_108run_Oct29/0_inverted_models/09_smesh"

#Node numbers for TOMO2D Model
nx = 951
nz = 151

#Node numbers for cropped models
cropx = [0, 190000.0]
cropz = [0, 30000.0]

""" Plotting Files """

# Make a colormap from GMT file
cpt = gmtColormap('/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/rgb.kw.Nov09.cpt')#rgb.celw.Jan26.cpt')#
cpt_convert = colors.LinearSegmentedColormap('cpt', cpt)

#Bathymetry file
bathy = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/000_DATA/2_RealData/17_110StartMdl_THEEND/0_models/model_50m/fbathy"

#Import x,z OBS locations
dirr = "/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/00_TOMO2DCONVERTED/2_EmedSynthetic/0_FinalSourceGeoms/"
f = "Synthetic_SourceGeom_Final_True"
src_loc = np.genfromtxt(dirr+f)

#Save as plot model parameters, turning the optional model plot on
plt_mdl = [cpt_convert, bathy, src_loc]

######################
### STARTING MODEL ###
######################

#TOMO2D Model; better results if the model is structured
vp_s_detrended, bathy_tomo, nodes_s, xdim_s, zdim_s, spacing_s = readtomo2dmodel(model_starting, nx, nz,  
                                                                                 sigma=[10,10],
                                                  interpolation="nearest", scaling = "no", upscale=4.0, 
                                                  max_vel="no", min_vel="no", plt_mdl=False)

plot_onemdl(vp_s_detrended, nodes_crop, spacing_s, vrange=[0,9000], unit="[km]", fontscale=3.0, size=[72,36], 
            cmap=cpt_convert, bathy=bathy, src_loc=src_loc)

"""~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BREAK ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""

##############################
### SYNTHETIC EMED VERSION ###
##############################

### RUN CROPPED VERSION ###
import time
tt = time.time()

### ~~~ Bathymetry ~~~ ###

#Import true bathymetry
f = "/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/66_pick9.2_65run_Aug02/_startmdl/kw_data/bathy_kw.dat"
bath = np.genfromtxt(f)
bathy = bath[:,1]
#Import OBS coordinates
f = "/home/chris/Documents/1_data/EasternMed/Chris_Williams/locations/OBS_Bathy.dat"
OBS = np.genfromtxt(f)

### ~~~ Models ~~~ ###

#True Model
model = "/home/chris/Documents/1_data/EasternMed/mdls/new_velocity_model_1kmby1km_Chris.txt"
#Starting Model (For synthetic generation, Keep as Full Model then Crop)
model_starting = "/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/66_pick9.2_65run_Aug02/0_inverted_models/09_smesh_2x_1z"

#Node numbers for TOMO2D Model
nx = 801 
nz = 141

#FOR ALL MODELS: BE CAREFUL HOW YOU CROP
#Node numbers for MODELS CAREFUL
cropx = [60000, 250000.0]
cropz = [0, 30000.0]

#cropxb = [180000, 250000.0]
#cropzb = [0, 30000.0]

#THERE IS SOMETHING FUNDAMENTALLY WEIRD GOING ON HERE WITH THE BAHTYMETRY FILE, IM SLIGHTLY OFF.
#I'M TOO FUCKING FED UP WITH IT NOW TO FIX IT.

######################
### STARTING MODEL ###
######################

#TOMO2D Model
#If the TOMO2D model is unstructured and ugly, its better to upscale it using edit_smesh rather than in this code.
#vp_s, bathy_tomo, nodes_s, xdim_s, zdim_s, spacing_s = readtomo2dmodel(model_starting, nx, nz,
#                                                  interpolation="cubic", scaling = "no", upscale=5.0,#6.25, 
#                                                  max_vel=10000.0, min_vel=1450.0)

#Fully upscale the bathymetry file
xcorr_b = np.arange(-30, 230+(260./300.), (260./300.))
xcorr_b_update = np.arange(-30, 230, 0.05)
tmp_bathy = interpolate.interp1d(xcorr_b, bathy*1000.0)
bathyfull_upscale = tmp_bathy(xcorr_b_update)
bathyfull_upscale = np.append(bathyfull_upscale, bathyfull_upscale[len(bathyfull_upscale)-1])

#########################################
### MY BATHYMETRY FUDGE FOR SYNTHETIC ###
#########################################

#MAKE OBS POINTS TO INTERPOLATE OVER
#obs = OBS[:,1:3]
#obs = np.insert(obs, 0, (-30,obs[0,1]), axis=0) #Insert start point
#obs = np.insert(obs, 0, (-15,obs[0,1]), axis=0) #Insert start point
#obs = np.insert(obs, len(obs), (230,obs[len(obs)-1,1]), axis=0) #Insert end point
#obs = np.insert(obs, len(obs), (222,obs[len(obs)-1,1]), axis=0) #Insert end point
#X-COORDINATES & INTERPOLATE
#xcorr = np.arange(-30, 230, 0.05) #SPACING & MODEL START STOP
#tmp_bath = interpolate.interp1d(obs[:,0],obs[:,1]*1000.0,kind="quadratic")
#BATHY_FUDGE = tmp_bath(xcorr)
#BATHY_FUDGE = np.append(BATHY_FUDGE, BATHY_FUDGE[len(BATHY_FUDGE)-1])

#Generic Model
vp_s, nodes_s, xdim_s, zdim_s, spacing_s = readgenericmodel(model, size=[260,35], sigma=[40,90], bathy=BATHY_FUDGE, 
                                                  interpolation="nearest", scaling="no", upscale=20.0,#16.0, 
                                                  max_vel=10000.0, min_vel=1450.0)

#Crop the generic model
vp_s, croppedbathy, nodes_crop = cropmodel(vp_s, spacing_s, cropx, cropz, bathy=BATHY_FUDGE)
#Crop the bathy model
#dummy, croppedbathy, nodes_crop = cropmodel(vp_s, spacing_s, cropxb, cropzb, bathy=bathyfull_upscale)

#Plot the model
plotmodel(vp_s, nodes_crop, cropx, cropz, title="Starting Velocity Model", unit="velocity")

##################
### TRUE MODEL ###
##################
vp, nodes, xdim, zdim, spacing = readgenericmodel(model, size=[260,35], bathy=BATHY_FUDGE, s=0, wvcu = "on",
                                                  interpolation="nearest", scaling="no", upscale=20.0,#16.0,
                                                  max_vel=10000.0, min_vel=1450.0)

#Crop the model
vp_crop, nodes_crop = cropmodel(vp, spacing, cropx, cropz)

#Plot the model
plotmodel(vp_crop, nodes_crop, cropx, cropz, title="True Velocity", unit="velocity")

##################
### BATHYMETRY ###
##################
#bathy = generatebathymetry(bath, OBS, 60.0, 250.0, spacing_s, plot='true', bathy_1=croppedbathy)

#######################
### SHARPENED MODEL ###
#######################
#sharpened = sharpenmodel(vp, bathy_tomo, spacing, alpha = 120.0, max_vel=8500.0, min_vel=1450.0)

#Plot the model
#plotmodel(sharpened, nodes_crop, cropx, cropz, title="Sharpened True Velocity", unit="velocity")

#####################
### DENSITY MODEL ###
#####################

rho = gardnerrelation(vp_s, croppedbathy, spacing_s, max_rho=5000.0, min_rho=1000.0)

#Plot the model
plotmodel(rho, nodes_crop, cropx, cropz, title="Starting Density", unit="density", cmap='gist_stern_r',
         vmin=1000,vmax=3400)

#########################
### ATTENUATION MODEL ###
#########################
q = qgeneration(vp_s, croppedbathy, spacing_s, max_q=10000.0, min_q=50.0)

#Plot the model
plotmodel(q, nodes_crop, cropx, cropz, title="Starting Attenuation", unit="attenuation", cmap='coolwarm')

"""~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BREAK ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""

##########################
### WRITE BINARY FILES ###
##########################

#Output file directory
direct = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/008_2020RealDataTests/043_MiscTestsFGProg1_SmthCnstL1/_finalmodels/" 

#Velocity - true
#a = np.array(vp_s.T,"float32")
#output_file = open(direct+'vp', 'wb')
#a.tofile(output_file)
#output_file.close()

#Velocity - true
#a = np.array(vp_s.T,"float32")
#output_file = open(direct+'vp_smooth_detrended_preupsc', 'wb')
#a.tofile(output_file)
#output_file.close()

#Velocity - starting
#a = np.array(vp_s.T,"float32")
#output_file = open(direct+'vp_smooth', 'wb')
#a.tofile(output_file)
#output_file.close()

#Density
#a = np.array(rho.T,"float32")
#output_file = open(direct+'rho', 'wb')
#a.tofile(output_file)
#output_file.close()

#Attenuation
#a = np.array(q.T,"float32")
#output_file = open(direct+'qp', 'wb')
#a.tofile(output_file)
#output_file.close()

#Bathymetry
#a = np.array(croppedbathy,"float32")
#output_file = open(direct+'fbathy_FirstLastOBSCrop', 'wb')
#a.tofile(output_file)
#output_file.close()