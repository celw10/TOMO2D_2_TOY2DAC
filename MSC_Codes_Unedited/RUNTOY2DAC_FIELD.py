#! /usr/bin/env python

#Run me in the terminal

################################################################
# Version 2.0
# Christopher Williams, Memorial University of Newfoundland
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Whats New in version 2.0?
# -> RUNTOY2DAC Built up from version 1.2.2 mod to run real data
#    inversions.
# -> There is no loger a forward modeling/inversion loop
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Whats new in version 2.1?
# -> I've added an output to save the gradient for each 
#    inversion group.
# Whats new in version 2.1.1?
# -> Added an option to have a multiscale deadzone above which  
#    the odel will not be updated (other than regularization)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Whats New in version 2.2?
# -> Added the optionality of a multiscale data weighting 
#    function.
# -> Altered how the laplace constant is applied to suite the
#    1% 10% 100% approach (i.e. precent of data passed at given
#    time).
# -> Changed the README.out file to incorporate more of the 
#    user defined settings.
# -> Changed the INVNAMES.out file to update after every
#    inversion group.
# -> Added the proper parallel implimentation...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Whats New in version 2.3?
# -> Automatically run a forward model following the inversion
################################################################

import os
import time
tt = time.time()
tt_inv = time.time()

####################
### DEFINE FLAGS ###
####################

#Define the number of cores to be run in parallel, must be str
num_cores = "16"

#Output multiscale parameters to README.out file, 0=off, 1=om
readme_out = 1

#Output inversion names to to INVNAMES.out file, 0=off, 1=om
inv_names_out = 1

#Output gradient models, 0=no, 1=yes
output_grad = 0

#Forward model after the inversion, 0=no, 1=yes
fwd_mdl_flag = 1

##########################
### DEFINE DIRECTORIES ###
##########################

#Toy2dac executable directory
toy2dac_dir = "/home/celw10/2_TOY2DAC/bin/toy2dac"

#Where everything will be placed
relative_path = "00RESULTS/006_Nov_26_19/"

#Make a directory for all final velocity models (Keep Consistent)
vel_models = relative_path+"0_inverted_models/"
os.mkdir(vel_models)

#Make a directory for all final optimization outputs (Keep Consistent)
opt_out = relative_path+"1_optimization_optput/"
os.mkdir(opt_out)

#Make a directory for all final optimization outputs (Keep Consistent)
grad_out = relative_path+"2_gradient_optput/"
os.mkdir(grad_out)

####################
### DEFINE FILES ###
####################
    
#Define optimization output
opt_mtd_out = "iterate_PLB.dat"

#Define starting model
start_mdl = "vp_smooth"

#Define the attenuation and rho file, starting model copied to param_vp_final to start inversion
mdl = ["param_vp_final ", " qp ", "rho "]

#Define the attenuation and rho file for forward modeling
mdl_fw = ["param_vp_final ", " qp ", "rho "]

#############################
### DEFINE INV/FWD PARAMS ###
#############################

#Inversion nodes [x,z]
nodes_inv = [3801, 601] 

#Inversion model spacing
space_inv = 50.0

#Define the laplace constant for forward modeling
fwd_LC = 0.0

#Forward modeling offset
fwd_OFF = "acqui_100"

#Forward modeling nodes [x,z]
nodes_fwd = [3801, 601] 

#Forward modeling spcaing
space_fwd = 50.0

##############################
### DEFINE INVERSION LOOPS ###
##############################

#Frequency loop - input as groups
freq=[[2.0,2.25,2.5],[2.0,2.5,3.0],[2.0,2.5,3.0,3.5],[2.0,3.0,3.5,4.0],[2.0,3.0,3.5,4.0,4.5],[2.0,3.0,4.0,4.5,5.0],[2.0,3.0,4.0,4.5,5.0,5.5],[2.0,3.0,4.0,5.0,5.5,6.0],[2.0,3.0,4.0,5.0,5.5,6.0,6.5],[2.0,3.0,4.0,5.0,6.0,6.5,7.0]]

#Time-damping loop - define by % then offset
tau=[[5.0,5.0,5.0,5.0,5.0],[2.0,2.0,2.0,2.0,2.0],[1.0,1.0,1.0,1.0,1.0],[0.5,0.5,0.5,0.5,0.5]]

#Offset loop
offset=["acqui_015","acqui_030","acqui_045","acqui_060","acqui_100"]

#Dataset loop
dataset=[[["DM_FG01_LC1_015","DM_FG01_LC1_030","DM_FG01_LC1_045","DM_FG01_LC1_060","DM_FG01_LC1_100"],["DM_FG01_LC2_015","DM_FG01_LC2_030","DM_FG01_LC2_045","DM_FG01_LC2_060","DM_FG01_LC2_100"],["DM_FG01_LC3_015","DM_FG01_LC3_030","DM_FG01_LC3_045","DM_FG01_LC3_060","DM_FG01_LC3_100"],["DM_FG01_LC4_015","DM_FG01_LC4_030","DM_FG01_LC4_045","DM_FG01_LC4_060","DM_FG01_LC4_100"]],[["DM_FG02_LC1_015","DM_FG02_LC1_030","DM_FG02_LC1_045","DM_FG02_LC1_060","DM_FG02_LC1_100"],["DM_FG02_LC2_015","DM_FG02_LC2_030","DM_FG02_LC2_045","DM_FG02_LC2_060","DM_FG02_LC2_100"],["DM_FG02_LC3_015","DM_FG02_LC3_030","DM_FG02_LC3_045","DM_FG02_LC3_060","DM_FG02_LC3_100"],["DM_FG02_LC4_015","DM_FG02_LC4_030","DM_FG02_LC4_045","DM_FG02_LC4_060","DM_FG02_LC4_100"]],[["DM_FG03_LC1_015","DM_FG03_LC1_030","DM_FG03_LC1_045","DM_FG03_LC1_060","DM_FG03_LC1_100"],["DM_FG03_LC2_015","DM_FG03_LC2_030","DM_FG03_LC2_045","DM_FG03_LC2_060","DM_FG03_LC2_100"],["DM_FG03_LC3_015","DM_FG03_LC3_030","DM_FG03_LC3_045","DM_FG03_LC3_060","DM_FG03_LC3_100"],["DM_FG03_LC4_015","DM_FG03_LC4_030","DM_FG03_LC4_045","DM_FG03_LC4_060","DM_FG03_LC4_100"]],[["DM_FG04_LC1_015","DM_FG04_LC1_030","DM_FG04_LC1_045","DM_FG04_LC1_060","DM_FG04_LC1_100"],["DM_FG04_LC2_015","DM_FG04_LC2_030","DM_FG04_LC2_045","DM_FG04_LC2_060","DM_FG04_LC2_100"],["DM_FG04_LC3_015","DM_FG04_LC3_030","DM_FG04_LC3_045","DM_FG04_LC3_060","DM_FG04_LC3_100"],["DM_FG04_LC4_015","DM_FG04_LC4_030","DM_FG04_LC4_045","DM_FG04_LC4_060","DM_FG04_LC4_100"]],[["DM_FG05_LC1_015","DM_FG05_LC1_030","DM_FG05_LC1_045","DM_FG05_LC1_060","DM_FG05_LC1_100"],["DM_FG05_LC2_015","DM_FG05_LC2_030","DM_FG05_LC2_045","DM_FG05_LC2_060","DM_FG05_LC2_100"],["DM_FG05_LC3_015","DM_FG05_LC3_030","DM_FG05_LC3_045","DM_FG05_LC3_060","DM_FG05_LC3_100"],["DM_FG05_LC4_015","DM_FG05_LC4_030","DM_FG05_LC4_045","DM_FG05_LC4_060","DM_FG05_LC4_100"]],[["DM_FG06_LC1_015","DM_FG06_LC1_030","DM_FG06_LC1_045","DM_FG06_LC1_060","DM_FG06_LC1_100"],["DM_FG06_LC2_015","DM_FG06_LC2_030","DM_FG06_LC2_045","DM_FG06_LC2_060","DM_FG06_LC2_100"],["DM_FG06_LC3_015","DM_FG06_LC3_030","DM_FG06_LC3_045","DM_FG06_LC3_060","DM_FG06_LC3_100"],["DM_FG06_LC4_015","DM_FG06_LC4_030","DM_FG06_LC4_045","DM_FG06_LC4_060","DM_FG06_LC4_100"]],[["DM_FG07_LC1_015","DM_FG07_LC1_030","DM_FG07_LC1_045","DM_FG07_LC1_060","DM_FG07_LC1_100"],["DM_FG07_LC2_015","DM_FG07_LC2_030","DM_FG07_LC2_045","DM_FG07_LC2_060","DM_FG07_LC2_100"],["DM_FG07_LC3_015","DM_FG07_LC3_030","DM_FG07_LC3_045","DM_FG07_LC3_060","DM_FG07_LC3_100"],["DM_FG07_LC4_015","DM_FG07_LC4_030","DM_FG07_LC4_045","DM_FG07_LC4_060","DM_FG07_LC4_100"]],[["DM_FG08_LC1_015","DM_FG08_LC1_030","DM_FG08_LC1_045","DM_FG08_LC1_060","DM_FG08_LC1_100"],["DM_FG08_LC2_015","DM_FG08_LC2_030","DM_FG08_LC2_045","DM_FG08_LC2_060","DM_FG08_LC2_100"],["DM_FG08_LC3_015","DM_FG08_LC3_030","DM_FG08_LC3_045","DM_FG08_LC3_060","DM_FG08_LC3_100"],["DM_FG08_LC4_015","DM_FG08_LC4_030","DM_FG08_LC4_045","DM_FG08_LC4_060","DM_FG08_LC4_100"]],[["DM_FG09_LC1_015","DM_FG09_LC1_030","DM_FG09_LC1_045","DM_FG09_LC1_060","DM_FG09_LC1_100"],["DM_FG09_LC2_015","DM_FG09_LC2_030","DM_FG09_LC2_045","DM_FG09_LC2_060","DM_FG09_LC2_100"],["DM_FG09_LC3_015","DM_FG09_LC3_030","DM_FG09_LC3_045","DM_FG09_LC3_060","DM_FG09_LC3_100"],["DM_FG09_LC4_015","DM_FG09_LC4_030","DM_FG09_LC4_045","DM_FG09_LC4_060","DM_FG09_LC4_100"]],[["DM_FG10_LC1_015","DM_FG10_LC1_030","DM_FG10_LC1_045","DM_FG10_LC1_060","DM_FG10_LC1_100"],["DM_FG10_LC2_015","DM_FG10_LC2_030","DM_FG10_LC2_045","DM_FG10_LC2_060","DM_FG10_LC2_100"],["DM_FG10_LC3_015","DM_FG10_LC3_030","DM_FG10_LC3_045","DM_FG10_LC3_060","DM_FG10_LC3_100"],["DM_FG10_LC4_015","DM_FG10_LC4_030","DM_FG10_LC4_045","DM_FG10_LC4_060","DM_FG10_LC4_100"]]]
 
#Gradient Preconditioner Loop
epsillon=[1e-2,1e-3,1e-4,1e-5,1e-6]

#Iterations loop
iterations=[[5],[5,5],[5,5,5,5],[5,5,5,5,5]]

#Lambda (Tikhonov Regularization) loop
#Negative lambda activates gradient smoothing
lamb = [-1.0,-1.0,-1.0,-1.0,-1.0]
#Lambda x, and z loop
lamb_x = [1.6,1.6,1.6,1.6,1.6]
lamb_z = [0.8,0.8,0.8,0.8,0.8]

#Data Weight File loop
data_weight_file = ["DW_0", "DW_1", "DW_2", "DW_3", "DW_4"]

#Data Weighting loop
data_weight = ["2 1000", "2 1000", "2 1000", "2  1000", "2 1000"]

##################################
### BEGIN MULTISCALE INVERSION ###
##################################

#Method of replacing lines in text files
def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()
    
#Define README.out file structure
def mult_readme_out(freq, tau, offset, epsillon, iterations, lamb, lamb_x, lamb_z):
    #Save multiscale inversion summary file
    f = open(relative_path+"README.out","w+")
    f.write("~~~~ MULTISCALE INVERSION PARAMETERS ~~~~ \n")
    f.write("Frequencies: " + str(freq)+"\n")
    f.write("Laplace-Fourier Constants: " + str(tau)+"\n")
    f.write("Maximum Offsets: " + str(offset)+"\n")
    f.write("Gradient Preconditioner Constants: " + str(epsillon)+"\n")
    f.write("Maximum Number of non-linear iterations: " + str(iterations)+"\n")
    f.write("Lambda for Tikhonov Regularization (negative activates gradient smoothing): " + str(lamb)+"\n")
    f.write("X-Weight for Tikhonov Regularization (x-weight of gradient smoothing as a function of wavelength): " + str(lamb_x)+"\n")
    f.write("Z-Weight for Tikhonov Regularization (z-weight of gradient smoothing as a function of wavelength): " + str(lamb_z)+"\n")
    f.write("Multiscale Deadzone values: " + str(deadzone)+"\n")
    f.write("Multiscale data weighting values (see data_weight_file for actual data weighting values): " + str(data_weight)+"\n")
    f.write("Total time to complete multiscale Inversion: {0}s".format(time.time()-tt)+"\n")
    f.close()

#Define INVNAMES.out file structure
def mult_inv_names_out(all_output_names):
    f = open(relative_path+"INVNAMES.out","w+")
    for i in range(len(all_output_names)):
        f.write(str(all_output_names[i])+"\n")
    f.close
    
#Model upscaling function for forward modeling final velocity model (Do this for density and qp prior to)
def upscale_toy2dac_mdl(smesh, nx, nz, upscale):
    
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "Mesh x-nodes: ", nx
    print "Mesh z-nodes: ", nz

    #Upscale model
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
    print "Mesh x-nodes have been increased to: ", nx
    print "Mesh z-nodes have been increased to: ", nz
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

    #Output scaled to meteres for TOY2DAC - doing so prior results in a memory error
    return smesh

#List of all output names
all_output_names = []

# TOY2DAC INPUT
#Set to inversion
replace_line('toy2dac_input', 0, '1 ! mode of the code (0 = MODELING, 1 = INVERSION)'+'\n')

# FREQUENCY MANAGEMENT
#Copy smooth model to param_vp_final for first iteration
os.system("cp "+start_mdl+" param_vp_final")
#Copy the inversion freq_management file to freq_management (SET TO freq_management_inv)
os.system("cp freq_management_inv freq_management")

# FDFD FILE
#Define the inversion model nodes
replace_line('fdfd_input', 0, str(nodes_inv[1])+" "+str(nodes_inv[0])+' ! nz,nx'+'\n')
#Define the ivnersion model spacing
replace_line('fdfd_input', 1, str(space_inv)+' ! h'+'\n')
#Define the inversion models
replace_line('fdfd_input', 2, str(mdl[0])+str(mdl[1])+str(mdl[2])+' ! files name'+'\n')

#FREQUENCIES
for f in range(len(freq)):
    #Set the number of frequencies in a group
    replace_line('freq_management', 0, str(len(freq[f]))+' ! Nfreq'+'\n')
    #Set the frequencies
    replace_line('freq_management', 1,  ", ".join(map(str, freq[f]))+'\n')
    
    #EXPONENTIAL TIME DAMPING
    for t in range(len(tau)):
        #Choose the damping value below, controlled by iterations.
        
        #LOOP THROUGH OFFSET, ITERATIONS, GRADIENT PRECONDITIONING
        for o in range(len(iterations[t])):
            #Define the Laplace Constant
            replace_line('fdfd_input', 7, str(tau[t][o])+' ! laplace constant'+'\n')
            #Define the offset
            replace_line('toy2dac_input', 2, str(offset[o])+' ! acquisition file name'+'\n')
            #Define the gradient preconditioner
            replace_line('fwi_input', 12, str(epsillon[o])+' ! threshold parameter for the preconditioner'+'\n')
            #Define the number of iterations
            replace_line('fwi_input', 9, str(iterations[t][o])+' ! number of nonlinear iterations'+'\n')
            #Define lambda
            replace_line('fwi_input', 14, str(lamb[o])+' ! lambda for Tikhonov regularization'+'\n')
            #Define lambda x and z
            replace_line('fwi_input', 15, str(lamb_x[o])+' '+str(lamb_z[o])+' ! lambda_x, lambda_z (directional weight for Tikhonov regularization)'+'\n')
            #Define deadzone
            replace_line('fwi_input', 21, str(data_weight_file[o])+' ! file for weighting'+'\n')
            #Define data weight
            replace_line('fwi_input', 22, str(data_weight[o])+' ! number of values in weighting file and spatial step value'+'\n')
            #Define dataset
            replace_line('fwi_input', 0, str(dataset[f][t][o])+'  !name of obs data'+'\n')
            
            #RUN TOY2DAC
            os.system("mpirun "+"-np "+num_cores+" "+toy2dac_dir)

            ###################
            ### SAVE OUTPUT ###
            ###################

            #Make a label for filenames
            label="FG"+str(f)+"_TAU"+str(t)+"_Off"+str(offset[o])
            
            #Copy inverted velocity model to specified directory
            os.system("cp param_vp_final "+vel_models+label)

            #Copy optimization routine output to specified directory
            os.system("cp "+opt_mtd_out+" "+opt_out+label)
            
            #Copy gradient output to specified directory
            if output_grad == 1:
                os.system("cp gradient "+grad_out+label)
            
            #Append name in loop
            all_output_names.append(label)

            #Output README.out
            if readme_out == 1:
                mult_readme_out(freq, tau, offset, epsillon, iterations, lamb, lamb_x, lamb_z)

            #Output INVNAMES.out
            if inv_names_out == 1:
                mult_inv_names_out(all_output_names)
                
print
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "| SUCCESSFULLY COMPLETED INVERSION"
print "| RUN TIME: {0}s".format(time.time()-tt_inv)
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print 

Inversion_Time = time.time()-tt_inv
                
###########################
### BEGIN FORWARD MODEL ###
###########################

if fwd_mdl_flag == 1:

    #Sleep for 30 seconds
    time.sleep(10)
    tt_fwd = time.time()

    # FREQUENCY MANAGEMENT
    #Remove the inversion freq_management file
    os.system("rm freq_management")
    #Copy the forward modeling freq_management file to freq_management (SET TO freq_management_fwd)
    os.system("cp freq_management_fwd freq_management")

    # FDFD FILE
    #Define the inversion model nodes
    replace_line('fdfd_input', 0, str(nodes_fwd[1])+" "+str(nodes_fwd[0])+' ! nz,nx'+'\n')
    #Define the ivnersion model spacing
    replace_line('fdfd_input', 1, str(space_fwd)+' ! h'+'\n')
    #Define the inversion models
    replace_line('fdfd_input', 2, str(mdl_fw[0])+str(mdl_fw[1])+str(mdl_fw[2])+' ! files name'+'\n')
    #Define the Laplace Constant (Zero for all forward modelings)
    replace_line('fdfd_input', 7, str(fwd_LC)+' ! laplace constant'+'\n')

    # TOY2DAC INPUT
    #Define the offset
    replace_line('toy2dac_input', 2, str(fwd_OFF)+' ! acquisition file name'+'\n')
    #Set to forward modeling
    replace_line('toy2dac_input', 0, '0 ! mode of the code (0 = MODELING, 1 = INVERSION)'+'\n')

    #RUN TOY2DAC
    os.system("mpirun "+"-np "+num_cores+" "+toy2dac_dir)

print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "| SUCCESSFULLY COMPLETED FORWARD MODEL"
print "| RUN TIME: {0}s".format(time.time()-tt_fwd)
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print 

print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "| JOB SUMMARY"
print "| RUN TIME - TOTAL: {0}s".format(time.time()-tt)
print "| RUN TIME - INVERSION: "+str(Inversion_Time)
print "| RUN TIME - FORWARD: {0}s".format(time.time()-tt_fwd)
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print 
