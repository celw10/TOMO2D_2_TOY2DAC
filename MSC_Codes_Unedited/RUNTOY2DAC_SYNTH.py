#! /usr/bin/env python

#Run me in the terminal

################################################################
# Version 1.3 mod 
# Christopher Williams, Memorial University of Newfoundland
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Whats New in version 1.3?
# -> Version tailored to synthetic inversions
# -> Changes made to match with the latest real-data RUNTOY2DAC
#    (version 2.2)
# -> Properly impliment parallelism
# -> Update README.out and INVNAMES.out
# -> Optionally output the gradient 
# -> Added additional multiscale options
################################################################

import os
import time
tt = time.time()

####################
### DEFINE FLAGS ###
####################

#Define the number of cores to be run in parallel, must be str
num_cores = "8"

#Output multiscale parameters to README.out file, 0=off, 1=om
readme_out = 1

#Output inversion names to to INVNAMES.out file, 0=off, 1=om
inv_names_out = 1

#Output gradient models, 0=no, 1=yes
output_grad = 1

##########################
### DEFINE DIRECTORIES ###
##########################

#Toy2dac executable directory
toy2dac_dir = "/home/celw10/2_TOY2DAC/bin/toy2dac"

#Where everything will be placed
relative_path = "00RESULTS/003_MultiSparse_Smth_Tau_S1/"

#Make a directory for all final velocity models
vel_models = relative_path+"0_inverted_models/"
os.mkdir(vel_models)

#Make a directory for all final optimization outputs
opt_out = relative_path+"1_optimization_optput/"
os.mkdir(opt_out)

#Make a directory for all gradient output
grad_out = relative_path+"2_gradient_optput/"
os.mkdir(grad_out)

####################
### DEFINE FILES ###
####################
    
#Define optimization output
opt_mtd_out = "iterate_PLB.dat"

#Define starting model
start_mdl = "vp_smooth"

#Velocity Models 
mdl = ["vp_Marmousi_exact " "qp " "rho ", "param_vp_final " "qp " "rho "]

##############################
### DEFINE INVERSION LOOPS ###
##############################

#Frequency loop - input as groups
freq = [[2.0,2.25,2.5],[2.0,2.5,3.0],[2.0,2.5,3.0,3.5],[2.0,3.0,3.5,4.0],[2.0,3.0,3.5,4.0,4.5],[2.0,3.0,4.0,4.5,5.0],[2.0,3.0,4.0,4.5,5.0,5.5],[2.0,3.0,4.0,5.0,5.5,6.0],[2.0,3.0,4.0,5.0,5.5,6.0,6.5],[2.0,3.0,4.0,5.0,6.0,6.5,7.0],[2.0,3.0,4.0,5.0,6.0,6.5,7.0,7.5],[2.0,3.0,4.0,5.0,6.0,7.0,7.5,8.0],[2.0,3.0,4.0,5.0,6.0,7.0,7.5,8.0,8.5],[2.0,3.0,4.0,5.0,6.0,7.0,8.0,8.5,9.0],[2.0,3.0,4.0,5.0,6.0,7.0,8.0,8.5,9.0,9.5],[2.0,3.0,4.0,5.0,6.0,7.0,7.5,8.0,9.0,9.5,10.0],[2.0,3.0,4.0,5.0,6.0,7.0,7.5,8.0,9.0,9.5,10.0,10.5],[2.0,3.0,4.0,5.0,6.0,7.0,7.5,8.0,9.0,10.0,10.5,11.0],[2.0,3.0,4.0,5.0,6.0,7.0,7.5,8.0,9.0,10.0,10.5,11.0,11.5],[2.0,3.0,4.0,5.0,6.0,7.0,7.5,8.0,9.0,10.0,11.0,11.5,12.0]]

#Time-damping loop
tau = [0.0] #[1.0, 0.5, 0.1, 0.0]

#Offset loop
offset = ["acqui_4km", "acqui_8km", "acqui_12km", "acqui_16km"]

#Gradient Preconditioner Loop
epsillon = [1e-2,1e-3,1e-4,1e-5]

#Iterations loop
iterations = [[5,5,5,5]] #[[5],[5,5],[5,5,5],[5,5,5,5]]

#Lambda (Tikhonov Regularization) loop
#Negative lambda activates gradient smoothing
lamb = [0.0,0.0,0.0,0.0]
#Lambda x, and z loop
lamb_x = [0.8, 0.6, 0.4, 0.2]
lamb_z = [0.4, 0.3, 0.2, 0.1]

#Deadzone loop
deadzone = [0.0,0.0,0.0,0.0]

#Data Weighting loop
data_weight = ["2 20000", "2 20000", "2 20000", "2  20000"]

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
    
#Copy smooth model to param_vp_final for first iteration
os.system("cp "+start_mdl+" param_vp_final")

#List of all output names
all_output_names = []

#FREQUENCIES
for f in range(len(freq)):
    #Set the number of frequencies in a group
    replace_line('freq_management', 0, str(len(freq[f]))+' ! Nfreq'+'\n')
    #Set the frequencies
    replace_line('freq_management', 1,  ", ".join(map(str, freq[f]))+'\n')
    
    #EXPONENTIAL TIME DAMPING
    for t in range(len(tau)):
        #Define the constant tau
        replace_line('fdfd_input', 7, str(tau[t])+' ! laplace constant'+'\n')
        
        #LOOP THROUGH OFFSET, ITERATIONS, GRADIENT PRECONDITIONING
        for o in range(len(iterations[t])):
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
            replace_line('fwi_input', 5, str(deadzone[o])+' ! deadzone'+'\n')
            #Define data weight
            replace_line('fwi_input', 22, str(data_weight[o])+' ! number of values in weighting file and spatial step value'+'\n')
            
            #LOOP THROUGH FORWARD MODELING AND INVERSION
            for m in range(2):
                #Define models
                replace_line('fdfd_input', 2, str(mdl[m])+' ! files name'+'\n')
                #Define modeling or inversion
                replace_line('toy2dac_input', 0, str(m)+' ! mode of the code (0 = MODELING, 1 = INVERSION)'+'\n')

                #RUN TOY2DAC
                os.system("mpirun "+"-np "+num_cores+" "+toy2dac_dir)

            ###################
            ### SAVE OUTPUT ###
            ###################
            
            #One label per loop
            label="FG"+str(f)+"_TAU"+str(tau[t])+"_Off"+str(offset[o])
      
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
print "|	        SUCCESSFULLY COMPLETED ALL JOBS            |"
print "|	    	        FINAL RUN TIME:                    |"
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "|      JOBS COMPLETED AT: {0}s".format(time.time()-tt),      "|"
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print 
