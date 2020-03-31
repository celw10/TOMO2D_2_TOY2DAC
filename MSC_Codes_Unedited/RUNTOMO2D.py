#! /usr/bin/env python

#Run me in the terminal

################################################################
# Version 1.1  
# Christopher Williams, Memorial University of Newfoundland
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Whats the purpose of RUNTOY2DAC?
# -> Gives the user the ability to automatically run a multi-
#    scale tomographic inversion creating organized, consistent
#    output.
# -> NOTE: This is only for refractions. We are trying to make 
#    sure our modeled first breaks are half a wavelength from 
#    the observed first breaks
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# What must be defined?
# -> In version 1.1, directories, files, and parameters must be 
#    defined.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Whats New in version 1.1 ?
# -> The code will now break if the inversion encounters an 
#    error.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Whats upcoming?
# -> The option to output an image of the output model after 
#    each inversion group
# -> The option to output a cycle-skip risk assesment after each
#    inversion group/or final model.
################################################################

import os
import time
tt = time.time()

####################
### DEFINE FLAGS ###
####################

#Coming soon...

##########################
### DEFINE DIRECTORIES ###
##########################

#Path to main directory for output
relative_path = "40_pick8_37run_July16/"

#Path to geometry files
geompath_inv = "picks_inv/"
geompath_fwd = "picks_fwd/"

#Inversion geometries
invgeoms = ["pg.0_02.all", "pg.0_04.all", "pg.0_06.all", "pg.0_08.all", "pg.0_10.all", "pg.0_12.all", "pg.0_14.all", "pg.0_16.all", "pg.0_18.all", "pg.0_20.all"]

#Forward modeling geometries
fwdgeoms = ["pg.geom.0_02.all", "pg.geom.0_04.all", "pg.geom.0_06.all", "pg.geom.0_08.all", "pg.geom.0_10.all", "pg.geom.0_12.all", "pg.geom.0_14.all", "pg.geom.0_16.all", "pg.geom.0_18.all", "pg.geom.0_20.all"]

#Make a directory for inversion output
vels_out = relative_path+"0_inverted_models/"
os.mkdir(vels_out)

#Make a directory for logfile output
logs_out = relative_path+"1_log_files/"
os.mkdir(logs_out)

#Make a directory for raypaths output
rays_out = relative_path+"2_raypaths/"
os.mkdir(rays_out)

#Make a directory for modeled traveltimes output
fwdtts_out = relative_path+"3_modeled_traveltimes/"
os.mkdir(fwdtts_out)

####################
### DEFINE FILES ###
####################
    
#Define starting model
start_mdl = "horavg_kwsmesh.dat"

#Define velocity correlation file
vcorr = "vcorrcelw.dat"

#########################
### DEFINE PARAMETERS ###
#########################

#Iterations
it = "5"

#Depth kernel weighting factor - <1 makes first arrivials priority
weight = "1"

#Velocity smoothing - > produces smoother
wsv = "50"

#Depth smoothing - > produces smoother
wsd = "50"

#Max velocity perturbation (%)
max_dv = "10" 

#Tolerance for LSQR algorithm 
lsqr_tol = "1e-4"

#Set inversion/fwd modeling parameters* 
#Forward star order for the graph method
#Maximum Segment length (clen), and number of interpolation points (nintp)
#Tolerance levels for iterations, [conjugate gradient (GC), Brent minimization (BM)] using the bending method
#"fwd*x/fwd*z/clen/nintp/tolCG/tolBM"
inv_fwd_params = "10/10/1/6/1e-5/1e-5"

#######################
### BEGIN INVERSION ###
#######################

tt = time.time()

#Copy smooth model to param_vp_final for first iteration
os.system("cp "+start_mdl+" inv_out.smesh.0.1")

#Number of inversion "groups" based on maximum offset
num_lines = len(invgeoms)
assert len(invgeoms) == len(fwdgeoms)

#Invert for each maximum offset group
for o in range(num_lines):

    #Inversion geometry file for current maximum offset group
    current_igeom = relative_path + geompath_inv + invgeoms[o]
 
    #Run Inversion
    os.system("tt_inverse -Minv_out.smesh."+str(o)+".1 -G"+current_igeom+" -N"+inv_fwd_params+" -Llog_group_"+str(o)+" -W"+weight+" -l -V -Q"+lsqr_tol+" -I"+it+" -SV"+wsv+" -SD"+wsd+" -CVvcorrcelw.dat -TV"+max_dv)
    
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "|		      INVERSION COMPLETE		   |"
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    
    #Forward Modeling geometry file for current maximum offset group
    current_fgeom = relative_path + geompath_fwd + fwdgeoms[o]

    #Forward model the dataset
    os.system("tt_forward -Minv_out.smesh."+str(o)+".1 -G"+current_fgeom+" -V -N"+inv_fwd_params+" -Rrays_group_"+str(o)+" > fwdtts_group_"+str(o))
    
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "|		   FORWARD MODELING COMPLETE		   |"
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    
    #Add a number to keep files in order (label)
    if o < 10:
        lbl = "0"+str(o)+"_"
    elif o < 100:
        lbl = +str(o)+"_"
    
    #Copy the velocity model 
    os.system("cp inv_out.smesh."+str(o)+".1 "+vels_out+lbl+"smesh")
    os.system("cp inv_out.smesh."+str(o)+".1 inv_out.smesh."+str(o+1)+".1")
    
    #Copy the log file
    os.system("mv log_group_"+str(o)+" "+logs_out+lbl+"log")
    
    #Copy the raypaths file
    os.system("mv rays_group_"+str(o)+" "+rays_out+lbl+"rays")
    
    #Copy the forward modeled traveltimes
    os.system("mv fwdtts_group_"+str(o)+" "+fwdtts_out+lbl+"fwdtts")

    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "|		COMPLETED ONE OFFSET GROUP, "+str(o),"		|"
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#Clean up
os.system("rm inv_out.smesh."+str(o)+".1")
os.system("rm inv_out.smesh."+str(o+1)+".1")

print
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "|	        SUCCESSFULLY COMPLETED ALL JOBS        |"
print "|	    	        FINAL RUN TIME:                |"
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "|        JOBS COMPLETED AT: {0}s".format(time.time()-tt),"	    |"
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print 