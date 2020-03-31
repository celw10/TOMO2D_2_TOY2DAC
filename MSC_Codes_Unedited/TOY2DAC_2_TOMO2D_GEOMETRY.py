#! /usr/bin/env python

#Run me in the terminal or jupyter notebook 

"""
TOY2DAC_to_TOMO2D_GEOMETRY - Version 1.0
Subroutine to convert from TOY2DAC geometry to TOMO2D Geometry
Version 1.0: Chris Williams, March 2019
    Original ipython notebook implimentation. 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Notes:
    (1) TOY2DAC acquisition files may only be converted into TOMO2D geometry files. There no first break information
contained in a TOY2DAC file. 
    (2) TOY2DAC uses acquisition files. TOMO2D geometry files are based on first break refraction/reflection picks.
We are going to pretend that there is a refracted arrival observed at each location specified in the TOY2DAC file.
Therefore all the modes are flagged as refractions in the converted TOMO2D geometry file.
    (3) This process also produces two geometry files required for the cycle skipping analysis. A source geometry
file and a tx.in file. 
"""

import numpy as np
import os

###################
### IMPORT FILE ###
###################

#Toy2dac geometry file
dirr = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/000_DATA/2_RealData/17_110StartMdl_THEEND/01_acqui_t0/"
f = "acqui_100"
f = dirr + f
acqui_toy = np.genfromtxt(f)

#Define output director
dirr = "/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/00_TOMO2DCONVERTED/2_EmedSynthetic/0_FinalSourceGeoms/"

#Define output geometry file name
geomnameout = "delete_me"

#Define output source geometry file name
sourcegeomout = "Synthetic_SourceGeom_Final_True_ReflectionCrop"

################
### OPTIONAL ###
################

#Shift sources if source and receiver lie on top of eachother (meters) (Runs into an error if so)
shiftx = -12922
shiftz = 0

##########################
### SCAN GEOMETRY FILE ###
##########################

#Count sources, doccument indicies
obsidx = []
count = 0 
for i in range(len(acqui_toy[:,0])):
    if acqui_toy[i,4] == 0:
        count += 1
        obsidx.append(i)
obsidx.append(len(acqui_toy[:,0]))
        
print "Number of sources: ", count

#Extract number of traces per source
obstraces = []
for i in range(len(obsidx)-1):
    obstraces.append(obsidx[i+1]-obsidx[i]-1)

print "Minimum number of traces per source: ", min(obstraces)
print "Maximum number of traces per source: ", max(obstraces)
print "Mean number of traces per source: ", float(sum(obstraces))/float(count)
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

#################################
### MAKE TOMO2D Geometry FILE ###
#################################

#Make blank TOMO2D geometry matrix
geom_tomo = np.zeros((len(acqui_toy[:,0])+1,6)).astype(object)
srccount = 0

#Loop through all points
for i in range(len(acqui_toy[:,0])+1):
    #Top line defines all shots
    if i == 0: 
        geom_tomo[i,0],geom_tomo[i,1],geom_tomo[i,2],geom_tomo[i,3],geom_tomo[i,4],geom_tomo[i,5] = \
        count, float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN')
    #Define location of source (optional shift for Marmousi acquisition)
    elif i-1 in obsidx:
        geom_tomo[i,0],geom_tomo[i,1],geom_tomo[i,2],geom_tomo[i,3],geom_tomo[i,4],geom_tomo[i,5] = \
        "s", (acqui_toy[i-1,1]+shiftx)/1000.0, (acqui_toy[i-1,0]+shiftz)/1000.0, float('NaN'), \
        obstraces[srccount], float('Nan')
        srccount += 1
    #Define subsiquent trace locations
    else:
        geom_tomo[i,0],geom_tomo[i,1],geom_tomo[i,2],geom_tomo[i,3],geom_tomo[i,4],geom_tomo[i,5] = \
        "r", acqui_toy[i-1,1]/1000.0, acqui_toy[i-1,0]/1000.0, 0, 0, 0
        
#################################
### MAKE TX.IN Geometry FILE ###
#################################

#################################
### MAKE SOURCE GEOMETRY FILE ###
#################################

#Blank source geometry file, stored as [x,z]
source_geom = np.zeros((srccount,2))
j=0

#Extract x and z position of each source
for i in range(len(geom_tomo[:,0])):
    if geom_tomo[i,0] == "s":
        source_geom[j,0], source_geom[j,1] = geom_tomo[i,1], geom_tomo[i,2]
        j+=1
        
#####################################
### SAVE THE TOMO2D GEOMETRY FILE ###
#####################################

#Setup text file with TOMO2D geometry
#np.set_printoptions(threshold=np.inf)
modeltmp = dirr+"tmp"
modelout = dirr+geomnameout
np.savetxt(modeltmp, geom_tomo, fmt='%5s', delimiter='\t')
np.savetxt(dirr+sourcegeomout, source_geom, fmt='%5s', delimiter='\t')

#Remove all "nan" from file
with open(modeltmp) as infile, open(modelout, "w") as outfile:
    for line in infile:
        outfile.write(line.replace("nan", ""))
        
#Delete temporary file
os.system("rm "+dirr+"tmp")