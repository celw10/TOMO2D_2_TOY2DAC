#! /usr/bin/env python

"""
Author: Christopher Williams
Last updated: April 3rd, 2019
Description: Convert a TOY2DAC velocity model to TOMO2D format. 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Notes:
    (1) TOMO2D operates in km and TOY2DAC operates in m. Will be working in km in this workbook.
"""

import numpy as np
import os

def read_binary(binary_file, nodes=None, reshape=False):
    
    """ 
    Author: Christopher Williams
    Last updated: April 3rd, 2019
    Description: read and reshape a binary file. 
    ~~~
    binary_file: Binary format model as param*final directly output from Toy2dac.
    nodes: Number of nodes in x- and z-dimension as [nx, ny]
    """
    
    with open(binary_file, 'rb') as f:
        tmp = np.fromfile(f, dtype=np.float32)
    
    if reshape == True:
        out = np.reshape(tmp, [nodes[0], nodes[1]])
        return out.T
    
    else:
        return tmp

#####################
### TOY2DAC MODEL ###
#####################

#Absolute path to model
fv = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/000_DATA/3_SynthEmed/4_SyntheticEMNov/model_50m/vp"

#Absolute path to bathymetry
fb = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/000_DATA/3_SynthEmed/4_SyntheticEMNov/model_50m/fbathy"

#Path for the output TOMO2D Model
path = "/home/chris/2_TOMO2D/TOMO2D_CELW/00CELW/00_TOMO2DCONVERTED/2_EmedSynthetic/"

#Name for the output model
modelnameout = "DUMMY"

#Define nodes [nx,nz]
nodes = [3801,601]

#Homogeneous model spacing in kilometers
spacing= 0.05

#Velocity of water and air [water,air] (km/s)
vels = [1.5,0.343]

#Sig Figs (remember, working in kilometers)
sigfigs = 7

#LINE 1
line1 = []
line1.append(nodes[0]), line1.append(nodes[1]), line1.append(vels[0]), line1.append(vels[1])

#LINE 2
line2 = []
for i in range(nodes[0]):
    line2.append(round(i*spacing,sigfigs))
    
#LINE 3
with open(fb, 'rb') as f:
    bathy = np.fromfile(f, dtype=np.float32)
line3=[]
for i in range(len(bathy)):
    line3.append(bathy[i]/1000.0)
    
#LINE 4
line4 = []
for i in range(nodes[1]):
    line4.append(round(i*spacing,sigfigs))

#LINE 5
with open(fv, 'rb') as f:
    data = np.fromfile(f, dtype=np.float32)
smesh = np.reshape(data, [nodes[0], nodes[1]])

#Setup text file with TOMO2D model
np.set_printoptions(threshold=np.inf)
modeltmp = path+"tmp"
modelout = path+modelnameout
with open(modeltmp, 'w+') as out:
    out.write(", ".join(map(str,line1[0:]))+"\n")
    out.write(", ".join(map(str,line2[0:]))+"\n")
    out.write(", ".join(map(str,line3[0:]))+"\n")
    out.write(", ".join(map(str,line4[0:]))+"\n")
    line5=[]
    for i in range(nodes[0]):
        for j in range(nodes[1]):
            #Don't include water velocities
            if smesh[i,j] > 1500:
                line5.append(round(smesh[i,j]/1000.0,sigfigs))
        out.write(", ".join(map(str,line5[1:-1]))+"\n")
        line5=[]

#Remove all commas from file
with open(modeltmp) as infile, open(modelout, "w") as outfile:
    for line in infile:
        outfile.write(line.replace(",", ""))
        
#Delete temporary file
os.system("rm "+path+"tmp")