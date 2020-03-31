#RUN ME IN JUPYTER NOTEBOOK, SEPARATE OUT THE FUNCTIONS FROM THE MAIN CODE AT THE END

def InversionStatistics(iterations, normmisfit, absolutemisfit, ngrad, path, min_max,
                        figsize=[26,13], convcriterion=0.1, fontscalar=2.0, leg_shift=[0,0],
                        xtickidx=None, label=None, line=None, linelbl=None):

    def make_patch_spines_invisible(ax):
        ax.set_frame_on(True)
        ax.patch.set_visible(False)
        for sp in ax.spines.values():
            sp.set_visible(False)

    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = figsize[0]
    fig_size[1] = figsize[1]
    plt.rcParams["figure.figsize"] = fig_size

    fig, host = plt.subplots()
    fig.subplots_adjust(right=0.75)

    par1 = host.twinx()
    par2 = host.twinx()

    # Offset the right spine of par2.  The ticks and label have already been
    # placed on the right by twinx above.
    par2.spines["right"].set_position(("axes", 1.2))
    # Having been created by twinx, par2 has its frame off, so the line of its
    # detached spine is invisible.  First, activate the frame but make the patch
    # and spines invisible.
    make_patch_spines_invisible(par2)
    # Second, show the right spine.
    par2.spines["right"].set_visible(True)

    #p1, = host.semilogy(np.arange(len(normmisfit)), absolutemisfit, "bo--", linewidth=2.5*fontscalar,
    #                markersize=5*fontscalar, label="Total Cost Function")
    p1, = host.plot(np.arange(len(normmisfit)), absolutemisfit, "bo--", linewidth=2.5*fontscalar,
                    markersize=5*fontscalar, label="Total Cost Function")
    p2, = par1.plot(np.arange(len(absolutemisfit)), normmisfit, "ro--", linewidth=2.5*fontscalar, 
                    markersize=4*fontscalar, label="Relative Cost Function")
    p3, = par2.plot(np.arange(len(ngrad)), ngrad, "g-.", linewidth=1*fontscalar,
                    label="Number of Computed Grads")

    if xtickidx == None:
        host.set_xlim(min(iterations), max(iterations))
    else:
        host.set_xlim(min(iterations), max(iterations))
        host.set_xticks(xtickidx)
        host.set_xticklabels(range(0,len(xtickidx),1))
    host.set_ylim(min_max[0][0], min_max[1][0])
    par1.set_ylim(min_max[0][1], min_max[1][1])
    par2.set_ylim(min_max[0][2], min_max[1][2])

    host.set_xlabel("Iteration", fontsize = 24*fontscalar)
    host.set_ylabel("Total Cost Function", fontsize = 24*fontscalar)
    par1.set_ylabel("Relative Cost Function", fontsize = 24*fontscalar)
    par2.set_ylabel("Number of Computed Gradients", fontsize = 24*fontscalar)
    host.tick_params(axis='x', labelsize=20*fontscalar)
    host.tick_params(axis='y', labelsize=20*fontscalar)
    par1.tick_params(axis='y', labelsize=20*fontscalar)
    par2.tick_params(axis='y', labelsize=20*fontscalar)

    host.yaxis.label.set_color(p1.get_color())
    par1.yaxis.label.set_color(p2.get_color())
    par2.yaxis.label.set_color(p3.get_color())

    tkw = dict(size=12*fontscalar, width=2)
    host.tick_params(which='both', axis='y', colors=p1.get_color(), **tkw)
    par1.tick_params(which='both', axis='y', colors=p2.get_color(), **tkw)
    par2.tick_params(which='both', axis='y', colors=p3.get_color(), **tkw)
    host.tick_params(which='both', axis='x', **tkw)

    lines = [p1, p2, p3]

    host.legend(lines, [l.get_label() for l in lines], loc='center right', fontsize = 20*fontscalar)
    
    plt.title("Inversion Statistics", fontsize=28*fontscalar)
    
    #Add a number to keep figures in order
    if len(normmisfit)-1 < 10:
        idx = "000"+str(len(normmisfit)-1)+"_"
    elif len(normmisfit)-1 < 100:
        idx = "00"+str(len(normmisfit)-1)+"_"
    elif len(normmisfit)-1 < 1000:
        idx = "0"+str(len(normmisfit)-1)+"_"
    else:
        idx = str(len(normmisfit)-1)+"_"
        
    #Plot lines showing where select frequency groups begin
    if line != None:
        i=0
        while i < len(line):
            if line[i] < len(normmisfit):
                #Figure out the label
                if linelbl != None:
                    linelabel = linelbl[i]
                    host.axvline(x=line[i], color="black", ls="--")
                    host.text(line[i]+1, min_max[1][1]*0.275, linelabel, #min_max[1][0]#Second term important, how far up the line to position
                              rotation=90, fontsize=16*fontscalar)
                #else:
                    #host.axvline(x=line[i], color="black", ls="--")
                i+=1
            else: 
                i=len(line)
    
    #Save output
    #os.mkdir(path+"4_invstats/")
    if label != None:
        plt.gcf().savefig(path+"7_invstats/"+str(idx)+str(label), bbox_inches = 'tight', pad_inches = 0.1)
    else:
        plt.gcf().savefig(path+"7_invstats/"+str(idx)+"INVSTATS_OUT.png", bbox_inches = 'tight', pad_inches = 0.1)
    
    plt.tight_layout()
    
    #Delete current figure
    fig.clear(); del fig
    plt.close()

    #plt.show()
    
    #print 'Plot completion time: {0}s'.format(time.time()-tt)
    #print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    
import numpy as np
from matplotlib import pyplot as plt, cm, colors
from mpl_toolkits.mplot3d import Axes3D
import time
import os
tt = time.time()

####################
### DEFINE FLAGS ###
####################

#0 = just output the final figure, 1 = output an image for each iteration (to produce gifs with MakeMovie)
figs_out = 0

##########################
### DEFINE DIRECTORIES ###
##########################

#Where to look for output (unaltered from RUNTOY2DAC_V***)
relative_path = "/home/chris/4_TOY2DAC_V2.5_2018_07_13/EM_L1/3_EML1_Torngat/009_Thesis_Inversions/2_EMReal/0_Manu_Results_Final_043/"

##################################
### PLOT OF INVERSION STASTICS ###
##################################

# README: Columns for output file
# ~~~ Column: (0) Nonlinear iteration number, (1) non normalized cost function value, (2) norm of the gradient, 
# ~~~ (3) the relative value of the cost function, normalized by the initial value, (4) size of the steplength, 
# ~~~ (5) number of linesearch iterations for determining steplength, (6) the total number of gradient computation

#Read INVNAMES.out
f = open(relative_path+"INVNAMES.out","r")
model_names = f.read().splitlines()
f.close()

#Read in inversion output
invstats = []
filelengths = []
flag  = ""
for i in range(len(model_names)):
    f = relative_path+"1_optimization_optput/"+str(model_names[i])
    num = -2
    tmp = open(f, 'r')
    for line in tmp:
        line_data = line.split()
        if len(line_data) == 7:
            num += 1
            #We only want the very first zero iteration, not subsequent zero iterations
            if num > -1 and flag=="":
                flag=0
                niter, fk, gk, fk_f0, alpha, nls, ngrad = line.split()
                invstats.append([int(niter), float(fk), float(gk), float(fk_f0), float(alpha), int(nls), int(ngrad)])
                init_cost = float(fk)
            elif num > 0 and flag!="":
                niter, fk, gk, fk_f0, alpha, nls, ngrad = line.split()
                invstats.append([int(niter), float(fk), float(gk), float(fk_f0), float(alpha), int(nls), int(ngrad)])
    #Account for the very first zero iteration
    if flag==0:
        filelengths.append(num+1)
        flag=1
    else:
        filelengths.append(num+filelengths[i-1])
        
        ###########################################################
        #### CONTINUE MY ATTEMPT AT THIS WHEN I HAVE MORE TIME ####
        # MAJOR ISSUE IS ALLIGNMENT, I NEED TO START OFF AFTER THE FIRST GROUP AS WE DON'T START 
        # THE OTHER GIFS AT THE STARTING MODEL
        # I ALSO NEED TO ACCOUNT FOR SKIPPED ITERATIONS DUE TO AN UNSATISIFIED MODEL UPDATE CRITERIA
        
#Save as matrix
InvStats = np.zeros((len(invstats), 7))
for i in range(len(invstats)):
    InvStats[i,:] = invstats[i]
    
#Join inversion stats files over separate inversion groups
for i in range(len(filelengths)-1):
    #Cumulative iterations over inversion groups
    InvStats[filelengths[i]:filelengths[i+1],0] = InvStats[filelengths[i]:filelengths[i+1],0] + InvStats[filelengths[i]-1,0]
    #Cumulative number of gradients
    InvStats[filelengths[i]:filelengths[i+1],6] = InvStats[filelengths[i]:filelengths[i+1],6] + InvStats[filelengths[i]-1,6]
    
#Figure out the minimum and maximum values for each plot
min_max = [[],[]]
#Total cost function
min_max[0].append(min(InvStats[:,1]))
min_max[1].append(max(InvStats[:,1]))
#Normalized cost function
min_max[0].append(0.95)
min_max[1].append(1.001) #This wont be higher than 1.0
#Total gradient computations
min_max[0].append(min(InvStats[:,6]))
min_max[1].append(max(InvStats[:,6]))

#In case you want to plot up some lines showing frequency groups
#xpoints=[]
#for i in range(11,len(filelengths),12):
#    xpoints.append(filelengths[i])
xpoints=[]
for i in range(0,len(filelengths)-1,3):
    xpoints.append(filelengths[i]-1)

if figs_out == 0:
    #Producce a single plot
    InversionStatistics(InvStats[:,0], InvStats[:,3], InvStats[:,1], InvStats[:,6], relative_path, min_max,
                   fontscalar=1.0, xtickidx=None)
    
elif figs_out == 1:
    idx=0
    #Produce a plot for each iteration to animate the inversions progress
    for i in range(1,len(InvStats[:,0])+1,1):
        if i > filelengths[idx]:
            idx+=1
            label=model_names[idx]
        else:
            label=model_names[idx]
        InversionStatistics(InvStats[:,0], InvStats[0:i,3], InvStats[0:i,1], InvStats[0:i,6], relative_path, min_max,
                           fontscalar=1.5, leg_shift=[1,1], xtickidx=None, label=label, line=xpoints,
                           linelbl=None)