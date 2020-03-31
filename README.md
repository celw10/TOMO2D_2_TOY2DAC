A collection of Python algorithms to format seismic data for traveltime tomography in TOMO2D and full-waveform inversion in TOY2DAC.
These codes are for an MSc. thesis titled "Wide-Angle, Full Waveform Inversion with a Sparse Ocean-Bottom Seismometer Dataset, Imaging the Cyprus Arc from the Eratosthenes Seamount to the Hecataeus Rise".
The current version of TOMO2D_2_TOY2DAC contains the unedited codes used in this project. 
Current efforts are focused on improving these codes and formatting them in a program.

All codes are provided as .py files and ran using Python 2.7. 
The only software required to execute some of these files not provided standard through the Anaconda installation of Python 2.7 is segyio, available here: https://github.com/equinor/segyio.git
Some functions are preferably ran in a Jupyter notebook. Others are preferably ran from the terminal. The preference is denoted at the top of each .py file. 
Brief descriptions for each module are provided as follows.

AMPLITUDE_PROCESS.py : Apply an amplitude bulk shift to the observed datsaet and automatically reduce the amplitude of noisy traces.

DATAMODELING_TDSEISMIC.py : Inverse Fourier transform a forward modelled seismic dataset from TOY2DAC view the data as a time-domain shot gather.

EDITACQUI.py : Edit a TOY2DAC acquisition file.
FWI_SG_ASSESS.py : Inverse Fourier transform a forward modelled seismic dataset from TOY2DAC. Optionally overlay the observed data on this dataset along with first break picks.

INVERSTION_STATS.py : Plot the cost function output from an appropiate RUNTOY2DAC executable.

MAKEFIGURES_RUNTOY2DAC.py : Make figures from the output of an appopiate RUNTOY2DAC executable.

MAKEFIGURES.py : A vast collection of plotting codes used to make different figures for the thesis/manuscript.

MATCH_REFLECTION_REFRACTION.py : Match the geometries of the seismic refraction dataset with a coincident seismic reflection line. Then crop the reflection line to the refraction line, and output as a GLOBE CLARITAS velocity file. This enables saving the data as segy where it may be imported into PETREL.


QC_FD_DATA_MODELLING.py : Plot frequency domain traces for the observed, predicted, and modelled datasets. 

RUNTOMO2D.py : Run a multiscale tomographic inversion.

RUNTOY2DAC_FIELD.py : Run multsicale FWI for a field datset.

RUNTOY2DAC_SYNTH.py : Run multiscale FWI for a synthetic datset.

TOMO2D_MDL_ASSMT_CycleSkip.py : Assess a starting model by analyzing the pick residuals for cycle skipping at a given frequency. 

TOMO2D_MDL_ASSMT_FBPicks.py : Plot first break picks overlain the real data.

TOMO2D_2_TOY2DAC_GEOM.py : Convert a TOMO2D geometry file to TOY2DAC format.

TOY2DACFILES.py : Generate the frequency domain data for FWI from the time domain OBS data.

TOY2DAC_MAKE_MODELS.py : Construct TOY2DAC formatted velocity models from a TOMO2D formatted velocity model.

TXIN_2_TOMO2D_GEOM.py : Construct a TOMO2D geometry file from the output of first break picking with PLOTSEC.

TOY2DAC_2_TOMO2D_GEOMETRY.py : Convert a TOY2DAC geometry file to TOMO2D format for forward modelling, ray tracing, etc.

TOY2DAC_2_TOMO2D_MODEL.py : Convert a TOY2DAC formatted velocity file to TOMO2D format.
