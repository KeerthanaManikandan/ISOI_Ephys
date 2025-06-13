# Read Me
This repository contains all codes necessary to process and store RS data collected from ISOI and linear electrode array simultaneously.

## Main scripts
1. `ephysImaging_compilation_vFinal.m`: Master script that performs temporal, spatial, laminar analysis
2. `controlAnalysis.m`: Script that performs spatial and control analyses.
3. `ephysImaging_plottingFinal.m`: Plots all compiled data for spatial, temporal, laminar analysis and other figures.
4. `animalStateAnalysis.m`: Script to calculate edge frequencies and powers from EEG spectrogram. 
5. `crossModal_FC_correlations.m`:Script to correlate recording-wise FC map to average FC maps for two example runs (check supplementary figure 1) 
6. `compare_LongitudinalFCMaps.m`:


## Main functions
1. `getMonkeyParams_Imaging_Ephys.m`: Initializes all necessary variables such as file name,experiment dates, lens combination, blood vessel maps, file locations for one monkey.
2. `getAllData_Imaging_Ephys.m`: Stores/retrieves imaging and LFP data. This function also identifies bad time segments, channels, and transition channels into cortex.
3. `saveLFPSingleProbe.m`: Sub-function within `getAllData_Imaging_Ephys.m` to extract and store electrophysiological data
4. `getCrossCorrROI.m`: Performs cross-correlations between rs-ISOI and LFP, plots the temporal profile of cross-correlations, and also the ROI at peak negative correlations for all frequencies and layers.
5. `getCrossCorrFOV.m`: Generates and saves the cross-modal maps for all frequencies and layers.
6. `getInfraSlowPowerLFP.m`: Calculates infraslow powers from electrophysiological data. 

## Dependent functions
1. `getRSConnectivityMaps.m`:  Generates FC map averaged across multiple RS runs for a specific animal. 
2. `getMonkeyParamsRS.m`: Initializes variables such as experiment date, runs, reference, folder locations. 
3. `getPreProcessedDataRestingState.m`: Processes and stores rs-ISOI data from experiments involving rs-ISOI only. 
4. `getMasks.m`:Loads the masks for a specific run.
5. `calculateSeedSignal.m:** Calculates the mean seed signal to generate FC map.
6. `plotCorrMap.m`:Calculates the FC map for a seed given a seed signal.
7. `removeBadTimesFromSpec.m`: This function removes the bad times determined from the spectrogram which have not been removed previously.
8. `plotLagProfiles.m`: Plots the median cross-correlations as a function of lag and returns the peak negative correlations and lag. 

### Necessary folders
- 01_Sunin: To process imaging data.
- neuroshare: To process neural data recorded from Ripple
- nonlinear: To perform image registration.
- chronux_2_12: To calculate spectrogram 

### Abbreviations
- rs-ISOI: Resting state intrinsic signal optical imaging
- LFP: Local Field Potential
- FC: Functional Connectivity 
- RS: Resting state
- ROI: Region of interest
- FOV: Field of view