# Arousal waves
MATLAB code and data files related to ["Global waves synchronize the brain's functional systems with fluctuating arousal."](https://advances.sciencemag.org/content/7/30/eabf2709)

## Overview
* **dependencies** folder includes code for reading and writing CIFTI files, pulled from the FieldTrip toolbox
* **surface_files** folder includes Conte69 atlas files for standard 32k and 4k (for optical flow analysis) meshes
* **supporting_files** folder includes:
  * Handy structure-specific masks for neocortex, thalamus, striatum, cerebellum, and hippocampus, based on the standard 91,282 grayordinate representation
  * List of subject IDs used in the paper (originally used in [Chen et al. 2020 Neuroimage](https://www.sciencedirect.com/science/article/pii/S1053811920301944))
  * Yeo 7-network parcellation spanning cerebral cortex, thalamus, striatum, and cerebellum (see paper for relevant citations)
* **output_files** folder includes:
  * HCP_RV_phasemap.mat - vertex/voxelwise maps of phase shifts in relation to respiratory time series (averaged across subjects and sessions)
  * RV_dynamics_sm3_zm.mat - the average BOLD signal changes mapped to a canonical respiratory variation cycle
  * FC_diffusion_coords.mat - a single file that contains principal functional connectivity coordinates in neocortex (from [Margulies et al. 2016 PNAS](https://www.pnas.org/content/113/44/12574)), thalamus, striatum, and cerebellum, as shown in the paper

## Additional dependencies (all freely available)
### Software packages
[Chronux](http://chronux.org/) \
[HCP Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench) \
[SUIT](http://www.diedrichsenlab.org/imaging/suit.htm) (for cerebellum flat map visualization only)
### Datasets
[HCP S1200 release](https://www.humanconnectome.org/study/hcp-young-adult/article/s1200-group-average-data-release) \
[NeuroTycho](http://neurotycho.org/anesthesia-and-sleep-task)

## Getting started

Add files in "supporting_files" to path

Begin with physio_phase_mapping.m, which outputs the phase-locking value result as well as intermediate files to use in other scripts (e.g., coherence_spectra.m)

All code for monkey analyses is contained within scripts having "monkey" prefix



**Please direct any questions to ryanvraut@gmail.com**
