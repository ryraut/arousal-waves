# Arousal waves
MATLAB code and data files related to Raut et al. 2021 Science Advances

## Dependencies (all freely available)
### Software packages:
Chronux: http://chronux.org/ \
[HCP Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench) \
[SUIT](http://www.diedrichsenlab.org/imaging/suit.htm) (for cerebellum flat map visualization only)

### Datasets:
[HCP S1200 release](https://www.humanconnectome.org/study/hcp-young-adult/article/s1200-group-average-data-release) \
[NeuroTycho](http://neurotycho.org/anesthesia-and-sleep-task)

## Overview:
* **Dependencies** folder includes code for reading and writing CIFTI files, pulled from the FieldTrip toolbox
* **surface_files** folder includes Conte69 atlas files for standard 32k and 4k (for optical flow analysis) meshes
* **supporting_files** folder includes:
  * Handy structure-specific masks for neocortex, thalamus, striatum, cerebellum, and hippocampus, based on the standard 91,282 grayordinate representation
  * List of subject IDs used in the paper (originally used in [Chen et al. 2020 Neuroimage](https://www.sciencedirect.com/science/article/pii/S1053811920301944))
  * Yeo 7-network parcellation spanning cerebral cortex, thalamus, striatum, and cerebellum (see paper for relevant citations)
* **output_files** folder includes:
  * vertex/voxelwise maps of phase shifts in relation to respiratory time series (averaged across subjects and sessions)
  * the average BOLD signal changes mapped to a canonical respiratory variation cycle
  * a single file that contains principal functional connectivity coordinates in neocortex (from [Margulies et al. 2016 PNAS](https://www.pnas.org/content/113/44/12574)), thalamus, striatum, and cerebellum, as shown in the paper

## Getting started

Add files in "supporting_files" to path

Begin with physio_phase_mapping.m, which outputs the phase-locking value result as well as intermediate files to use in other scripts (e.g., coherence_spectra.m)





**Please direct any questions to ryanvraut@gmail.com**
