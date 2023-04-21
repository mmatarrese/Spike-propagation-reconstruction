# Spike-propagation-reconstruction
This repository is part of the study "Spike Propagation Mapping Reveals Effective Connectivity and Predicts Surgical Outcome in Epilepsy" (Brain, 2023) by Matarrese M.A.G. et al. (https://doi.org/10.1093/brain/awad118)

The functions written for this project are made integrable with the open-source application Brainstorm (Tadel et al., 2011). Many of the files necessary for the correct execution of the code are Brainstorm outputs. Specifically, in this study, we used an Unconstrained Source Estimation with a dSPM kernel. The resulting matrix "source_map", to be given as input to this code, must be a Full Result matrix, i.e., kernel x data [3NxM] where N is the number of sources of the forward problem (solved with OpenMEEG) and M is the number of time samples of the analyzed window. The anatomical regions were also defined with Brainstorm. In particular, starting from the specific anatomy of the patient, we calculated an MNI parcellation (in this study we used the AAL3) and build the volumetric scouts on the forward model built for the localization of the sources.

contact: mag.matarrese@gmail.com
