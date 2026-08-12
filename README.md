# nCAF Paper Code
## Overview
This repository contains the code and datasets used to generate the figures in the preprint [Millisecond-scale, single-neuron credit assignment in a songbird](https://www.biorxiv.org/content/10.64898/2026.07.31.740148v1). 

## System Requirements
This code is intended to run on MATLAB release R2025a or later (on any compatible operating system). It has been tested in MATLAB R2025a on Ubuntu 24.04 LTS. No non-standard hardware is required. 

The following MATLAB toolboxes are required:
* Curve Fitting Toolbox
* Image Processing Toolbox
* Parallel Computing Toolbox
* Signal Processing Toolbox
* Statistics and Machine Learning Toolbox

## Installing and Running the Code
No installation is required; this repository can be cloned and each script run in MATLAB. 
> [!WARNING]
> This repository relies on Git LFS to store preprocessed data. Make sure Git LFS is installed using:
> ```bash
> git lfs install
> ```

As this repository relies on submodules, cloning should be performed using the following command:
```bash
git clone --recurse-submodules https://github.com/FeeLab/ncaf_paper_code.git
```
If you already cloned the repo without submodules, initialize them using:
```bash
git submodule update --init --recursive
```
Note that as the repository contains the preprocessed data, this cloning process may take a few minutes. 

The scripts to generate the figures are found in the `figure_generation` directory. There is a corresponding script for each main and extended data figure. Quantative results quoted in the manuscript are printed to the terminal by the script for the associated figure. Note that although most scripts run quickly, the simulations in `fig3.m` are computationally intensive and may take >24 hours to perform on a typical machine. 

## Dependencies and External Submodules
This project uses the following external repositories as submodules:
* [SpikeGLX_Datafile_Tools](https://github.com/jenniferColonell/SpikeGLX_Datafile_Tools/tree/5b684c32991523866beab7875cd93490550802b3) - Used for interacting with SpikeGLX files.
* [npy-matlab](https://github.com/kwikteam/npy-matlab/tree/7e61a884f7f5cf4f2382174efc8fe8f82e6e7e10) - Used for loading numpy files in MATLAB.

This project also uses the following external files, stored in the `external` directory:
* [colorcet](https://colorcet.com/) - Used for perceptually uniform colormaps.
* [MATLAB SMOTE](https://www.mathworks.com/matlabcentral/fileexchange/75401-synthetic-minority-over-sampling-technique-smote) - Used for balancing classes in the LMAN localizer classifier.
* [ndSparse](https://www.mathworks.com/matlabcentral/fileexchange/29832-n-dimensional-sparse-arrays) - Used for creating >2-dimensional sparse arrays. 



