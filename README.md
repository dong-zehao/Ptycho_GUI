# Ptycho_GUI README

## Overview

Ptycho_GUI is a MATLAB-based graphical user interface designed for performing ptychographic reconstructions from 4D-STEM datasets. The software facilitates data preprocessing, parameter selection, and iterative reconstruction using GPU acceleration. **Note: Ptycho_GUI is in a very early stage!**

**Mostly for personal/in-group use, maybe more features will be added in the future for generalization**

![figoverview](img\Overview.png)

## Requirements

- MATLAB (Tested with MATLAB 2023a and later)
- Required MATLAB toolboxes:
  - Parallel Computing Toolbox
  - Curve Fitting Toolbox
  - Image Processing Toolbox
  - Optimization Toolbox
  - Signal Processing Toolbox
- Ptychography scripts and dependencies are all included in the `Ptycho_GUI` directory

## Installation

1. Clone or download the `Ptycho_GUI` repository.
2. Open MATLAB and navigate to the `Ptycho_GUI\dist` folder.
3. Run `Ptycho_GUI.mlappinstall` to install Ptycho_GUI in MATLAB.

## Example Workflow

### Left panel: tunable parameters

1. **Path initialization**: Select appropriate paths and click "Test."
   - **Ptycho path**: path to `ptycho_GUI`
   - **Data path**: Choose the dataset folder containing 4D-STEM data (multiple datasets can be found).
2. **Prepare data**: 
   - Select your camera (detector) parameters (saved in `Ptycho_GUI\+configs\xxxx.mat`), including voltage, alpha, rbf, stepsize, rotation, ADU, which you can check and change in the `Show Hyper Param.` tab.
   - Select which dataset (scan number) to be prepared, and the corresponding defocus value (minus values represent over-focus). Note that multiple datasets can be prepared on one click.
   - Adjust binning and padding parameters for CBED, then click "Prepare."

3. **Set reconstruction parameters**: Define which dataset (scan number), thickness, layers, iterations, and other parameters. (details in below).
4. **Run reconstruction**: Click "Reconstruct" and monitor progress in the right panel.
5. **Adjust hyperparameters (if needed)**: Use the "Show Hyper Param." option to fine-tune settings.
6. **Analyze results**: Inspect the reconstructed phase image and Fourier error convergence plot.

### Right panel: Log outputs

1. **Log output**: Information outputs
2. **Figure outputs**: 
   - Fourier error vs. iteration number plot
   - Real-time reconstructed multi-modal probe
   - Real-time reconstructed phase image of object

## User Interface and Parameters Guide

### 1. Test Path and Camera

- **Ptycho path**: Select the directory to `ptycho_GUI`.
- **Data path**: Choose the dataset folder containing 4D-STEM datasets.
- **Detector**: Select the detector parameters saved in `+config\detectorname.mat`(e.g., K3_ARM300).
- **4D folder pattern**: Define the folder structure for loading multiple datasets.
- **File name**: filename for the 4D dataset.
- **Test button**: Verifies if the selected paths and camera settings are correct.

### 2. Prepare Data

- **Np_binto**: Set the binning factor for reducing dataset size.
- **Np_padto**: Define the padding size for diffraction patterns (padding zeros).
- **Rotate CBED?**: Choose whether to rotate the convergent beam electron diffraction (CBED) patterns.
- **Save dp?**: Select whether to save the `dp.hdf5` files
- **Scan number**: Defines the datasets to be prepared.
- **df**: Set the defocus value for each scan (minus values represent over-focus).
- **Prepare button**: Prepare the data for ptychographic reconstruction.

### 3. Ptycho Reconstruction

- **Scan number**: Defines the datasets to be reconstructed.
- **Total thickness**: sample thickness for multislice reconstruction.
- **N_layers**: Number of layers in multislice simulations.
- **N_iterations**: Number of iterations for the reconstruction.
- **Grouping**: batch size for reconstructions.
- **beta-LSQ**: similar to the learning rate in deep learning (1.0 is usually fine).
- **GPU id**: Choose the GPU device for reconstruction.
- **Np_presolve**: Set the region of CBED to be reconstructed (normally equals Np_padto).
- **Reconstruct button**: Begins the iterative reconstruction process.

### Log Output

- Displays real-time progress of the reconstruction, including iteration count, Fourier error, and estimated time remaining.
- Plots Fourier error convergence.
- Shows reconstructed phase image of the sample and amplitude image of the probe.

### Hyper-parameters (Optional)

Clicking "Show Hyper Param." reveals additional reconstruction settings:

- **Acc. Volt**: Acceleration voltage of the electron beam (e.g., 300 keV).
- **Alpha**: Convergence semi-angle in mrad.
- **rbf**: radius of bright-field disk in pixel.
- **Step size**: scanning step size in angstroms.
- **blur (PSF)**: Point spread function blur factor.
- **Probe mode**: Number of probe modes used.
- **reg_layer**: Regularization factor for multislice reconstuctions.
- **Remove ambiguity?**: Toggle for ambiguity removal in phase retrieval.
- **Niter probe upd.**: Number of iteration to start probe update.
- **Niter save**: Interval for saving intermediate results.
- **Niter pos. corr.**: Number of iteration to start position correction.
- **Niter plot**: Interval for plotting results in Log output panel.

## Contact

For questions or bug reports, please contact Zehao Dong at THU.