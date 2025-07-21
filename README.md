<img width="400" height="200" alt="lps_horizontal" src="https://github.com/user-attachments/assets/6893530c-6f89-4567-ae8b-7f54092dfa1b" />

# Harmonized Looping Star MRI (lps)

by David Frey (djfrey@umich.edu)

This repository contains vendor-agnostic pulse sequence programming and reconstruction for Looping Star MRI. This sequence has been extensively tested and used to acquire fMRI data on GE MR750 scanners. However, the sequence is still being developed for use on Siemens scanners.

## Table of contents
1. [Table of contents](#table-of-contents)
2. [LpS Theory](#lps-theory)
   - [Pulse sequence basics](#pulse-sequence-basics)
   - [Reconstruction model](#reconstruction-model)
4. [Pulse sequence development](#pulse-sequence-development)
   - [Getting started](#getting-started)
   - [Generating the LpS sequence](#generating-the-lps-sequence)
   - [3D GRE calibration sequence](#3d-gre-calibration-sequence)
5. [Reconstruction](#reconstruction)
   - [Getting started](#getting-started)
   - [Converting scanner data](#converting-scanner-data)
   - [Reconstruction demo](#reconstruction-demo)
6. [References](#references)

## LpS Theory

### Pulse sequence basics
Looping Star (LpS) is a silent MRI pulse sequence based on gradient recalled zero-echo-time (ZTE) radial imaging, which multiplexes the readout of multiple echoes in order to achieve efficient k-space sampling while maintaining low gradient amplitude and slew rate.

A single block (TR) of a basic LpS pulse sequence is shown in the figure below:


### Reconstruction model

For a single TR, the acquisition of LpS data can be modeled as:

$$b = (A_\text{in} + A_\text{out}) x + \epsilon$$

where:
- $$x \in \mathbb{C}^N$$: discrete volumetric image
- $$N \in \mathbb{Z}$$: number of voxels in single volumetric image
- $$b \in \mathbb{C}^{MQ}$$: acquired multi-channel k-space data containing overlapping echoes
- $$\epsilon \in \mathbb{C}^{MQ}$$: acquisition noise
- $$M \in \mathbb{Z}$$: number of acquired k-space samples per TR
- $$Q \in \mathbb{Z}$$: number of coil channels for parallel imaging
- $$A_\text{in},A_\text{out} \in \mathbb{C}^{MQ \times N}$$: multi-channel k-space encoding matrices for the in-spoke and out-spoke data, respectively, i.e.:

$$A_\text{in/out} = (I_Q \otimes F_\text{in/out})\ S$$

- $$I_Q \in \mathbb{R}^{Q \times Q}$$: rank $$Q$$ identity matrix
- $$F_\text{in/out} \in \mathbb{C}^{M \times N}$$: Fourier encoding matrix for the in- or out-spoke (typically substituted for an nuFFT operator), i.e.:

$$F_\text{in/out} = \exp{ (-j 2 \pi k_\text{in/out} r^T ) }$$

- $$k_\text{in/out} \in \mathbb{R}^{M \times 3}$$: non-cartesian 3D k-space sample coordinates for the in- or out-spoke
- $$r \in \mathbb{R}^{N \times 3}$$: cartesian image space sampling coordinates
- $$S \in \mathbb{C}^{NQ \times N}$$: sensitivity encoding operator, i.e.:

$$
S = \begin{bmatrix}
   \text{diag}(s_1) \\
   \text{diag}(s_2) \\
   \vdots \\
   \text{diag}(s_Q)
   \end{bmatrix}
$$

- $$s_q \in \mathbb{C}^N$$: volumetric coil sensitivity map of the $$q$$ th channel

Additionally, we often account for the ill-conditioning at the edges of k-space due to echo overlap by applying a Fermi-shaped filter to the data which downweights high frequencies, i.e.:

$$A_\text{in/out} = (I_Q \otimes (H_\text{in/out} F_\text{in/out}))\ S$$

- $$H_\text{in/out} \in \mathbb{C}^{M \times M}$$: Fermi-shaped k-space filter for the in- or out-spoke, i.e.:

$$H_\text{in/out} = \text{diag}\left( \frac{1}{1 + \exp{\left(2 \pi \frac{|k_\text{in/out}| - \alpha_\text{cutoff}k_\text{max}}{\alpha_\text{rolloff}k_\text{max}} \right)}} \right)$$

- $$|k_\text{in/out}| \in \mathbb{R}^M$$: non-cartesian 3D k-space sample coordinate radii for the in- or out-spoke
- $$k_\text{max} \in \mathbb{R}$$: maximum k-space sample coordinate radii
- $$\alpha_\text{cutoff}$$, $$\alpha_\text{rolloff} \in (0,1)$$: Fermi filter cutoff and rolloff parameters, respectively

## Pulse sequence development

### Getting started

#### Install the repository
Clone this repo from the command line using `git clone git@github.com:fmrifrey/lps.git` to install the code into your working directory. Then, from MATLAB, add the base directory to your path using `addpath ./lps`. Additionally, if you would like to run 3D GRE calibration sequences, you can also add the [3D GRE subfolder](gre3d/) using `addpath ./lps/gre3d`. More information can be found in the [3D GRE calibration sequence](#3d-gre-calibration-sequence) section.

#### Install required packages
Run the script [`lps_update_packages.m`](lps_update_packages.m). This will install or update several required packages into the working directory. One can pick and choose the packages to update by running the script section by section:

##### General:
| Package | Use |
| --- | --- |
| [Pulseq v1.5](https://github.com/pulseq/pulseq) | Vendor agnostic pulse sequence development |
| [MIRT](https://github.com/JeffFessler/mirt) | General matlab functionality |

##### GE development only:
| Package | Use |
| --- | --- |
| [PulCeq v7](https://github.com/HarmonizedMRI/PulCeq/tree/tv7) | Conversion of .seq files to .pge files |
| [Toppe v7](https://github.com/toppeMRI/toppe/tree/develop) | Conversion of .seq files to .pge files |
| \*Orchestra 2.1 or higher | Reading in raw GE data |

\*note: <i>Orchestra requires the GE SDK license and can be installed from [WeConnect](https://weconnect.gehealthcare.com/s/feed/0D53a00008pQ1Q8CAK). Once Orchestra has been installed locally, you can set its path as an environment variable called `ORCHESTRA_PATH_MATLAB` (i.e. on Linux, add the line `export ORCHESTRA_PATH_MATLAB="/path/to/orchestra"` to your `.bashrc` file).</i>

#### Run the demo script:
The script [`lps_psd_demo.m`](lps_psd_demo.m) contains demo code to generate a simple LpS dataset along with a 3D GRE calibration sequence on GE scanners.

### Generating the LpS sequence

To generate the LpS sequence files, run the function [`lps_write_seq.m`](lps_write_seq.m). This function can be ran without any input arguments to generate files based on GE scanner defaults.

### 3D GRE calibration sequence

### Converting scanner data

## Reconstruction

### Getting started

### Reconstruction demo

## References
