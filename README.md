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

A single block (TR) of a basic LpS pulse sequence is shown in the figure below. 

### Reconstruction model

## Pulse sequence development

### Getting started

#### Install the repository
Clone this repo from the command line using `git clone git@github.com:fmrifrey/lps.git` to install the code into your working directory. Then, from MATLAB, add the base directory to your path using `addpath ./lps`. Additionally, if you would like to run 3D GRE calibration sequences, you can also add the [3D GRE subfolder](gre3d/) using `addpath ./lps/gre3d`. More information can be found in the [3D GRE calibration sequence](#3d-gre-calibration-sequence) section.

#### Install required packages
Run the script [`update_packages.m`](update_packages.m). This will install or update several required packages into the working directory. One can pick and choose the packages to update by running the script section by section:

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
The script [`psd_demo.m`](psd_demo.m) contains demo code to generate a simple LpS dataset along with a 3D GRE calibration sequence on GE scanners.

### Generating the LpS sequence

To generate the LpS sequence files, run the function [`lps_write_seq.m`](lps_write_seq.m). This function can be ran without any input arguments to generate files based on GE scanner defaults.

### 3D GRE calibration sequence

### Converting scanner data

## Reconstruction

### Getting started

### Reconstruction demo

## References
