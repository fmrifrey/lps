<img width="400" height="200" alt="lps_horizontal" src="https://github.com/user-attachments/assets/6893530c-6f89-4567-ae8b-7f54092dfa1b" />


# Harmonized Looping Star MRI (lps)

by David Frey (djfrey@umich.edu)

This repository contains vendor-agnostic pulse sequence programming and reconstruction for Looping Star MRI. This sequence has been extensively tested and used to acquire fMRI data on GE MR750 scanners. However, <b>this sequence is still being developed for use on Siemens scanners</b>.

## Table of contents
1. [Table of contents](#table-of-contents)
2. [Getting started](#getting-started)
3. [Pulse sequence development](#pulse-sequence-development)
   - [PSD demo script](#psd-demo-script)
   - [Generating the LpS sequence files](#generating-the-lps-sequence-files)
   - [3D GRE calibration sequence](#3d-gre-calibration-sequence)
5. [Reconstruction](#reconstruction)
   - [Converting scanner data](#converting-scanner-data)
   - [Recon demo script](#reconstruction-demo)
   - [Reconstructing the calibration data (SENSE map estimation)](#reconstructing-the-calibration-data-sense-map-estimation)
7. [References](#references)


## Getting started

### Install the repository
Clone this repo from the command line using `git clone git@github.com:fmrifrey/lps.git` to install the code into your working directory. Then, from MATLAB, add the base directory to your path using `addpath ./lps`. Additionally, if you would like to run 3D GRE calibration sequences, you can also add the [3D GRE subfolder](gre3d/) using `addpath ./lps/gre3d`. More information can be found in the [3D GRE calibration sequence](#3d-gre-calibration-sequence) section.

### Install required packages
Run the script [`lps_update_packages.m`](lps_update_packages.m). This will install or update several required packages into the working directory. One can pick and choose the packages to update by running the script section by section:

#### General:
| Package | Use |
| --- | --- |
| [Pulseq v1.5](https://github.com/pulseq/pulseq) | Vendor agnostic pulse sequence development |
| [MIRT](https://github.com/JeffFessler/mirt) | General matlab functionality |
| [PISCO](https://github.com/ralobos/PISCO) | Sensitivity map estimation |

#### GE development only:
| Package | Use |
| --- | --- |
| [PulCeq v7](https://github.com/HarmonizedMRI/PulCeq/tree/tv7) | Conversion of .seq files to .pge files |
| [Toppe v7](https://github.com/toppeMRI/toppe/tree/develop) | Conversion of .seq files to .pge files |
| \*Orchestra 2.1 or higher | Reading in raw GE data |

\*note: <i>Orchestra requires the GE SDK license and can be installed from [WeConnect](https://weconnect.gehealthcare.com/s/feed/0D53a00008pQ1Q8CAK). Once Orchestra has been installed locally, you can set its path as an environment variable called `ORCHESTRA_PATH_MATLAB` (i.e. on Linux, add the line `export ORCHESTRA_PATH_MATLAB="/path/to/orchestra"` to your `.bashrc` file).</i>

#### Siemens development only:
| Package | Use |
| --- | --- |
| [mapVBVD](https://github.com/pehses/mapVBVD) | Reading in raw Siemens data |

## Pulse sequence development

#### PSD demo script:
The script [`lps_psd_demo.m`](lps_psd_demo.m) contains demo code to generate a simple LpS dataset along with a 3D GRE calibration sequence on GE scanners.

### Generating the LpS sequence files

To generate the LpS sequence files, run the function [`lps_write_seq.m`](lps_write_seq.m). This function can be ran without any input arguments to generate files based on GE scanner defaults.

Several parameters can be set when running [`lps_write_seq.m`](lps_write_seq.m):

| Parameter name | Type | Default | Description |
| --- | --- | --- | --- |
| `sys` | struct | GE system | Pulseq system structure (see [`get_scanner_defaults.m`](+lpsutl/get_scanner_defaults.m)) |
| `dir` | string | current directory | output directory name |
| `tr` | float | 100 | repetition time (ms) |
| `fov` | float | 20 | field of view (cm) |
| `N_nom` | int | 128 | nominal matrix size (i.e. $$k_\text{max} = \frac{N_\text{nom}}{\text{fov}}$$) |
| `dummyshots` | int | 20 | number of dummy shots at beginning of sequence (to reach steady state) |
| `nrep` | int | 1 | number of repeated projections (i.e. $$PV = n_\text{rep} \times n_\text{int} \times n_\text{prj}$$) |
| `nint` | int | 1 | number of in-plane rotations (interleaves) (i.e. $$PV = n_\text{rep} \times n_\text{int} \times n_\text{prj}$$) |
| `nprj` | int | 16 | number of thru-plane rotations (i.e. $$PV = n_\text{rep} \times n_\text{int} \times n_\text{prj}$$) |
| `nechoes` | int | 16 | number of echoes to acquire (GREs only, not including FID) |
| `nspokes` | int | 23 | number of spokes/rf excitations per plane |
| `psd_rf_wait` | float | 100e-6 | RF–gradient delay (s), GE scanner-specific |
| `psd_grd_wait` | float | 100e-6 | ADC–gradient delay (s), GE scanner-specific |
| `coil` | str | 'xrm' | GE coil model (for PNS check) |
| `dwell` | float | 4 | ADC sample rate (us) |
| `t_seg` | float | 1120 | spoke length (us) |
| `t_rf` | float | 16 | rf pulse width (us) |
| `fa` | float | 4 | rf flip angle (degrees) |
| `C` | $$L \times 3$$ matrix; float | default 2D looping star coefficients | k-space trajectory Fourier series basis coefficient matrix (for research purposes only) |
| `plotseq` | bool | 0 | option to plot the sequence |
| `pislquant` | int | 1 | number of TRs to acquire during prescan |
| `writepge` | bool | 1 | option to write out .pge file to run on GE scanners |

Running [`lps_write_seq.m`](lps_write_seq.m) will create a few files in the directory specified by `dir`:
- `lps.seq`: .seq file for Pulseq
- `lps.pge`: .pge file for GE interpreter (if specified by `writepge`)
- `seq_args.h5`: .h5 file containing the sequence arguments for reconstruction

### 3D GRE calibration sequence

This repository also contains code for a 3D GRE pulse sequence. To generate the 3D GRE sequence files, run the function [`gre3d/gre3d_write_seq.m`](gre3d/gre3d_write_seq.m). This function can be ran without any input arguments to generate files based on GE scanner defaults.

Several parameters can be set when running [`gre3d/gre3d_write_seq.m`](gre3d/gre3d_write_seq.m):

| Parameter name | Type | Default | Description |
| --- | --- | --- | --- |
| `sys` | struct | GE system | Pulseq system structure (see [`get_scanner_defaults.m`](+lpsutl/get_scanner_defaults.m)) |
| `dir` | string | current directory | output directory name |
| `te` | float | min | echo time (ms), use 'min' for minimum |
| `tr` | float | 30 | repetition time (ms), use 'min' for minimum |
| `fov` | float | 20 | field of view (cm) |
| `slabfrac` | float | 0.7 | excitation slab width (fraction of fov) |
| `fa` | float | 4 | rf flip angle (degrees) |
| `rfspoil` | bool | 1 | option to do rf phase spoiling |
| `fatsat` | bool | 1 | option to do fat saturation |
| `fatsatChemShift` | float | 3.5 | off resonance frequency of fat (ppm) |
| `N` | int | 128 | 3D matrix size |
| `Nacs` | int | 32 | width of fully sampled (ACS) region at center of k-space |
| `Ry` | int | 2 | Y CAIPI acceleration factor outside of ACS region |
| `Rz` | int | 2 | Z CAIPI acceleration factor outside of ACS region |
| `delta` | int | 1 | CAIPI odd/even shift |
| `peorder` | string | snake | phase encode ordering scheme (snake or spiral) |
| `dwell` | float | 4 | ADC sample rate (us) |
| `dt` | float | 20e-6 | ADC sampling rate (s) |
| `dummyshots` | int | 20 | number of dummy shots at beginning of sequence (to reach steady state) |
| `writepge` | bool | 1 | option to write out .pge file to run on GE scanners |
| `pislquant` | int | 1 | number of TRs to acquire during prescan |
| `plotseq` | bool | 0 | option to plot the sequence |
| `psd_rf_wait` | float | 100e-6 | RF–gradient delay (s), GE scanner-specific |
| `psd_grd_wait` | float | 100e-6 | ADC–gradient delay (s), GE scanner-specific |
| `coil` | str | 'xrm' | GE coil model (for PNS check) |

Running [`gre3d/gre3d_write_seq.m`](gre3d/gre3d_write_seq.m) will create a few files in the directory specified by `dir`:
- `gre3d.seq`: .seq file for Pulseq
- `gre3d.pge`: .pge file for GE interpreter (if specified by `writepge`)
- `seq_args.h5`: .h5 file containing the sequence arguments for reconstruction


## Reconstruction

### Converting scanner data

This repository works with custom, readable h5 file formats which makes reconstructing the data externally much easier without having to rely on vendor-specific recon functions.

[`lps_convert_data.m`](lps_convert_data.m) will convert the raw scanner data (scanarchive from GE, or twix file from Siemens) into the desired .h5 file format for reconstruction. The function will read in the specified raw vendor file and attempt to locate the corresponding `seq_args.h5` file in the same directory. The output .h5 file will have the following heirarchical structure:

```
lps raw data structure
   ├─── ncoil  % integer describing total number of acquired channels
   ├─── kdata  % acquired k-space data (M x nint x nprj x nrep x ncoil)
   |      ├─── real
   |      └─── imag
   ├─── ktraj  % k-space trajectory (M x nint x nprj x nrep x 3)
   |      ├─── spoke_in
   |      └─── spoke_out
   └─── seq_args  % sequence arguments structure
```

The gre3d sequence has a similar function, [`gre3d/gre3d_convert_data.m`](gre3d/gre3d_convert_data.m). The output .h5 file will have the following heirarchical structure:

```
gre3d raw data structure
   ├─── ncoil  % integer describing total number of acquired channels
   ├─── kdata  % acquired k-space data (N x N x N x ncoil)
   |      ├─── real
   |      └─── imag
   ├─── msk  % k-space sampling locations (N x N x N)
   └─── seq_args  % sequence arguments structure
```

### Reconstruction demo

[`lps_recon_demo.m`](lps_recon_demo.m) contains demo code for reconstructing image data from an lps raw data file using regularized CG-SENSE with a spatio-temporal roughness penalty. Several parameters control the reconstruction, and can be set at the top of the script:

| Parameter name | Type | Default | Description |
| --- | --- | --- | --- |
| `fname` | string | ./rawdata.h5 | file name for converted raw data |
| `fname_smaps` | string | ./smaps.h5 | name of input .h5 smaps file (see [section below](#reconstructing-the-calibration-data) for info on calculating SENSE maps using PISCO) |
| `Q` | int | 4 | number of coils to compress to |
| `N` | int | empty | reconstruction matrix size (leave empty to use N_nom) |
| `niter` | int | 30 | number of iterations of CG |
| `ints2use` | array, int | empty | number of interleaves to use in reconstruction (leave empty to use all) |
| `ints2use` | array, int | empty | number of interleaves to use in reconstruction (leave empty to use all) |
| `prjs2use` | array, int | empty | number of projections to use in reconstruction (leave empty to use all) |
| `reps2use` | array, int | empty | number of repetitions to use in reconstruction (leave empty to use all) |
| `P` | int | empty | number of projections to bin per volume (leave empty to use all ints and prjs) |
| `dcf_init` | bool | 0 | option to initialize solution with density-compensated adjoint operation |
| `use_parfor` | bool | 1 | option to parallelize over block computations (sense-maps and volumes) |
| `fermi_cutoff` | float (range 1-2) | 1 | fermi voxel basis function cutoff (frac of nominal resolution, 1 corresponds to no overlapping echoes, 2 corresponds to fully overlapping) |
| `fermi_rolloff` | float (range 0-1) | 0.1 | fermi voxel basis function rolloff (frac of nominal resolution) |
| `beta` | float | 2^16 | Tikhonov regularization parameter |
| `debug` | bool | 0 | toggle debug mode |

### Reconstructing the calibration data (+SENSE map estimation)

[`gre3d/gre3d_recon_demo.m`](gre3d/gre3d_recon_demo.m) contains demo code for reconstructing image data from a gre3d raw data file using regularized CG-SENSE with a spatial roughness penalty. Several parameters control the reconstruction, and can be set at the top of the script:

| Parameter name | Type | Default | Description |
| --- | --- | --- | --- |
| `fname` | string | ./raw_data.h5 | name of input .h5 data file |
| `fname_smaps` | string | ./smaps.h5 | name of input or output .h5 smaps file |
| `estimate_smap` | bool | 1 | option to estimate sensitivity maps using [PISCO](https://mr.usc.edu/download/pisco/) |
| `Q` | int | 8 | number of coils to compress to |
| `niter` | int | 30 | number of iterations of CG |
| `beta` | float | 0.5 | Tikhonov regularization parameter |

## References

[1] Wiesinger, F., Menini, A., & Solana, A. B. (2019). Looping star. Magnetic resonance in medicine, 81(1), 57-68.

[2] Xiang, H., Fessler, J. A., & Noll, D. C. (2024). Model‐based reconstruction for looping‐star MRI. Magnetic resonance in medicine, 91(5), 2104-2113.

[3] Lobos, R. A., Chan, C. C., & Haldar, J. P. (2023). New theory and faster computations for subspace-based sensitivity map estimation in multichannel MRI. IEEE transactions on medical imaging, 43(1), 286-296.