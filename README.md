<img width="400" height="200" alt="lps_horizontal" src="https://github.com/user-attachments/assets/6893530c-6f89-4567-ae8b-7f54092dfa1b" />

# Harmonized Looping Star MRI (lps)

by David Frey (djfrey@umich.edu)

This repository contains vendor-agnostic pulse sequence programming and reconstruction for Looping Star MRI. This sequence has been extensively tested and used to acquire fMRI data on GE MR750 scanners. However, the sequence is still being developed for use on Siemens scanners.

## Table of contents
1. [Table of contents](#table-of-contents)
2. [Pulse sequence development](#pulse-sequence-development)
   - [Getting started](#getting-started)
   - [Generating the LpS sequence](#generating-the-lps-sequence)
   - [3D GRE calibration sequence](#3d-gre-calibration-sequence)
3. [Reconstruction](#reconstruction)
   - [Getting started](#getting-started)
   - [Converting scanner data](#converting-scanner-data)
   - [Reconstruction demo](#reconstruction-demo)
4. [LpS Theory](#lps-theory)
   - [Pulse sequence basics](#pulse-sequence-basics)
   - [Reconstruction model](#reconstruction-model)
      - [Basic model](#basic-model)
      - [Echo-in/out filtering](#echo-inout-filtering)
      - [Volume-wise acquisition model](#volume-wise-acquisition-model)
      - [Solving the reconstruction problem](#solving-the-reconstruction-problem)
5. [References](#references)

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

### Generating the LpS sequence files

To generate the LpS sequence files, run the function [`lps_write_seq.m`](lps_write_seq.m). This function can be ran without any input arguments to generate files based on GE scanner defaults.

Several parameters can be set when running the function:
| Parameter name | Type | Description |
| --- | --- | --- |


### 3D GRE calibration sequence

### Converting scanner data

## Reconstruction

### Getting started

### Reconstruction demo

## LpS Theory

### Pulse sequence basics

Looping Star (LpS) is a silent MRI pulse sequence based on gradient recalled zero-echo-time (ZTE) radial imaging, which multiplexes the readout of multiple echoes in order to achieve efficient k-space sampling while maintaining low gradient amplitude and slew rate.

A single block (TR) of a basic LpS pulse sequence is shown in the figure below:
<img width="1000" height="500" alt="lps_tr" src="https://github.com/user-attachments/assets/20f5e807-1702-4753-8586-c24addb6e136" />

Each RF pulse in the first half of the sequence excites a population of spins, which are then dephased to the outer edge of k-space, $$k_\text{max}$$, before another population is excited by the next RF pulse. Each RF pulse excites a single spoke in k-space.


#### Overlapping echoes

When $$\frac{N}{\text{FOV}} > k_{\text{max}}$$, some of signal encoded from the echo moving into the sampling region (in-spoke) will overlap with the signal encoded by the out-spoke. If not accounted for, this will result in high signal from the center of k-space appearing on the outside of k-space. It is also known that this signal overlap results in ill-conditioning of the encoding matrix. However, if correctly accounted for, this can allow one to improve the image resolution without needing to acquire more data.

### Reconstruction model

#### Basic model

The acquisition of the $$p^\text{th}$$ 3D k-space projection can be modeled as:

$$b_p = (A_{\text{in},p} + A_{\text{out},p})\ x(p\text{ TR}) + \epsilon_p, \quad p = 0,1,...VP-1$$

where:
- $$x(t): \mathbb{R} \to \mathbb{C}^{N^3}$$: discrete volumetric image at time $$t$$
- $$N^3 \in \mathbb{Z}$$: number of <i>reconstructed</i> voxels in single volumetric image
- $$b_p \in \mathbb{C}^{MQ}$$: acquired multi-channel k-space data containing overlapping echoes from the $$p^\text{th}$$ 3D k-space projection
- $$P \in \mathbb{Z}$$: number of 3D k-space projections per volume
- $$V \in \mathbb{Z}$$: total number of acquired volumes
- $$\epsilon_p \in \mathbb{C}^{MQ}$$: acqusition noise from the $$p^\text{th}$$ 3D k-space projection
- $$M \in \mathbb{Z}$$: number of acquired k-space samples per 3D k-space projection
- $$Q \in \mathbb{Z}$$: number of coil channels for parallel imaging
- $$A_{\text{in},p},A_{\text{out},p} \in \mathbb{C}^{MQ \times N^3}$$: multi-channel k-space encoding matrices for the in-spoke and out-spoke data of the $$p^\text{th}$$ 3D k-space projection, respectively, i.e.:

$$A_{\text{in/out},p} = (I_Q \otimes F_{\text{in/out},p})\ S$$

- $$I_Q \in \mathbb{R}^{Q \times Q}$$: rank $$Q$$ identity matrix
- $$F_{\text{in/out},p} \in \mathbb{C}^{M \times N^3}$$: \*Fourier encoding matrix for the in- or out-spoke of the $$p^\text{th}$$ 3D k-space projection, i.e.:

$$F_{\text{in/out},p} \approx \exp{ (-j 2 \pi\ k_{\text{in/out},p}\ r^T ) }$$

\*note: <i>the true DTFT matrix is typically substituted for an non-uniform FFT (nuFFT) operator to improve speed, hence the approximation</i>

- $$k_{\text{in/out},p} \in \mathbb{R}^{M \times 3}$$: non-cartesian 3D k-space sample coordinates for the in- or out-spoke of the $$p^\text{th}$$ 3D k-space projection, i.e.:

$$k_{\text{in/out},p} = k_{\text{in/out},0} R_p^T$$

- $$R_p \in \mathbb{R}^{3 \times 3}$$: 3D golden angle rotation matrix for $$p^\text{th}$$ 3D k-space projection
- $$r \in \mathbb{R}^{N^3 \times 3}$$: 3D cartesian image space sampling coordinates
- $$S \in \mathbb{C}^{N^3Q \times N^3}$$: sensitivity encoding operator, i.e.:

$$
S = \begin{bmatrix}
   \text{diag}(s_1) \\
   \text{diag}(s_2) \\
   \vdots \\
   \text{diag}(s_Q)
   \end{bmatrix}
$$

- $$s_q \in \mathbb{C}^{N^3}$$: discrete volumetric coil sensitivity map of the $$q^\text{th}$$ channel

#### Echo-in/out filtering

Additionally, we often account for the ill-conditioning at the edges of k-space due to echo overlap by applying a Fermi-shaped filter to the data which downweights the signal high frequencies, i.e.:

$$A_{\text{in/out},p} = (I_Q \otimes (H_\text{in/out} F_{\text{in/out},p}))\ S$$

where:
- $$H_\text{in/out} \in \mathbb{C}^{M \times M}$$: Fermi-shaped k-space filter for the in- or out-spoke, i.e.:

$$H_\text{in/out} = \text{diag}\left( \frac{1}{1 + \exp{\left(2 \pi \frac{|k_{\text{in/out},0}| - \mu_\text{cutoff}k_\text{max}}{\sigma_\text{rolloff}k_\text{max}} \right)}} \right)$$

- $$|k_{\text{in/out},0}| \in \mathbb{R}^M$$: non-cartesian 3D k-space sample coordinate radii for the in- or out-spoke of the $$0^\text{th}$$ k-space projection
- $$k_\text{max} \in \mathbb{R}$$: maximum k-space sample coordinate radii
- $$\mu_\text{cutoff},\sigma_\text{rolloff} \in (0,1)$$: Fermi filter cutoff and rolloff parameters, respectively

#### Volume-wise acquisition model

For multi-volumetric imaging, the acquisition can be modeled by using a block matrix representation which bins together $$P$$ projections per volume and reconstructs each volume as a seperate volumetric image:

$$\mathbf{b} = \mathbf{A} \mathbf{x} + \mathbf{\epsilon}$$

where:
- $$\mathbf{x} \in \mathbb{C}^{N^3V}$$: discrete volumetric image timeseries, i.e.:

$$
\mathbf{x} =
\begin{bmatrix}
   x(0) \\
   x(P\text{ TR}) \\
   \vdots \\
   x((V-1)P\text{ TR})
\end{bmatrix}
$$

- $$\mathbf{b} \in \mathbb{C}^{MQPV}$$: multi-volume, multi-projection, multi-channel k-space data, i.e.:

$$
\mathbf{b} =
\begin{bmatrix}
   \begin{bmatrix}
      b_0 \\
      b_1 \\
      \vdots \\
      b_{P-1}
   \end{bmatrix} \\
   \begin{bmatrix}
      b_P \\
      b_{P+1} \\
      \vdots \\
      b_{2P-1}
   \end{bmatrix} \\
   \vdots \\
   \begin{bmatrix}
      b_{(V-1)P} \\
      b_{(V-1)P+1} \\
      \vdots \\
      b_{VP-1}
   \end{bmatrix}
\end{bmatrix}
$$

- $$\mathbf{\epsilon} \in \mathbb{C}^{MQPV}$$: acquisition noise
- $$\mathbf{A} \in \mathbb{C}^{MQPV \times N^3V}$$: multi-volume, multi-projection, multi-channel acquisition model, i.e.:

$$
\mathbf{A} = 
\begin{bmatrix}
   \begin{bmatrix}
      A_{\text{in},0} + A_{\text{out},0} \\
      A_{\text{in},1} + A_{\text{out},1} \\
      \vdots \\
      A_{\text{in},P-1} + A_{\text{out},P-1}
   \end{bmatrix} & \mathbf{0} \\
   \mathbf{0} & \begin{bmatrix}
      A_{\text{in},P} + A_{\text{out},P} \\
      A_{\text{in},P+1} + A_{\text{out},P+1} \\
      \vdots \\
      A_{\text{in},2P-1} + A_{\text{out},2P-1}
   \end{bmatrix} \\
   && \ddots \\
   &&& \begin{bmatrix}
      A_{\text{in},(V-1)P} + A_{\text{out},(V-1)P} \\
      A_{\text{in},(V-1)P+1} + A_{\text{out},(V-1)P+1} \\
      \vdots \\
      A_{\text{in},VP-1} + A_{\text{out},VP-1}
   \end{bmatrix} \\
\end{bmatrix}
$$

\*note: <i>here, we assume linearity of the encoding matrices in order to vertically concatenate the projection-wise encoding matrices into each volume. However, in the case of nuFFT approximation which is non-linear, each volume-wise operator must be constructed by first concatenating the k-space sampling coordinates of each projection when computing the nuFFT operator for each volume</i>

#### Solving the reconstruction problem

We can solve for the best image in terms of $\ell2$ norm error by solving the following Least Squares problem:

$$ \mathbf{x}_* = \underset{\mathbf{x}}{\text{argmin}} \text{ ||}\mathbf{A} \mathbf{x} - \mathbf{b}\text{||}_2^2 $$

If the total number of k-space samples is less than the number of reconstructed voxels (i.e. $$MPQ \lt N^3$$), the problem above is underdetermined, and has infinitely many solutions of the form:

$$ \mathbf{x}_* = \mathbf{A}^+ \mathbf{b} + \mathbf{z}$$

where:
- $$\mathbf{A}^+ \in \mathbb{C}^{N^3V \times MQPV}$$: right Moore-Penrose psuedoinverse of $$\mathbf{A}$$
- $$\mathbf{z} \in \mathcal{N}(\mathbf{A})$$: any vector in the $$(N^3 - MQP)V$$-dimensional nullspace of $$\mathbf{A}$$

In other words, since k-space is undersampled, there are infinitely many ways to fill in the missing frequency components. The nullspace of $$\mathbf{A}$$ represents the missing frequency components. $$\mathbf{A}^+ \mathbf{b}$$ is called the minimum-norm solution, which zero-fills the missing frequency components. This is the solution that contains aliasing which gradient-based methods will converge to. However, with the right choice of $$z$$, the true, un-aliased image is also an equally likely solution in terms of $$\ell2$$ error.

Similarly, we can still acquire enough samples such that $$MPQ \geq N^3$$ without fully representing all the needed frequency components. This is especially the case in LpS when the edges of 3D k-space are undersampled while the origin is oversampled. In this case, the problem is fully or overdetermined, but ill-conditioned. The problem has an exact solution, $$\mathbf{x}_* = \mathbf{A}^+ \mathbf{b}$$, but is most likely to be dominated by the overfitting of noise. Other solutions that are equally as likely (such as the true image) only differ in data consistency due to subtle differences in the noise profile. To better condition our problem, we can promote or penalize certain solutions based on prior assumptions about the true image by adding a regularization term:

$$ \mathbf{x}_* = \underset{\mathbf{x}}{\text{argmin}} \text{ ||}\mathbf{A} \mathbf{x} - \mathbf{b}\text{||}_2^2 + g(\mathbf{x})$$

where:
- $$g(x): \mathbb{C}^{N^3V} \to \mathbb{R}$$: regularization function

For example, if $$g(x) = \text{||}(F \otimes I_{N^3})x\text{||}_1$$ where $$F \in \mathbb{C}^{V \times V}$$ is a discrete temporal Fourier Transform operator, solutions which have a sparse temporal frequency representation will be promoted.

The problem can then be solved using gradient based methods (i.e. CG, FISTA, PGM), which depend on the choice of regularization function, its convexity/linearity and smoothness.


## References
