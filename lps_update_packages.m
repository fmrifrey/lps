% gets the current compatible packages for running code from this repo
% by David Frey (djfrey@umich.edu)



%% Pulseq
fprintf('updating pulseq... ')
system('[ -d "./pulseq" ] && rm -rf ./pulseq');
system('git clone --branch v1.5.0 git@github.com:pulseq/pulseq.git 2> /dev/null');
addpath pulseq/matlab
fprintf('done.\n')

%% MIRT
fprintf('updating MIRT... ')
system('[ -d "./MIRT" ] && rm -rf ./MIRT');
system('git clone git@github.com:JeffFessler/MIRT.git 2> /dev/null');
run ./MIRT/setup.m
run ./MIRT/ir_mex_build.m
fprintf('done.\n')

%% ISMRMRD
fprintf('updating ismrmrd... ')
system('[ -d "./ismrmrd" ] && rm -rf ./ismrmrd');
system('git clone git@github.com:ismrmrd/ismrmrd.git 2> /dev/null');
addpath ./ismrmrd/matlab
fprintf('done.\n')

%% PulCeq (for GE pulse sequence development)
fprintf('updating PulCeq... ')
system('[ -d "./PulCeq" ] && rm -rf ./PulCeq');
system('git clone --branch tv7 git@github.com:HarmonizedMRI/PulCeq.git 2> /dev/null');
addpath PulCeq/matlab
fprintf('done.\n')

%% toppe (for GE pulse sequence development)
fprintf('updating toppe... ')
system('[ -d "./toppe" ] && rm -rf ./toppe');
system('git clone --branch develop git@github.com:fmrifrey/toppe.git 2> /dev/null');
addpath toppe
fprintf('done.\n')

%% Orchestra (for reading in data from GE scanners)
fprintf('adding orchestra to path... ')
addpath(getenv('ORCHESTRA_PATH_MATLAB')); % ORCHESTRA_PATH_MATLAB environment variable should point to orchestra folder
fprintf('done.\n')

%% mapVBVD (for reading in data from Siemens scanners)
fprintf('updating mapVBVD... ')
system('[ -d "./mapVBVD" ] && rm -rf ./mapVBVD');
system('git clone git@github.com:pehses/mapVBVD.git 2> /dev/null');
addpath mapVBVD
fprintf('done.\n')