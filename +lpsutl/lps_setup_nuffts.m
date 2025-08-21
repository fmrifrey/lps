function [Fs_in,Fs_out,b] = lps_setup_nuffts(kdata,k_in,k_out,seq_args,varargin)
% sets up the volume-wise nufft objects given kspace data, trajectory, and
% desired indexing parameters
% by David Frey (djfrey@umich.edu)
%
% inputs:
% prjs2use - number of projections to use (leave empty for all)
% ints2use - number of interleaves to use (leave empty for all)
% reps2use - number of repetitions to use (leave empty for all)
% volwidth - rotations per volume (number of prjs/ints/reps to use per vol)
% N - reconstructed matrix size
% stride - number of rotations overlapping in each volume (not yet
% implemented)
% (also reads arguments from seq_args.mat file in current directory)
%
% outputs:
% Fs_in - cell array of nufft operators per volume (echo in)
%   if spokewise = true, Fs_in will be spokewise and framewise
% Fs_out - cell array of nufft operators per volume (echo out)
%   same applies for spokewise = true
% b - formatted kspace measurements
%
%
    
    % set default arguments
    arg.ints2use = [];
    arg.prjs2use = [];
    arg.reps2use = [];
    arg.volwidth = [];
    arg.N = seq_args.N_nom; % recon matrix size
    arg.rmspoke1 = true; % option to remove first spoke data
    arg.rmspokeN = true; % option to remove last spoke data
    arg.stride = 0; % still need to implement

    % parse arguments
    arg = vararg_pair(arg,varargin);
    
    % determine loop indicies to use
    if isempty(arg.ints2use)
        arg.ints2use = seq_args.nint; % use all unique in-plane rotations
    end
    if isempty(arg.prjs2use)
        arg.prjs2use = seq_args.nprj; % use all unique thru-plane rotations
    end
    if isempty(arg.reps2use)
        arg.reps2use = seq_args.nrep; % use all repetitions
    end

    % index desired loop indicies
    dt = seq_args.sys.gradRasterTime*1e6; % (us)
    nseg = round(seq_args.t_seg/dt);
    kdata = kdata(nseg*arg.rmspoke1+1:end-nseg*arg.rmspokeN,1:arg.ints2use,1:arg.prjs2use,1:arg.reps2use,:);
    k_in = k_in(nseg*arg.rmspoke1+1:end-nseg*arg.rmspokeN,1:arg.ints2use,1:arg.prjs2use,1:arg.reps2use,:);
    k_out = k_out(nseg*arg.rmspoke1+1:end-nseg*arg.rmspokeN,1:arg.ints2use,1:arg.prjs2use,1:arg.reps2use,:);

    % split data and trajectory into 3D volumes
    if isempty(arg.volwidth)
        arg.volwidth = arg.ints2use*args.prjs2use; % each rep is a vol
    end
    if mod(arg.ints2use*arg.reps2use*arg.prjs2use, arg.volwidth)
        error('total planes (%d) must be divisible by rpv (%d)', ...
            arg.ints2use*arg.reps2use*arg.prjs2use, arg.volwidth)
    else
        nvol = arg.ints2use*arg.reps2use*arg.prjs2use / arg.volwidth;
    end
    nc = size(kdata,5);
    kdata = reshape(kdata,[],nvol,nc);
    k_in = reshape(k_in,[],nvol,3);
    k_out = reshape(k_out,[],nvol,3);

    % convert trajectory to spatial frequencies
    omega_in = 2*pi*seq_args.fov/arg.N * k_in;
    omega_out = 2*pi*seq_args.fov/arg.N * k_out;

    % format the kspace measurements volume-wise
    b = reshape(kdata,[],nvol,nc);
    b = permute(b,[1,3,2]);

    % initialize volume-wise nufft operators
    Fs_in = cell(nvol,1);
    Fs_out = cell(nvol,1);

    % create spherical mask
    [Xgrd,Ygrd,Zgrd] = ndgrid(linspace(-1,1,arg.N), ...
        linspace(-1,1,arg.N), ...
        linspace(-1,1,arg.N));
    msk = (Xgrd.^2 + Ygrd.^2 + Zgrd.^2) <= 1;
    
    % set up volume-wise dcf objects
    nufft_args = {arg.N*ones(1,3), 6*ones(1,3), 2*arg.N*ones(1,3), ...
        arg.N/2*ones(1,3), 'table', 2^10, 'minmax:kb'};

    % loop through volumes
    for ivol = 1:nvol
        % get trajectory for current vol
        omegav_in = reshape(omega_in(:,ivol,:),[],3);
        omegav_out = reshape(omega_out(:,ivol,:),[],3);

        % assemble nufft object for current vol
        Fs_in{ivol} = Gnufft(msk, [omegav_in, nufft_args]);
        Fs_out{ivol} = Gnufft(msk, [omegav_out, nufft_args]);
    end

end
