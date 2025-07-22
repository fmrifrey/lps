function [kdata,k_in,k_out,seq_args] = lps_convert_data(safile, h5file)
% gets data from scanarchive file, and formats the kspace trajectories and
% sequence arguments based on seq_args.mat, then writes formatted data to
% h5 file for external recon if specified
% by David Frey (djfrey@umich.edu)
%
% inputs:
% safile - name of scanarchive file to read in
% h5file - name of output h5 file to write formatted data to (leave empty
% to not write to file)
%
% outputs:
% kdata - kspace data (ndat x nint x nprj x nrep x ncoil)
% k_in - spoke-in kspace trajectory (ndat x nint x nprj x nrep x 3)
% k_out - spoke-out kspace trajectory (ndat x nint x nprj x nrep x 3)
% seq_args - struct containing pulse sequence arguments
%

    % get directory and file names
    d = dir(safile);
    sadir = d(1).folder;
    safile = d(1).name;

    % load in sequence arguments
    seq_args = lpsutl.loadh5struct([sadir,'/seq_args.h5']);

    % load archive
    archive = GERecon('Archive.Load', [sadir,'/',safile]);

    % skip past receive gain calibration TRs (pislquant)
    for n = 1:seq_args.pislquant
        currentControl = GERecon('Archive.Next', archive);
    end
    
    % loop through shots
    for i = 1:seq_args.nprj*seq_args.nint*seq_args.nrep
        currentControl = GERecon('Archive.Next', archive);
        if i == 1
            [ndat,nc] = size(currentControl.Data);
            kdata = zeros(ndat, nc, seq_args.nint, seq_args.nprj, seq_args.nrep);
        end
        [iint,iprj,irep] = ind2sub([seq_args.nint, seq_args.nprj, seq_args.nrep], i);
        kdata(:,:,iint,iprj,irep) = currentControl.Data;
    end
    kdata = permute(kdata,[1,3:5,2]); % n x nint x nprj x nrep x nc

    % generate kspace trajectory
    [~,~,~,k_in0,k_out0] = lpsutl.gen_lps_waveforms( ...
        'sys', seq_args.sys, ... % pulseq mr system object
        'fov', seq_args.fov, ... % fov (cm)
        'N', seq_args.N_nom, ... % nominal matrix size
        'nspokes', seq_args.nspokes, ... % number of lps spokes
        't_seg', seq_args.t_seg, ... % number of samples/segment
        't_rf', seq_args.t_rf, ... % number of samples/rf pulse
        'fa', seq_args.fa ... % rf flip angle (deg)
        );

    % convert to 3D
    k_in0 = padarray(k_in0,[0,1],0,'post');
    k_out0 = padarray(k_out0,[0,1],0,'post');
    
    % rotate the kspace trajectory & excitation grads
    k_in = zeros(ndat,3,seq_args.nint,seq_args.nprj);
    k_out = zeros(ndat,3,seq_args.nint,seq_args.nprj);
    for iprj = 1:seq_args.nprj
        for iint = 1:seq_args.nint
            R = lpsutl.rot_3dtga(iprj,iint);
            k_in(:,:,iint,iprj) = k_in0 * R';
            k_out(:,:,iint,iprj) = k_out0 * R';
        end
    end

    % repeat for repetitions
    k_in = repmat(k_in,[1,1,1,1,seq_args.nrep]); % [n x 3 x nint x iprj x nrep]
    k_out = repmat(k_out,[1,1,1,1,seq_args.nrep]);
    k_in = permute(k_in,[1,3:5,2]); % [n x nint x nprj x nrep x 3]
    k_out = permute(k_out,[1,3:5,2]);

    % save h5 file
    if nargin > 1 && ~isempty(h5file)

        if isfile(h5file)
            system(sprintf('rm %s',h5file));
        end

        % save kspace data
        lpsutl.saveh5struct(h5file, real(kdata), '/kdata/real');
        lpsutl.saveh5struct(h5file, imag(kdata), '/kdata/imag');

        % save sampling locations
        lpsutl.saveh5struct(h5file, k_in, '/ktraj/spoke_in');
        lpsutl.saveh5struct(h5file, k_out, '/ktraj/spoke_out');

        % save number of coils
        lpsutl.saveh5struct(h5file, nc, '/ncoil');

        % save sequence arguments
        lpsutl.saveh5struct(h5file, seq_args, '/seq_args');

    end

end

