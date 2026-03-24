function [kdata,k_in,k_out,seq_args] = lps_convert_data(fname, ofile)
% gets data from scanarchive file, and formats the kspace trajectories and
% sequence arguments based on seq_args.mat, then writes formatted data to
% h5 file for external recon if specified
% by David Frey (djfrey@umich.edu)
%
% inputs:
% ifile - name of input file to read in (.h5 for GE scanarchive or .dat for
%   Siemens TWIX)
% ofile - name of output h5 file to write formatted data to (leave empty
% to not write to file)
%
% outputs:
% kdata - kspace data (ndat x nint x nprj x nrep x ncoil)
% k_in - spoke-in kspace trajectory (ndat x nint x nprj x nrep x 3)
% k_out - spoke-out kspace trajectory (ndat x nint x nprj x nrep x 3)
% seq_args - struct containing pulse sequence arguments
%

    % get directory and file names
    fname = which(fname);
    if isempty(fname)
        error('file %s not found', fname);
    end
    [dir,~,ext] = fileparts(fname);

    % load in sequence arguments (generated from lps_write_seq.m)
    seq_args = lpsutl.loadh5struct([dir,'/seq_args.h5']);

    switch ext
        case '.h5' % GE ScanArchive
            % load archive
            archive = GERecon('Archive.Load', fname);

            % skip past receive gain calibration TRs (pislquant)
            for n = 1:seq_args.pislquant
                currentControl = GERecon('Archive.Next', archive);
            end

            % loop through shots
            for i = 1:seq_args.nprj*seq_args.nint*seq_args.nrep
                currentControl = GERecon('Archive.Next', archive);
                if i == 1
                    [ndat,nc] = size(currentControl.Data);
                    ndat = ndat/seq_args.nechoes;
                    kdata = zeros(ndat, seq_args.nechoes, nc, seq_args.nint, seq_args.nprj, seq_args.nrep);
                end
                [iint,iprj,irep] = ind2sub([seq_args.nint, seq_args.nprj, seq_args.nrep], i);
                kdata(:,:,:,iint,iprj,irep) = reshape(currentControl.Data,ndat,seq_args.nechoes,nc);
            end
            kdata = permute(kdata,[1,2,4:6,3]); % n x necho x nint x nprj x nrep x nc
        
        case '.dat' % SIEMENS Twix
            twix = mapVBVD(fname);
            kdata = twix{end}.image.unsorted();
            ndat = size(kdata,1); % data points/projection
            nc = size(kdata,2); % number of coils
            kdata = permute(kdata(:,:,seq_args.pislquant+1:end),[1,3,2]);
            kdata = reshape(kdata, ... % n x nint x nprj x nrep x nc
                ndat, seq_args.nechoes, seq_args.nint, seq_args.nprj, seq_args.nrep, nc);

        otherwise
            error('invalid file extension: %s', ext)
    end

    % generate kspace trajectory
    [~,~,~,~,k_in0,k_out0] = lpsutl.gen_lps_waveforms( ...
        'sys', seq_args.sys, ... % pulseq mr system object
        'fov', seq_args.fov, ... % fov (cm)
        'N', seq_args.N_nom, ... % nominal matrix size
        'nspokes', seq_args.nspokes, ... % number of lps spokes
        'nechoes', seq_args.nechoes, ... % number of echoes
        't_seg', seq_args.t_seg, ... % number of samples/segment
        't_rf', seq_args.t_rf, ... % number of samples/rf pulse
        'C', seq_args.C, ... % option to use satellite trajectories
        'fa', seq_args.fa ... % rf flip angle (deg)
        );
    
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
    k_in = repmat(k_in,[1,1,1,1,seq_args.nrep]); % [n*necho x 3 x nint x iprj x nrep]
    k_out = repmat(k_out,[1,1,1,1,seq_args.nrep]);
    k_in = permute(k_in,[1,3:5,2]); % [n*necho x nint x nprj x nrep x 3]
    k_out = permute(k_out,[1,3:5,2]);

    % save h5 file
    if nargin > 1 && ~isempty(ofile)

        if isfile(ofile)
            system(sprintf('rm %s',ofile));
        end

        % save kspace data
        lpsutl.saveh5struct(ofile, real(kdata), '/kdata/real');
        lpsutl.saveh5struct(ofile, imag(kdata), '/kdata/imag');

        % save sampling locations
        lpsutl.saveh5struct(ofile, k_in, '/ktraj/spoke_in');
        lpsutl.saveh5struct(ofile, k_out, '/ktraj/spoke_out');

        % save number of coils
        lpsutl.saveh5struct(ofile, nc, '/ncoil');

        % save sequence arguments
        lpsutl.saveh5struct(ofile, seq_args, '/seq_args');

    end

end

