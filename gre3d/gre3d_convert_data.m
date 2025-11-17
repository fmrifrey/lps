function [kdata,msk,seq_args] = gre3d_convert_data(vendor,infile,h5file)
% gets data from scanarchive file, and formats the kspace trajectories and
% sequence arguments based on seq_args.mat, then writes formatted data to
% h5 file for external recon if specified
% by David Frey (djfrey@umich.edu)
%
% inputs:
% vendor - 'ge' or 'siemens'
% infile - name of scanarchive file to read in
% h5file - name of output h5 file to write formatted data to (leave empty
% to not write to file)
%
% outputs:
% kdata - kspace data (Nx x Ny x Nz x ncoil)
% msk - sampling mask (Nx x Ny x Nz)
% seq_args - struct containing pulse sequence arguments
%

    % get directory and file names
    d = dir(infile);
    sadir = d(1).folder;
    infile = d(1).name;
    
    % load in sequence arguments
    seq_args = lpsutl.loadh5struct([sadir,'/seq_args.h5']);
    
    switch lower(vendor)
        case 'ge'
            % load archive
            archive = GERecon('Archive.Load', [sadir,'/',infile]);

            % get phase encode indicies
            if strcmpi(seq_args.peorder,'snake')
                pe_idcs = gre3dutl.snake_caipi_idcs(seq_args.N, seq_args.Ry, ...
                    seq_args.Rz, seq_args.delta, seq_args.Nacs);
            elseif strcmpi(seq_args.peorder,'spout')
                pe_idcs = gre3dutl.spout_caipi_idcs(seq_args.N, seq_args.Ry, ...
                    seq_args.Rz, seq_args.delta, seq_args.Nacs);
            else
                error('invalid option for peorder');
            end
            npe = length(pe_idcs);

            % skip past receive gain calibration TRs (pislquant)
            for n = 1:seq_args.pislquant
                currentControl = GERecon('Archive.Next', archive);
            end

            % loop through shots
            for i = 1:npe
                currentControl = GERecon('Archive.Next', archive);
                if i == 1
                    [~,nc] = size(currentControl.Data);
                    d1 = zeros(seq_args.N, nc, seq_args.N, seq_args.N);
                    msk = zeros(seq_args.N*ones(1,3));
                end
                msk(:,pe_idcs(i,1),pe_idcs(i,2)) = 1;
                d1(:,:,pe_idcs(i,1),pe_idcs(i,2)) = currentControl.Data;
            end

            % get kspace data
            kdata = permute(d1,[1,3,4,2]); % Nx x Ny x Nz x nc

        case 'siemens'
            twix = mapVBVD(infile);
            data = twix{2}.image();

    end

    % save h5 file
    if nargin > 1 && ~isempty(h5file)

        if isfile(h5file)
            system(sprintf('rm %s',h5file));
        end

        % save kspace data
        lpsutl.saveh5struct(h5file, real(kdata), '/kdata/real');
        lpsutl.saveh5struct(h5file, imag(kdata), '/kdata/imag');

        % save sampling mask
        lpsutl.saveh5struct(h5file, msk, '/msk');

        % save number of coils
        lpsutl.saveh5struct(h5file, nc, '/ncoil');

        % save sequence arguments
        lpsutl.saveh5struct(h5file, seq_args, '/seq_args');

    end

end
