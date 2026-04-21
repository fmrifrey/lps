function [kdata,msk,seq_args] = gre3d_convert_data(fname,ofile)
% gets data from scanarchive file, and formats the kspace trajectories and
% sequence arguments based on seq_args.mat, then writes formatted data to
% h5 file for external recon if specified
% by David Frey (djfrey@umich.edu)
%
% inputs:
% fname - name of scanarchive file to read in
% ofile - name of output h5 file to write formatted data to (leave empty
% to not write to file)
%
% outputs:
% kdata - kspace data (Nx x Ny x Nz x ncoil)
% msk - sampling mask (Nx x Ny x Nz)
% seq_args - struct containing pulse sequence arguments
%

    % get directory and file names
    fname = which(fname);
    if isempty(fname)
        error('file %s not found', fname);
    end
    [dir,~,ext] = fileparts(fname);

    % load in sequence arguments
    seq_args = lpsutl.loadh5struct([dir,'/seq_args.h5']);

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
    msk = zeros(seq_args.N*ones(1,3));
    for i = 1:npe, msk(:,pe_idcs(i,1),pe_idcs(i,2)) = 1; end

    
    switch ext
        case '.h5' % GE ScanArchive
            % load archive
            archive = GERecon('Archive.Load', fname);
    
            % skip past receive gain calibration TRs (pislquant)
            for n = 1:seq_args.pislquant
                currentControl = GERecon('Archive.Next', archive);
            end

            % loop through shots
            for i = 1:npe
                currentControl = GERecon('Archive.Next', archive);
                if i == 1
                    [~,nc] = size(currentControl.Data);
                end
                d1(:,:,pe_idcs(i,1),pe_idcs(i,2)) = currentControl.Data;
            end
            
            % get kspace data
            kdata = permute(d1,[1,3,4,2]); % Nx x Ny x Nz x nc

        case '.dat' % Siemens Twix file
            twix = mapVBVD(fname);
            d1 = twix{end}.image.unsorted();
            nc = size(d1,2);
            kdata = zeros([size(msk),nc]);
            for i = 1:npe
                kdata(:,pe_idcs(i,1),pe_idcs(i,2),:) = reshape(d1(:,:,seq_args.pislquant+i), [], 1, 1, nc);
            end
    end

    % save h5 file
    if nargin > 1 && ~isempty(ofile)

        if isfile(ofile)
            system(sprintf('rm %s',ofile));
        end

        % save kspace data
        lpsutl.saveh5struct(ofile, real(kdata), '/kdata/real');
        lpsutl.saveh5struct(ofile, imag(kdata), '/kdata/imag');

        % save sampling mask
        lpsutl.saveh5struct(ofile, msk, '/msk');

        % save number of coils
        lpsutl.saveh5struct(ofile, nc, '/ncoil');

        % save sequence arguments
        lpsutl.saveh5struct(ofile, seq_args, '/seq_args');

    end

end
