%% GRE3D recon + sensitivity map estimation demo script
% by David Frey
%% set parameters
rec_args.fname = './raw_data.h5'; % input raw data .h5 file name (see gre3d_convert_data.m)
rec_args.fname_smaps = '../smaps.h5'; % smaps input file name
rec_args.estimate_smap = true; % option to estimate sensitivity maps from ACS data
rec_args.Q = 8; % number of coils to compress to for recon
rec_args.niter = 30; % number of iterations for CG
rec_args.beta = 0.5; % regularization parameter (spatial roughness)

%% load gre3d data from h5 file
fprintf('loading raw data...\n');
kdata = lpsutl.loadh5struct(rec_args.fname,'/kdata').real + ...
    1i*lpsutl.loadh5struct(rec_args.fname,'/kdata').imag; % raw kspace data
seq_args = lpsutl.loadh5struct(rec_args.fname,'/seq_args'); % sequence arguments
msk = lpsutl.loadh5struct(rec_args.fname).msk; % sampling mask
if isfile(rec_args.fname_smaps)
    smaps = lpsutl.loadh5struct(rec_args.fname_smaps).real + ...
        1i*lpsutl.loadh5struct(rec_args.fname_smaps).imag; % sensitivity maps
else
    smaps = []; % leave empty if no sensitity map file provided
end

%% estimate sense maps
if isempty(smaps) && rec_args.estimate_smap

    fprintf('estimating sensitivity maps...\n');
    
    % get ACS data
    if mod(seq_args.Nacs,2) % odd ACS lines
        acs_idcs = floor(seq_args.N)/2 + (-(seq_args.Nacs-1)/2:(seq_args.Nacs-1)/2);
    else % even
        acs_idcs = floor(seq_args.N)/2 + (-seq_args.Nacs/2:seq_args.Nacs/2-1);
    end
    acs_data = kdata(acs_idcs,acs_idcs,acs_idcs,:);

    % estimate sensitivity maps from ACS data
    smaps = PISCO_sensitivity_map_estimation(acs_data, seq_args.N*ones(1,3));

    % save the sensitivity maps
    lpsutl.saveh5struct(rec_args.fname_smaps,real(smaps),'/real');
    lpsutl.saveh5struct(rec_args.fname_smaps,imag(smaps),'/imag');

end

%% coil compression
if isempty(smaps)
    fprintf('no sensitivity maps provided -> compressing data...\n');
    rec_args.Q = 1;
    if ncoil > 1
        kdata = ir_mri_coil_compress(kdata,'ncoil',rec_args.Q);
    end
    smaps = ones([rec_args.N*ones(1,3),1]);
else
    fprintf('compressing data to %d virtual coils...\n', rec_args.Q);
    if isempty(rec_args.Q)
        rec_args.Q = size(kdata,4); % default - use all coils
    end

    % compress the data and sensitivity maps
    [kdata,~,Vt] = ir_mri_coil_compress(kdata,'ncoil',rec_args.Q);
    smaps = reshape(reshape(smaps, [], size(smaps,4)) * Vt, [size(smaps,1:3), rec_args.Q]);

    % resample sensitivity maps to recon size
    smaps = lpsutl.resample3d(smaps,seq_args.N*ones(1,3));
end

%% set up forward model
fprintf('setting up forward model...\n');
A = fatrix2('idim', seq_args.N*ones(1,3), ...
    'odim', [seq_args.N*ones(1,3),rec_args.Q], ...
    'forw', @(~,x) msk .* lpsutl.fftc(smaps .* x, 1:3), ...
    'back', @(~,y) numel(y)*sum(conj(smaps) .* lpsutl.ifftc(msk .* y, 1:3), 4) ...
    );

%% solve with CG
x0 = zeros(seq_args.N*ones(1,3));
x_star = qpwls_pcg1(x0, A, 1, kdata(:), 0, 'niter', rec_args.niter, 'stop_grad_tol', 1);