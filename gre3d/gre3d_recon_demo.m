%% set arguments
basedir = './'; % directory containing data
fname_kdata = 'raw_data.h5'; % name of input .h5 file (in basedir)
fname_smaps = 'smaps.h5'; % name of input (or output) smaps .h5 file
fname_out = sprintf('recon_%s.h5',...
    string(datetime('now','Format','yyyyMMddHHmm'))); % name of output recon .h5 file (in basedir)

% set recon parameters
rec_args.estimate_smap = 1; % option to estimate sensitivity maps from ACS data
rec_args.Q = 8; % number of coils to compress to
rec_args.niter = 30; % number of iterations for CG
rec_args.beta = 0.5; % regularization parameter (spatial roughness)

%% load gre3d data from g5 file
fprintf('loading raw data...\n');
str = lpsutl.loadh5struct(fullfile(basedir,fname_kdata));
kdata = str.kdata.real + 1i*str.kdata.imag;
msk = str.msk > 0;
seq_args = str.seq_args;

%% set up Fourier encoding operator
fprintf('setting up fft...\n');
F = fatrix2('idim', seq_args.N*ones(1,3), ...
    'odim', seq_args.N*ones(1,3), ...
    'omask', msk, ...
    'forw', @(~,x)lpsutl.fftc(x,[],1:3), ...
    'back', @(~,x)lpsutl.ifftc(x,[],1:3) ...
    );

%% load or estimate sense maps
if rec_args.estimate_smap

    fprintf('estimating sensitivity maps...\n');
    
    % get ACS data
    if mod(seq_args.Nacs,2) % odd ACS lines
        acs_idcs = floor(seq_args.N)/2 + (-(seq_args.Nacs-1)/2:(seq_args.Nacs-1)/2);
    else % even
        acs_idcs = floor(seq_args.N)/2 + (-seq_args.Nacs/2:seq_args.Nacs/2-1);
    end
    acs_data = kdata(acs_idcs,acs_idcs,acs_idcs,:);

    % estimate sensitivity maps from ACS data
    smaps = mri_sensemap_denoise(lpsutl.ifftc(acs_data,[],1:3)); % MIRT method - CG
    % smaps = bart('ecalib -b0 -m1', acs_data) % BART method - eSPIRIT (recommended)

    % save the sensitivity maps
    lpsutl.saveh5struct(fname_smaps,real(smaps),'/real');
    lpsutl.saveh5struct(fname_smaps,real(smaps),'/imag');

else

    fprintf('loading sensitivity maps...\n');

    % load the sensitivity maps
    s = lpsutl.loadh5struct(fname_smaps);
    smaps = s.real + 1i*s.imag;

end

% upsample smaps
smaps = lpsutl.resample3d(smaps,seq_args.N*ones(1,3));

%% perform coil compression and create SENSE operator
fprintf('coil compressing data...\n');
[kdata,~,Vr] = ir_mri_coil_compress(kdata,'ncoil',rec_args.Q);
smaps = reshape(reshape(smaps,[],size(smaps,4))*Vr,[seq_args.N*ones(1,3),rec_args.Q]);

fprintf('constructing SENSE operator...\n');
Ss = cell(rec_args.Q,1);
for i = 1:rec_args.Q
    smapi = smaps(:,:,:,i);
    Ss{i} = Gdiag(smapi);
end
S = lpsutl.block_col(Ss,0);

%% create the system matrix and regularization function
A = lpsutl.kronI(rec_args.Q,F,0)*S;
C = sqrt(rec_args.beta) * lpsutl.finite_diff(A.idim,1:3);
x0 = zeros(seq_args.N*ones(1,3));

%% solve with CG
fprintf('solving with CG...\n');
A_pcg = fatrix2( ...
    'idim', [prod(A.idim),1], ...
    'odim', [prod(A.odim),1], ...
    'forw', @(~,x) reshape(A*reshape(x,A.idim),[],1), ...
    'back', @(~,y) reshape(A'*reshape(y,A.odim),[],1) ...
    );
x_star = qpwls_pcg1(x0(:), A_pcg, 1, kdata(:), C, 'niter', rec_args.niter);
img_gre3d = reshape(x_star,A.idim);

%% save to h5 recon file
fprintf('saving to h5 file...\n');
fname = fullfile(basedir,fname_out);
if exist(fname, 'file')
    delete(fname);
end
lpsutl.saveh5struct(fname, real(img_gre3d), '/sol/real');
lpsutl.saveh5struct(fname, imag(img_gre3d), '/sol/imag');
lpsutl.saveh5struct(fname, seq_args, '/seq_args');
lpsutl.saveh5struct(fname, rec_args, '/rec_args');
fprintf('reconstruction completed --> saved to %s\n',fname_out);