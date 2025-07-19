%% set arguments
basedir = './'; % directory containing data
fname_kdata = 'raw_data.h5'; % name of input .h5 file (in basedir)
fname_smaps = '../smaps.h5'; % name of input smaps .h5 file (in basedir)
fname_out = sprintf('recon_%s.h5',...
    string(datetime('now','Format','yyyyMMddHHmm'))); % name of output recon .h5 file (in basedir)

% set recon parameters
rec_args.ncoil_comp = 8; % number of coils to compress to
rec_args.cutoff = 0.85; % kspace cutoff for echo-in/out filtering
rec_args.rolloff = 0.1; % kspace rolloff for echo-in/out filtering
rec_args.beta = 2^18; % regularization parameter for quadratic finite differencing penalty
rec_args.niter = 30; % number of iterations for CG
rec_args.M = []; % reconstruction matrix size (leave empty for cutoff * N)
rec_args.ints2use = []; % number of interleaves to use (leave empty for all)
rec_args.prjs2use = []; % number of projections to use (leave empty for all)
rec_args.reps2use = []; % number of repetitions to use (leave empty for all)
rec_args.volwidth = 32; % number of projections per volume (leave empty for all)
rec_args.initdcf = 0; % option to initialize with density compensated recon
rec_args.usepar = 1; % option to parallelize block matrix computations

%% load LpS data from h5 file
fprintf('loading raw data...\n');
str = lpsutl.loadh5struct(fullfile(basedir,fname_kdata));
kdata = str.kdata.real + 1i*str.kdata.imag;
k_in = str.ktraj.spoke_in;
k_out = str.ktraj.spoke_out;
seq_args = str.seq_args;
if isempty(rec_args.M)
    rec_args.M = 2*ceil(rec_args.cutoff*seq_args.N/2);
end

%% set up Fourier encoding operators and data
fprintf('setting up nufft operators...\n');
if isempty(rec_args.ints2use)
    rec_args.ints2use = seq_args.nint; % use all unique in-plane rotations
end
if isempty(rec_args.prjs2use)
   rec_args.prjs2use = seq_args.nprj; % use all unique thru-plane rotations
end
if isempty(rec_args.reps2use)
    rec_args.reps2use = seq_args.nrep; % use all repetitions
end
if isempty(rec_args.volwidth)
    rec_args.volwidth = rec_args.ints2use*rec_args.prjs2use; % each rep is a vol
end
[Fs_in,Fs_out,b] = lpsutl.lps_setup_nuffts(kdata,k_in,k_out,seq_args, ...
    'M', rec_args.M, ...
    'ints2use', rec_args.ints2use, ...
    'prjs2use', rec_args.prjs2use, ...
    'reps2use', rec_args.reps2use, ...
    'volwidth', rec_args.volwidth ...
    );
nvol = size(b,3);

%% determine dimensions of parallelization
if nvol > 1 && rec_args.usepar % parallelize over volumes
    par_vols = 1;
    par_coils = 0;
elseif nvol == 1 && rec_args.usepar % parallelize over coils
    par_vols = 0;
    par_coils = 1;
else % no parallelization
    par_vols = 0;
    par_coils = 0;
end

%% create the kspace echo-in/out Fermi filters
fprintf('setting up echo-in/out filters...\n');
[Hs_in,Hs_out] = lpsutl.lps_setup_filters(Fs_in,Fs_out, ...
    seq_args.N/rec_args.M*rec_args.cutoff, ... % kspace filter cutoff
    seq_args.N/rec_args.M*rec_args.rolloff ... % kspace filter rolloff
    );

%% load in the sensitivity maps and coil compress
fprintf('loading sensitivity maps...\n');
str = lpsutl.loadh5struct(fullfile(basedir,fname_smaps));
smaps = str.real + 1i*str.imag;

% compress data and get compression matrix
fprintf('coil compressing data...\n');
[tmp,~,Vr] = ir_mri_coil_compress(permute(b,[1,3,2]),'ncoil',rec_args.ncoil_comp);
b = permute(tmp,[1,3,2]);

% upsample smaps
smaps = lpsutl.resample3d(smaps,rec_args.M*ones(1,3));

% coil compress the smaps
smaps = reshape(reshape(smaps,[],size(smaps,4))*Vr,[rec_args.M*ones(1,3),rec_args.ncoil_comp]);

%% create the sensitivity encoding operator
fprintf('constructing SENSE operator...\n');
Ss = cell(rec_args.ncoil_comp,1);
for i = 1:rec_args.ncoil_comp
    smapi = smaps(:,:,:,i);
    Ss{i} = Gdiag(smapi(Fs_in{1}.imask(:)),'mask',Fs_in{1}.imask);
end
S = lpsutl.block_col(Ss,par_coils);

%% create the system matrix
fprintf('constructing full system matrix...\n');
As_in = cell(nvol,1);
As_out = cell(nvol,1);
for i = 1:nvol
    As_in{i} = kronI(rec_args.ncoil_comp, Hs_in{i} * Fs_in{i}) * S;
    As_out{i} = kronI(rec_args.ncoil_comp, Hs_out{i} * Fs_out{i}) * S;
end
if nvol == 1 % single volume
    A_in = As_in{1};
    A_out = As_out{1};
else
    A_in = lpsutl.block_diag(As_in,par_vols);
    A_out = lpsutl.block_diag(As_out,par_vols);
end
A = A_in + A_out; % sum of echo-in/out system matrices

%% initialize the solution
if rec_args.initdcf % initialize with dcf-nuFFT solution
    fprintf('computing dcf-nufft initialization...\n');
    WAs_in = cell(nvol,1);
    WAs_out = cell(nvol,1);
    if nvol == 1 % single volume
        W_in = lpsutl.dcf_pipe(Fs_in{1});
        W_out = lpsutl.dcf_pipe(Fs_out{1});
        WAs_in{1} = kronI(rec_args.ncoil_comp, W_in * Hs_in{1} * Fs_in{1}) * S;
        WAs_out{1} = kronI(rec_args.ncoil_comp, W_out * Hs_out{1} * Fs_out{1}) * S;
        WA_in = WAs_in{1};
        WA_out = WAs_out{1};
    elseif par_vols
        % compute dcfs in parallel
        parfor i = 1:nvol
            Wi_in = lpsutl.dcf_pipe(Fs_in{i});
            Wi_out = lpsutl.dcf_pipe(Fs_out{i});
            WAs_in{i} = kronI(rec_args.ncoil_comp, Wi_in * Hs_in{i} * Fs_in{i}) * S;
            WAs_out{i} = kronI(rec_args.ncoil_comp, Wi_out * Hs_out{i} * Fs_out{i}) * S;
        end
        % construct new weighted system matrix
        WA_in = lpsutl.block_diag(WAs_in{:},1);
        WA_out = lpsutl.block_diag(WAs_out{:},1);
    else % same but not in parallel
        for i = 1:nvol
            Wi_in = lpsutl.dcf_pipe(Fs_in{i});
            Wi_out = lpsutl.dcf_pipe(Fs_out{i});
            WAs_in{i} = kronI(rec_args.ncoil_comp, Wi_in * Hs_in{i} * Fs_in{i}) * S;
            WAs_out{i} = kronI(rec_args.ncoil_comp, Wi_out * Hs_out{i} * Fs_out{i}) * S;
        end
        WA_in = lpsutl.block_diag(WAs_in{:},0);
        WA_out = lpsutl.block_diag(WAs_out{:},0);
    end
    WA = WA_in + WA_out; % sum of echo-in/out system matrices
    x0 = WA' * b;
    x0 = ir_wls_init_scale(A,b,x0);
else % initialize with zeros
    fprintf('initializing solution with zeros...\n');
    x0 = zeros(A.idim);
end

%% make the regularizer (spatial roughness penalty)
fprintf('creating regularization function...\n');
qp = Reg1(true(rec_args.M*ones(1,3)),'beta',rec_args.beta);
% % code to determine a good regularization parameter:
% qpwls_psf(Av1, qp.C, beta, true(seq_args.N*ones(1,3)),1, ...
%     'loop', 1, 'dx', seq_args.fov/seq_args.N, 'dz', seq_args.fov/seq_args.N);
C = qp.C;
if nvol > 1
    C = kronI(nvol,C);
end

%% solve with CG
fprintf('solving with CG...\n');
% create system operator compatible with qpwls_pcg1 function
A_pcg = fatrix2( ...
    'idim', [prod(A.idim),1], ...
    'odim', [prod(A.odim),1], ...
    'forw', @(~,x) reshape(A*reshape(x,A.idim),[],1), ...
    'back', @(~,y) reshape(A'*reshape(y,A.odim),[],1) ...
    );
x_star = qpwls_pcg1(x0(:), A_pcg, 1, b(:), C, 'niter', rec_args.niter);
img_lps = reshape(x_star,A.idim);

%% save to h5 recon file
fprintf('saving to h5 file...\n');
fname = fullfile(basedir,fname_out);
if exist(fname, 'file')
    delete(fname);
end
lpsutl.saveh5struct(fname, real(img_lps), '/sol/real');
lpsutl.saveh5struct(fname, imag(img_lps), '/sol/imag');
lpsutl.saveh5struct(fname, seq_args, '/seq_args');
lpsutl.saveh5struct(fname, rec_args, '/rec_args');
fprintf('reconstruction completed --> saved to %s\n',fname_out);