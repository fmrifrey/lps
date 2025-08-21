%% set arguments
basedir = './'; % directory containing data
fname_kdata = 'raw_data.h5'; % name of input .h5 file (in basedir)
fname_smaps = '../smaps.h5'; % name of input smaps .h5 file (in basedir)
fname_out = sprintf('recon_%s.h5',...
    string(datetime('now','Format','yyyyMMddHHmm'))); % name of output recon .h5 file (in basedir)

% set recon parameters
rec_args.Q = 6; % number of coils to compress to
rec_args.mu_cutoff = 0.85; % kspace cutoff for echo-in/out filtering
rec_args.sig_rolloff = 0.1; % kspace rolloff for echo-in/out filtering
rec_args.beta = 2^18; % spatial roughness penalty
rec_args.beta_t = 0; % temporal roughness penalty
rec_args.niter = 30; % number of iterations for CG
rec_args.N = []; % reconstruction matrix size (leave empty for cutoff * N)
rec_args.ints2use = []; % number of interleaves to use (leave empty for all)
rec_args.prjs2use = []; % number of projections to use (leave empty for all)
rec_args.reps2use = []; % number of repetitions to use (leave empty for all)
rec_args.P = 32; % number of projections per volume (leave empty for all)
rec_args.initdcf = 0; % option to initialize with density compensated recon
rec_args.usepar = 1; % option to parallelize block matrix computations

%% load LpS data from h5 file
fprintf('loading raw data...\n');
str = lpsutl.loadh5struct(fullfile(basedir,fname_kdata));
kdata = str.kdata.real + 1i*str.kdata.imag;
k_in = str.ktraj.spoke_in;
k_out = str.ktraj.spoke_out;
seq_args = str.seq_args;
if isempty(rec_args.N)
    rec_args.N = 2*ceil(rec_args.mu_cutoff*seq_args.N_nom/2);
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
if isempty(rec_args.P)
    rec_args.P = rec_args.ints2use*rec_args.prjs2use; % each rep is a vol
end
[Fs_in,Fs_out,b] = lpsutl.lps_setup_nuffts(kdata,k_in,k_out,seq_args, ...
    'N', rec_args.N, ...
    'ints2use', rec_args.ints2use, ...
    'prjs2use', rec_args.prjs2use, ...
    'reps2use', rec_args.reps2use, ...
    'volwidth', rec_args.P ...
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
    seq_args.N_nom/rec_args.N*rec_args.mu_cutoff, ... % kspace filter cutoff
    seq_args.N_nom/rec_args.N*rec_args.sig_rolloff ... % kspace filter rolloff
    );

%% load in the sensitivity maps and coil compress
fprintf('loading sensitivity maps...\n');
str = lpsutl.loadh5struct(fullfile(basedir,fname_smaps));
smaps = str.real + 1i*str.imag;

%% compress data and get compression matrix
fprintf('coil compressing data...\n');
[tmp,~,Vr] = ir_mri_coil_compress(permute(b,[1,3,2]),'ncoil',rec_args.Q);
b = permute(tmp,[1,3,2]);

% upsample smaps
smaps = lpsutl.resample3d(smaps,rec_args.N*ones(1,3));

% coil compress the smaps
smaps = reshape(reshape(smaps,[],size(smaps,4))*Vr,[rec_args.N*ones(1,3),rec_args.Q]);

%% create the sensitivity encoding operator
fprintf('constructing SENSE operator...\n');
if rec_args.Q > 1
    Ss = cell(rec_args.Q,1);
    for i = 1:rec_args.Q
        smapi = smaps(:,:,:,i);
        Ss{i} = Gdiag(smapi(Fs_in{1}.imask(:)),'mask',Fs_in{1}.imask);
    end
    S = lpsutl.block_col(Ss,0);
end

%% create the system matrix
fprintf('constructing full system matrix...\n');
As_in = cell(nvol,1);
As_out = cell(nvol,1);
for i = 1:nvol % repeat forward encoding for each sense map
    if rec_args.Q > 1
        As_in{i} = lpsutl.kronI(rec_args.Q, Hs_in{i}*Fs_in{i}, par_coils) * S;
        As_out{i} = lpsutl.kronI(rec_args.Q, Hs_out{i}*Fs_out{i}, par_coils) * S;
    else
        As_in{i} = Hs_in{i}*Fs_in{i};
        As_out{i} = Hs_out{i}*Fs_out{i};
    end
end
if nvol == 1 % single volume
    A_in = As_in{1};
    A_out = As_out{1};
else % each volume is a block in a block diagonal matrix
    A_in = lpsutl.block_diag(As_in,par_vols);
    A_out = lpsutl.block_diag(As_out,par_vols);
end
A = A_in + A_out; % sum of echo-in/out system matrices

%% initialize the solution
if rec_args.initdcf % initialize with dcf-nuFFT solution
    fprintf('computing dcf-nufft initialization...\n');    
    WAs_in = cell(nvol,1);
    WAs_out = cell(nvol,1);
    if par_vols
        parfor i = 1:nvol % repeat forward encoding for each sense map
            W_in = lpsutl.dcf_pipe(Fs_in{i},3,par_vols);
            W_out = lpsutl.dcf_pipe(Fs_out{i},3,par_vols);
            WAs_in{i} = lpsutl.kronI(rec_args.Q, W_in * Hs_in{i}*Fs_in{i}, par_coils) * S;
            WAs_out{i} = lpsutl.kronI(rec_args.Q, W_out * Hs_out{i}*Fs_out{i}, par_coils) * S;
        end
    else
        for i = 1:nvol % repeat forward encoding for each sense map
            W_in = lpsutl.dcf_pipe(Fs_in{i},3,par_vols);
            W_out = lpsutl.dcf_pipe(Fs_out{i},3,par_vols);
            WAs_in{i} = lpsutl.kronI(rec_args.Q, W_in * Hs_in{i}*Fs_in{i}, par_coils) * S;
            WAs_out{i} = lpsutl.kronI(rec_args.Q, W_out * Hs_out{i}*Fs_out{i}, par_coils) * S;
        end
    end
    if nvol == 1 % single volume
        WA_in = WAs_in{1};
        WA_out = WAs_out{1};
    else % each volume is a block in a block diagonal matrix
        WA_in = lpsutl.block_diag(WWAs_in,par_vols);
        WA_out = lpsutl.block_diag(As_out,par_vols);
    end
    WA = WA_in + WA_out;
    x0 = WA' * b;
    x0 = ir_wls_init_scale(WA,b,x0);
else % initialize with zeros
    fprintf('initializing solution with zeros...\n');
    x0 = zeros(A.idim);
end

%% make the regularizer (spatial roughness penalty)
fprintf('creating regularization function...\n');
if nvol == 1 % spatial regularizer only
    C = sqrt(rec_args.beta) * lpsutl.finite_diff(A.idim,1:3);
else % spatial + temporal
    Cblocks = cell(4,1);
    for i = 1:3 % spatial
        Cblocks{i} = sqrt(rec_args.beta) * lpsutl.finite_diff(A.idim,i);
    end
    Cblocks{4} = sqrt(rec_args.beta_t) * lpsutl.finite_diff([A.idim,1],4);
    C = lpsutl.block_col(Cblocks,0);
end

%% solve with CG
fprintf('solving with CG...\n');
% create system operator compatible with qpwls_pcg1 function
A_pcg = fatrix2( ...
    'idim', [prod(A.idim),1], ...
    'odim', [prod(A.odim),1], ...
    'forw', @(~,x) reshape(A*reshape(x,[A.idim,1]),[],1), ...
    'back', @(~,y) reshape(A'*reshape(y,[A.odim,1]),[],1) ...
    );
x_star = qpwls_pcg1(x0(:), A_pcg, 1, b(:), 0, 'niter', rec_args.niter);
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