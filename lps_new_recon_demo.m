%% set arguments
basedir = './'; % directory containing data
fname_kdata = './raw_data.h5'; % name of input .h5 file (in basedir)
fname_smaps = '../smaps.h5'; % name of input smaps .h5 file (in basedir)
fname_stmaps = './stmaps.h5'; % name of input STM .h5 file (in basedir)
fname_out = sprintf('recon_%s.h5',...
    string(datetime('now','Format','yyyyMMddHHmm'))); % name of output recon .h5 file (in basedir)

% set recon parameters
rec_args.Q = 4; % number of coils to compress to
rec_args.fermi_cutoff = 0.85; % kspace cutoff for echo-in/out filtering
rec_args.fermi_rolloff = 0.15; % kspace rolloff for echo-in/out filtering
rec_args.beta = 2^18; % tikhonov regularization parameter
rec_args.beta_unfold = 2e12; % UNFOLD regularization parameter
rec_args.stm_thresh = 0.0005; % STM eigenvalue threshold (to control L(x))
rec_args.niter = 30; % number of iterations for CG
rec_args.N = []; % reconstruction matrix size (leave empty for cutoff * N)
rec_args.ints2use = []; % number of interleaves to use (leave empty for all)
rec_args.prjs2use = 1; % number of projections to use (leave empty for all)
rec_args.reps2use = []; % number of repetitions to use (leave empty for all)
rec_args.P = 1; % number of projections per volume (leave empty for all)
rec_args.initdcf = 0; % option to initialize with density compensated recon
rec_args.usepar = 0; % option to parallelize block matrix computations

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
    'msk', 'sphere', ...
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

%% compress data and get compression matrix
fprintf('coil compressing data...\n');
[tmp,~,Vr] = ir_mri_coil_compress(permute(b,[1,3,2]),'ncoil',rec_args.Q);
b = permute(tmp,[1,3,2]);

%% load in the sensitivity maps and coil compress
if ~isempty(fname_smaps) && isfile(fname_smaps)
    fprintf('loading sensitivity maps...\n');
    str = lpsutl.loadh5struct(fullfile(basedir,fname_smaps));
    smaps = str.real + 1i*str.imag;
    
    % upsample smaps
    smaps = lpsutl.resample3d(smaps,rec_args.N*ones(1,3));
    
    % coil compress the smaps
    smaps = reshape(reshape(smaps,[],size(smaps,4))*Vr,[rec_args.N*ones(1,3),rec_args.Q]);

else % compress data to 1 coil if no smaps are used
    fprintf('no sensitivity maps found --> compressing data to 1 coil...\n');
    rec_args.Q = 1;
    tmp = ir_mri_coil_compress(permute(b,[1,3,2]),'ncoil',rec_args.Q);
    b = permute(tmp,[1,3,2]);
    smaps = ones(rec_args.N*ones(1,3));
end

%% create the sensitivity encoding operator
fprintf('constructing SENSE operator...\n');
S = fatrix2('idim', rec_args.N*ones(1,3), ...
    'odim', [rec_args.N*ones(1,3),rec_args.Q], ...
    'forw', @(~,x) x .* smaps, ... % forward: [S1; S2; ... SQ] * x
    'back', @(~,x) sum(x .* conj(smaps), 4) ... % adjoint: [S1' S2' ... SQ'] * x
    );

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
fprintf('initializing solution with zeros...\n');
x0 = zeros(A.idim);

%% create regularizers
I = fatrix2( ...
    'idim', [rec_args.N*ones(1,3),nvol], ...
    'odim', [rec_args.N*ones(1,3),nvol], ...
    'forw', @(~,x) x, ...
    'back', @(~,x) x ...
    );
PsiF = fatrix2( ...
    'idim', [rec_args.N*ones(1,3)], ...
    