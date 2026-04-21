%% Looping Star demo recon script
% by David Frey
%% set parameters
rec_args.fname = './raw_data_sphere.h5'; % input raw data .h5 file name (see lps_convert_data.m)
rec_args.fname_smaps = '../smaps_sphere.h5'; % smaps input file name
rec_args.Q = 6; % number of compressed coils to use
rec_args.N = []; % recon image size
rec_args.echoes2use = []; % indices of echoes to include (empty = all)
rec_args.ints2use = []; % indices of interleaves to include (empty = all)
rec_args.prjs2use = []; % indices of projections to include (empty = all)
rec_args.reps2use = []; % indices of repetitions to include (empty = all)
rec_args.P = []; % number of projections to use per frame (empty = nint*nprj)
rec_args.niter = 30; % number of CG iterations
rec_args.dcf_init = true; % option to initialize solution with density compensated NUFFT
rec_args.use_parfor = true; % option to use parfor loop in frame/coil-wise NUFFTs
rec_args.fermi_cutoff = 1; % fermi voxel basis function cutoff (frac of nominal resolution bound)
rec_args.fermi_rolloff = 0.1; % fermi voxel basis function rolloff (frac of nominal resolution bound)
rec_args.beta = 2^14; % tikhonov regularization parameter
rec_args.debug = 0; % debug mode

%% load in data
fprintf('loading data...\n');
rec_args.fname = which(rec_args.fname);
kdata = lpsutl.loadh5struct(rec_args.fname,'/kdata').real + ...
    1i*lpsutl.loadh5struct(rec_args.fname,'/kdata').imag; % raw kspace data
ncoil = lpsutl.loadh5struct(rec_args.fname).ncoil;
seq_args = lpsutl.loadh5struct(rec_args.fname,'/seq_args'); % sequence arguments
ktraj_in = lpsutl.loadh5struct(rec_args.fname,'/ktraj').spoke_in; % spoke-in kspace trajectory
ktraj_out = lpsutl.loadh5struct(rec_args.fname,'/ktraj').spoke_out; % spoke-out kspace trajectory
if isfile(rec_args.fname_smaps)
    smaps = lpsutl.loadh5struct(rec_args.fname_smaps).real + ...
        1i*lpsutl.loadh5struct(rec_args.fname_smaps).imag; % sensitivity maps
else
    smaps = []; % leave empty if no sensitity map file provided
end

%% format the data and k-space trajectory
fprintf('formatting data and k-space trajectories...\n');
if isempty(rec_args.N)
    rec_args.N = seq_args.N_nom;
end
if isempty(rec_args.echoes2use)
    rec_args.echoes2use = 1:seq_args.nechoes;
end
if isempty(rec_args.ints2use)
    rec_args.ints2use = 1:seq_args.nint;
end
if isempty(rec_args.prjs2use)
    rec_args.prjs2use = 1:seq_args.nprj;
end
if isempty(rec_args.reps2use)
    rec_args.reps2use = 1:seq_args.nrep;
end
Ptotal = length(rec_args.ints2use)*length(rec_args.prjs2use)*length(rec_args.reps2use);
if isempty(rec_args.P)
    rec_args.P = length(rec_args.prjs2use)*length(rec_args.reps2use);
end

% select out desired indices
nechoes = length(rec_args.echoes2use);
ktraj_in = reshape(ktraj_in(:,rec_args.ints2use,rec_args.prjs2use,rec_args.reps2use,:),[],Ptotal,3);
ktraj_out = reshape(ktraj_out(:,rec_args.ints2use,rec_args.prjs2use,rec_args.reps2use,:),[],Ptotal,3);
kdata = reshape(kdata(:,rec_args.echoes2use,rec_args.ints2use,rec_args.prjs2use,rec_args.reps2use,:),[],nechoes,Ptotal,ncoil);

%% reformat into individual echoes/volumes
nvol = max(floor(Ptotal / rec_args.P),1);
Ptotal = nvol*rec_args.P;
ktraj_in = reshape(ktraj_in(:,1:Ptotal,:),[],nvol,3);
ktraj_out = reshape(ktraj_out(:,1:Ptotal,:),[],nvol,3);
kdata = reshape(kdata(:,:,1:Ptotal,:),[],nechoes,rec_args.P,nvol,ncoil);
kdata = reshape(permute(kdata,[1,3,2,4,5]),[],nechoes,nvol,ncoil);

%% handle sensitivity maps
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
    smaps = lpsutl.resample3d(smaps,rec_args.N*ones(1,3));
end

%% create the forward model and initialize the solution
fprintf('creating forward acquisition model...\n')
if rec_args.dcf_init
    [A,WA] = lpsutl.lps_model(ktraj_in, ktraj_out, smaps, rec_args, seq_args);
    fprintf('initializing solution with density compensated adjoint...\n')
    x0 = reshape(WA' * kdata, [rec_args.N*ones(1,3),nechoes,nvol]);
    x0 = norm(kdata(:),2) / norm(reshape(A*x0,[],1),2) * x0;
else % don't bother constructing WA
    A = lpsutl.lps_model(ktraj_in, ktraj_out, smaps, rec_args, seq_args);
    x0 = zeros([rec_args.N*ones(1,3),nechoes,nvol]);
end

%% solve the recon problem with CG
fprintf('solving recon problem with CG...\n')
C = fatrix2( ... % identity operator for tikhonov regularization
    'idim', rec_args.N*ones(1,3), ...
    'odim', rec_args.N*ones(1,3), ...
    'forw', @(~,x) sqrt(rec_args.beta) * x, ...
    'back', @(~,y) sqrt(rec_args.beta) * y ...
    );
if nechoes > 1 || nvol > 1
    C = kronI(nechoes,C);
    if nvol > 1
        C = kronI(nvol,C);
    end
end
x_star = qpwls_pcg1(x0,A,1,kdata(:),C,'niter',rec_args.niter);
x_star = reshape(x_star,size(x0));