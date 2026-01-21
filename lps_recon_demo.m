%% Looping Star demo recon script
% by David Frey
%% set parameters
rec_args.fname = './raw_data.h5'; % input raw data .h5 file name (see lps_convert_data.m)
rec_args.fname_smaps = './smaps.h5'; % smaps input file name
rec_args.Q = 4; % number of compressed coils to use
rec_args.N = 64; % recon image size
rec_args.ints2use = []; % indices of interleaves to include (empty = all)
rec_args.prjs2use = []; % indices of projections to include (empty = all)
rec_args.reps2use = []; % indices of repetitions to include (empty = all)
rec_args.P = 32; % number of projections to use per frame (empty = nint*nprj)
rec_args.niter = 30; % number of CG iterations
rec_args.dcf_init = false; % option to initialize solution with density compensated NUFFT
rec_args.use_parfor = true; % option to use parfor loop in frame/coil-wise NUFFTs
rec_args.fermi_cutoff = 0.5; % fermi voxel basis function cutoff (frac of nominal resolution)
rec_args.fermi_rolloff = 0.2; % fermi voxel basis function rolloff (frac of nominal resolution)
rec_args.beta = 2^14; % quadratic roughness penalty regularization parameter
rec_args.debug = 0; % debug mode

%% load in data
fprintf('loading data...\n');
rec_args.fname = which(rec_args.fname);
kdata = lpsutl.loadh5struct(rec_args.fname,'/kdata').real + ...
    1i*lpsutl.loadh5struct(rec_args.fname,'/kdata').imag; % raw kspace data
seq_args = lpsutl.loadh5struct(rec_args.fname,'/seq_args'); % sequence arguments
ktraj_in = lpsutl.loadh5struct(rec_args.fname,'/ktraj').spoke_in; % spoke-in kspace trajectory
ktraj_out = lpsutl.loadh5struct(rec_args.fname,'/ktraj').spoke_out; % spoke-out kspace trajectory
if isfile(rec_args.fname_smaps)
    rec_args.fname_smaps = which(rec_args.fname_smaps);
    smaps = lpsutl.loadh5struct(rec_args.fname_smaps).real + ...
        1i*lpsutl.loadh5struct(rec_args.fname_smaps).imag; % sensitivity maps
else
    smaps = []; % leave empty if no sensitity map file provided
end

%% format the data and k-space trajectory
fprintf('formatting data and k-space trajectories...\n');
if isempty(rec_args.ints2use)
    rec_args.ints2use = 1:seq_args.nint;
end
if isempty(rec_args.prjs2use)
    rec_args.prjs2use = 1:seq_args.nprj;
end
if isempty(rec_args.reps2use)
    rec_args.reps2use = 1:seq_args.nrep;
end
if isempty(rec_args.P)
    rec_args.P = seq_args.nint*seq_args.nprj;
end

% select out desired indices
Ptotal = length(rec_args.ints2use)*length(rec_args.prjs2use)*length(rec_args.reps2use);
ktraj_in = reshape(ktraj_in(:,rec_args.ints2use,rec_args.prjs2use,rec_args.reps2use,:),[],Ptotal,3);
ktraj_out = reshape(ktraj_out(:,rec_args.ints2use,rec_args.prjs2use,rec_args.reps2use,:),[],Ptotal,3);
kdata = reshape(kdata(:,rec_args.ints2use,rec_args.prjs2use,rec_args.reps2use,:),[],Ptotal,size(kdata,5));

% reformat into individual frames
nvol = floor(Ptotal / rec_args.P);
Ptotal = nvol*rec_args.P;
ktraj_in = reshape(ktraj_in(:,1:Ptotal,:),[],nvol,3);
ktraj_out = reshape(ktraj_out(:,1:Ptotal,:),[],nvol,3);
kdata = reshape(kdata(:,1:Ptotal,:),[],nvol,size(kdata,3));

%% handle sensitivity maps
if isempty(smaps)
    fprintf('no sensitivity maps provided -> compressing data...\n');
    rec_args.Q = 1;
    if size(kdata,3) > 1
        kdata = ir_mri_coil_compress(kdata,'ncoil',rec_args.Q);
    end
    smaps = ones([rec_args.N*ones(1,3),1]);
else
    fprintf('compressing data to %d virtual coils...\n', rec_args.Q);
    if isempty(rec_args.Q)
        rec_args.Q = size(kdata,3); % default - use all coils
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
    x0 = WA' * kdata;
else % don't bother constructing WA
    A = lpsutl.lps_model(ktraj_in, ktraj_out, smaps, rec_args, seq_args);
    x0 = zeros([rec_args.N*ones(1,3),nvol]);
end

%% solve the recon problem with CG
fprintf('solving recon problem with CG...\n')
C = kronI(nvol,Reg1(ones(rec_args.N*ones(1,3)),'beta',rec_args.beta).C1);
x_star = qpwls_pcg1(x0,A,1,kdata(:),C,'niter',rec_args.niter);
x_star = reshape(x_star,size(x0));
