function [A,WA] = lps_model(ktraj_in, ktraj_out, ...
    smaps, rec_args, seq_args)
% creates forward acquisition model for looping star reconstruction
% by David Frey (djfrey@umich.edu)
%
% inputs:
% ktraj_in - spoke-in kspace trajectory (n x nvol x 3)
% ktraj_out - spoke-out kspace trajectory (n x nvol x 3)
% smaps - sensitivity maps (N x N x N x Q)
% rec_args - reconstruction arguments structure created by lps_recon_demo.m
% seq_args - sequence arguments structure created by lps_write_seq.m
%
% outputs:
% A - forward acquisition model
% WA - density-weighted acquisition model (only to be used for
%   initialization)
%

    % get sizes
    M = size(ktraj_in,1);
    necho = length(rec_args.echos2use);
    nvol = size(ktraj_in,2);
    N = size(smaps,1);
    Q = size(smaps,4);

    % create spherical mask for k-space gridding
    [Xgrd,Ygrd,Zgrd] = ndgrid(linspace(-1,1,N), ...
        linspace(-1,1,N), ...
        linspace(-1,1,N));
    msk = (Xgrd.^2 + Ygrd.^2 + Zgrd.^2) <= 1;

    % create k-space masks
    kmsk_in = 2*seq_args.fov/rec_args.N * vecnorm(reshape(ktraj_in(:,1,:),[],3),2,2) ...
        <= min(1, seq_args.N_nom/rec_args.N * (2*rec_args.fermi_cutoff + rec_args.fermi_rolloff));
    kmsk_out = 2*seq_args.fov/rec_args.N * vecnorm(reshape(ktraj_out(:,1,:),[],3),2,2) ...
        <= min(1, seq_args.N_nom/rec_args.N * (2*rec_args.fermi_cutoff + rec_args.fermi_rolloff));

    % set up volume-wise nufft operators
    Fs_in = cell(nvol,1);
    Fs_out = cell(nvol,1);
    w_in = zeros(sum(1*kmsk_in),nvol);
    w_out = zeros(sum(1*kmsk_out),nvol);
    nufft_args = {rec_args.N*ones(1,3), 6*ones(1,3), 2*rec_args.N*ones(1,3), ...
        rec_args.N/2*ones(1,3), 'table', 2^10, 'minmax:kb'};
    for ivol = 1:nvol % loop through volumes and assemble nufft object for current vol
        omegav_in = 2*pi*seq_args.fov/rec_args.N * ...
            reshape(ktraj_in(kmsk_in,ivol,:),[],3);
        omegav_out = 2*pi*seq_args.fov/rec_args.N * ...
            reshape(ktraj_out(kmsk_out,ivol,:),[],3);
        Fs_in{ivol} = Gnufft(msk, [omegav_in, nufft_args]);
        Fs_out{ivol} = Gnufft(msk, [omegav_out, nufft_args]);
        if nargout > 1 % compute density compensation
            w_in(:,ivol) = lpsutl.dcf_pipe(Fs_in{ivol});
            w_out(:,ivol) = lpsutl.dcf_pipe(Fs_out{ivol});
        end
    end

    % create fermi windowing functions for voxel basis function
    ferm_fun = @(x,fwhm,tband) 1 ./ (1 + exp(8.788 * (abs(x) - fwhm/2)/tband)); % fermi windowing function for basis function
    ferm_omega = @(omega) ferm_fun(vecnorm(omega,2,2)/pi, ...
        seq_args.N_nom/rec_args.N * 2*rec_args.fermi_cutoff, ...
        seq_args.N_nom/rec_args.N * rec_args.fermi_rolloff);

    % create forward operation
    function y = A_fwd(x,use_dcf)
        t = tic;
        y = zeros(M,necho*nvol*Q);
        if rec_args.use_parfor
            parfor i = 1:necho*Q*nvol
                [e,v,q] = ind2sub([necho,nvol,Q],i);
                y(:,i) = Ai_fwd(x(:,:,:,e,v), smaps(:,:,:,q), ...
                    Fs_in{v}, Fs_out{v}, w_in(:,v), w_out(:,v), ...
                    kmsk_in, kmsk_out, ferm_omega, use_dcf);
            end
        else
            for i = 1:necho*Q*nvol
                [e,v,q] = ind2sub([necho,nvol,Q],i);
                y(:,i) = Ai_fwd(x(:,:,:,e,v), smaps(:,:,:,q), ...
                    Fs_in{v}, Fs_out{v}, w_in(:,v), w_out(:,v), ...
                    kmsk_in, kmsk_out, ferm_omega, use_dcf);
            end
        end
        y = reshape(y,M,necho,nvol,Q);
        if rec_args.debug > 0
            fprintf('forward op time = %.3fs\n', toc(t));
        end
    end

    % create adjoint operation
    function x = A_adj(y,use_dcf)
        t = tic;
        x_evq = zeros([N*ones(1,3),necho*nvol*Q]);
        if rec_args.use_parfor
            parfor i = 1:necho*Q*nvol
                [e,v,q] = ind2sub([necho,nvol,Q],i);
                x_evq(:,:,:,i) = Ai_adj(y(:,e,v,q), smaps(:,:,:,q), ...
                    Fs_in{v}, Fs_out{v}, w_in(:,v), w_out(:,v), ...
                    kmsk_in, kmsk_out, ferm_omega, use_dcf);
            end
        else
            for i = 1:necho*Q*nvol
                [e,v,q] = ind2sub([necho,nvol,Q],i);
                x_evq(:,:,:,i) = Ai_adj(y(:,e,v,q), smaps(:,:,:,q), ...
                    Fs_in{v}, Fs_out{v}, w_in(:,v), w_out(:,v), ...
                    kmsk_in, kmsk_out, ferm_omega, use_dcf);
            end
        end
        x_evq = reshape(x_evq,[N*ones(1,3),necho,nvol,Q]);
        x = sum(x_evq,6); % sum over coil dimension
        if rec_args.debug > 0
            fprintf('adjoint op time = %.3fs\n', toc(t));
        end
    end

    % create the forward model fatrix operator
    A = fatrix2( ...
        'idim', [N*ones(1,3),necho,nvol], ...
        'odim', [M,necho,nvol,Q], ...
        'forw', @(~,x) A_fwd(x,false), ...
        'back', @(~,y) A_adj(y,false) ...
        );

    % create the forward model fatrix operator
    WA = fatrix2( ...
        'idim', [N*ones(1,3),necho,nvol], ...
        'odim', [M,necho,nvol,Q], ...
        'forw', @(~,x) A_fwd(x,true), ...
        'back', @(~,y) A_adj(y,true) ...
        );

end

% create foward operation for a SINGLE VOLUME/COIL/ECHO
function y = Ai_fwd(x, smap, F_in, F_out, w_in, w_out, ...
    kmsk_in, kmsk_out, ferm_omega, use_dcf)

    % multiply by sense map
    x = smap .* x;
    
    % forward NUFFT
    y_in = F_in * x;
    y_out = F_out * x;

    % weigh by fermi filter
    y_in = ferm_omega(F_in.arg.arg{1}) .* y_in;
    y_out = ferm_omega(F_out.arg.arg{1}) .* y_out;

    % weigh by density compensation
    if use_dcf
        y_in = w_in .* y_in;
        y_out = w_out .* y_out;
    end

    % embed into k-space mask
    y = embed(y_in, kmsk_in) + embed(y_out, kmsk_out);

end

% create adjoint operation for a SINGLE VOLUME/COIL/ECHO
function x = Ai_adj(y, smap, F_in, F_out, w_in, w_out, ...
    kmsk_in, kmsk_out, ferm_omega, use_dcf)

    % get msk-indexed data
    y_in = y(kmsk_in);
    y_out = y(kmsk_out);

    % weigh by dcf
    if use_dcf
        y_in = w_in .* y_in;
        y_out = w_out .* y_out;
    end

    % weigh by fermi function
    y_in = ferm_omega(F_in.arg.arg{1}) .* y_in;
    y_out = ferm_omega(F_out.arg.arg{1}) .* y_out;

    % adjoint NUFFT
    x = F_in' * y_in + F_out' * y_out;
    
    % embed into image mask
    x = embed(x, F_in.arg.mask);
    
    % multiply by conjugated coil map
    x = conj(smap) .* x;

end