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
    nvol = size(ktraj_in,2);
    N = size(smaps,1);
    Q = size(smaps,4);

    % create spherical mask for k-space gridding
    [Xgrd,Ygrd,Zgrd] = ndgrid(linspace(-1,1,N), ...
        linspace(-1,1,N), ...
        linspace(-1,1,N));
    msk = (Xgrd.^2 + Ygrd.^2 + Zgrd.^2) <= 1;

    % set up volume-wise nufft operators
    Fs_in = cell(nvol,1);
    Fs_out = cell(nvol,1);
    w_in = zeros(M,nvol);
    w_out = zeros(M,nvol);
    nufft_args = {rec_args.N*ones(1,3), 6*ones(1,3), 2*rec_args.N*ones(1,3), ...
        rec_args.N/2*ones(1,3), 'table', 2^10, 'minmax:kb'};
    for ivol = 1:nvol % loop through volumes and assemble nufft object for current vol
        omegav_in = 2*pi*seq_args.fov/rec_args.N * ...
            reshape(ktraj_in(:,ivol,:),[],3);
        omegav_out = 2*pi*seq_args.fov/rec_args.N * ...
            reshape(ktraj_out(:,ivol,:),[],3);
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
        y = zeros(M,nvol*Q);
        if rec_args.use_parfor
            parfor i = 1:Q*nvol
                [v,q] = ind2sub([nvol,Q],i);
                Sx_vq = smaps(:,:,:,q) .* x(:,:,:,v);
                if use_dcf % density compensated
                    y(:,i) = w_in(:,v) .* ferm_omega(Fs_in{v}.arg.arg{1}) .* (Fs_in{v}*Sx_vq) + ...
                        w_out(:,v) .* ferm_omega(Fs_out{v}.arg.arg{1}) .* (Fs_out{v}*Sx_vq);
                else
                    y(:,i) = ferm_omega(Fs_in{v}.arg.arg{1}) .* (Fs_in{v}*Sx_vq) + ...
                        ferm_omega(Fs_out{v}.arg.arg{1}) .* (Fs_out{v}*Sx_vq);
                end
            end
        else
            for i = 1:Q*nvol
                [v,q] = ind2sub([nvol,Q],i);
                Sx_vq = smaps(:,:,:,q) .* x(:,:,:,v);
                if use_dcf % density compensated
                    y(:,i) = w_in(:,v) .* ferm_omega(Fs_in{v}.arg.arg{1}) .* (Fs_in{v}*Sx_vq) + ...
                        w_out(:,v) .* ferm_omega(Fs_out{v}.arg.arg{1}) .* (Fs_out{v}*Sx_vq);
                else
                    y(:,i) = ferm_omega(Fs_in{v}.arg.arg{1}) .* (Fs_in{v}*Sx_vq) + ...
                        ferm_omega(Fs_out{v}.arg.arg{1}) .* (Fs_out{v}*Sx_vq);
                end
            end
        end
        y = reshape(y,M,nvol,Q);
        if rec_args.debug > 0
            fprintf('forward op time = %.3fs\n', toc(t));
        end
    end

    % create adjoint operation
    function x = A_adj(y,use_dcf)
        t = tic;
        x_vq = zeros([N*ones(1,3),nvol*Q]);
        if rec_args.use_parfor
            parfor i = 1:Q*nvol
                [v,q] = ind2sub([nvol,Q],i);
                if use_dcf
                    x_vq(:,:,:,i) = embed( ...
                        Fs_in{v}'*(ferm_omega(Fs_in{v}.arg.arg{1}) .* (w_in(:,v) .* y(:,v,q))) + ...
                        Fs_out{v}'*(ferm_omega(Fs_out{v}.arg.arg{1}) .* (w_out(:,v) .* y(:,v,q))), ...
                        msk);
                else
                    x_vq(:,:,:,i) = embed( ...
                        Fs_in{v}'*(ferm_omega(Fs_in{v}.arg.arg{1}) .* y(:,v,q)) + ...
                        Fs_out{v}'*(ferm_omega(Fs_out{v}.arg.arg{1}) .* y(:,v,q)), ...
                        msk);
                end
                x_vq(:,:,:,i) = conj(smaps(:,:,:,q)) .* x_vq(:,:,:,i);
            end
        else
            for i = 1:Q*nvol
                [v,q] = ind2sub([nvol,Q],i);
                if use_dcf
                    x_vq(:,:,:,i) = embed( ...
                        Fs_in{v}'*(ferm_omega(Fs_in{v}.arg.arg{1}) .* (w_in(:,v) .* y(:,v,q))) + ...
                        Fs_out{v}'*(ferm_omega(Fs_out{v}.arg.arg{1}) .* (w_out(:,v) .* y(:,v,q))), ...
                        msk);
                else
                    x_vq(:,:,:,i) = embed( ...
                        Fs_in{v}'*(ferm_omega(Fs_in{v}.arg.arg{1}) .* y(:,v,q)) + ...
                        Fs_out{v}'*(ferm_omega(Fs_out{v}.arg.arg{1}) .* y(:,v,q)), ...
                        msk);
                end
                x_vq(:,:,:,i) = conj(smaps(:,:,:,q)) .* x_vq(:,:,:,i);
            end
        end
        x_vq = reshape(x_vq,[N*ones(1,3),nvol,Q]);
        x = sum(x_vq,5); % sum over coil dimension
        if rec_args.debug > 0
            fprintf('adjoint op time = %.3fs\n', toc(t));
        end
    end

    % create the forward model fatrix operator
    A = fatrix2( ...
        'idim', [N*ones(1,3),nvol], ...
        'odim', [M,nvol,Q], ...
        'forw', @(~,x) A_fwd(x,false), ...
        'back', @(~,y) A_adj(y,false) ...
        );

    % create the forward model fatrix operator
    WA = fatrix2( ...
        'idim', [N*ones(1,3),nvol], ...
        'odim', [M,nvol,Q], ...
        'forw', @(~,x) A_fwd(x,true), ...
        'back', @(~,y) A_adj(y,true) ...
        );

end