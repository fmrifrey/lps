function [xk,cost] = cg(A, x0, b, varargin)

    % set defaults
    arg.niter = 30;
    arg.tol = 1e-3;
    arg.chat = 0;

    % parse input arguments
    arg = vararg_pair(arg,varargin);

    % initialize variables
    r0 = b - A*x0;
    if size(A, 1) ~= size(A, 2) % if non-square, solve normal equations
        r0 = A' * r0;
    end
    rk = r0;
    pk = r0;
    xk = x0;
    rktrk = real(rk(:)'*rk(:));

    if nargout > 1 % calculate cost
        cost = zeros(arg.niter+1,1);
        cost(1) = norm(r0(:),2);
        if arg.chat
            fprintf('cg(): initial cost = %f\n', cost(1));
        end
    end

    % loop through iterations
    for k = 1:arg.niter
        fprintf('cg(): beginning iteration %d...\n', k);

        % check for convergence
        if rktrk < arg.tol
            cost = cost(1:k);
            if arg.chat
                fprintf('\tconvergence criteria met --> terminating early...\n')
            end
            return
        end

        % start timer
        timerk = tic;

        % calculate step size
        Apk = A*pk;
        if size(A, 1) ~= size(A, 2) % if non-square, solve normal equations
            Apk = A' * Apk;
        end
        pktApk = Apk(:)' * pk(:);
        alpha = rktrk / pktApk;

        % gradient step
        xk1 = xk + alpha * pk;
        rk1 = rk - alpha * Apk;

        % update aux variables
        rk1trk1 = real(rk1(:)'*rk1(:));
        beta = rk1trk1 / rktrk;
        pk1 = beta * pk + rk1;

        % update next iteration
        rk = rk1;
        xk = xk1;
        pk = pk1;
        rktrk = rk1trk1;

        if arg.chat % print elapsed time
            fprintf('\titeration time = %fs\n', toc(timerk));
        end

        if nargout > 1 % calculate cost
            cost(k+1) = norm(rk(:),2);
            if arg.chat
                fprintf('\tcost = %f\n', cost(1));
            end
        end

    end

end