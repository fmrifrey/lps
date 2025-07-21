function D = finite_diff(sz,dim,bc)
% function to return finite differencing operator over arbitrary dimension
% by David Frey
%
% inputs:
% sz - size of array to compute finite difference over
% dim - dimension(s) to compute finite difference over
% bc - boundary condition ('dirichlet' (default), 'neumann', or 'periodic')
%
% outputs:
% ob - finite differencing operator
%

    if nargin < 3 || isempty(bc)
        bc = 'dirichlet';
    end

    if length(dim) > 1
        Dblocks = cell(length(dim),1);
        for i = 1:length(dim) % loop through dimensions
            Dblocks{i} = lpsutl.finite_diff(sz,dim(i),bc);
        end
        % append operators as block column matrix
        D = lpsutl.block_col(Dblocks,0);
    else
        % set up fatrix operator
        D = fatrix2(...
            'idim', sz, ...
            'odim', sz, ...
            'forw', @(~,x) finite_diff_forw(x,dim,bc), ...
            'back', @(~,y) finite_diff_back(y,dim,bc) ...
            );
    end

end

function y = finite_diff_forw(x,dim,bc)

    % get dimensions
    nd = ndims(x);

    % permute the object so that desired dimension is last
    perm = [setdiff(1:nd, dim, 'stable'),dim];
    x_tmp = permute(x,perm);
    xv = reshape(x_tmp,[],size(x,dim)); % column vectorize

    % pad with boundary condition
    switch bc
        case 'dirichlet'
            yv = [
                xv(:,2:end) - xv(:,1:end-1), ...
                zeros(size(xv,1),1)
                ];
        case 'neumann'
            yv = [
                xv(:,2:end) - xv(:,1:end-1), ...
                xv(:,end) - xv(:,end-1)
                ];
        case 'periodic'
            yv = [
                xv(:,2:end) - xv(:,1:end-1), ...
                xv(:,1) - xv(:,end)
                ];
    end

    % reshape and permute back
    y_tmp = reshape(yv, size(x_tmp));
    [~,invperm] = sort(perm);
    y = permute(y_tmp,invperm);

end

function x = finite_diff_back(y,dim,bc)
    
    % get dimensions
    nd = ndims(y);

    % permute the object so that desired dimension is last
    perm = [setdiff(1:nd, dim, 'stable'),dim];
    y_tmp = permute(y,perm);
    yv = reshape(y_tmp,[],size(y,dim)); % column vectorize

    % take adjoint finite difference
    switch bc
        case 'dirichlet'
            xv = [
                -yv(:,1), ...
                yv(:,1:end-2) - yv(:,2:end-1), ...
                yv(:,end-1)
                ];
        case 'neumann'
            xv = [
                -yv(:,1), ...
                yv(:,1:end-3) - yv(:,2:end-2), ...
                yv(:,end-2) - yv(:,end-1) - yv(:,end), ...
                yv(:,end-1) + yv(:,end)
                ];
        case 'periodic'
            xv = [
                yv(:,end) - yv(:,1), ...
                yv(:,1:end-1) - yv(:,2:end)
                ];
    end

    % reshape and permute back
    x_tmp = reshape(xv, size(y_tmp));
    [~,invperm] = sort(perm);
    x = permute(x_tmp,invperm);

end