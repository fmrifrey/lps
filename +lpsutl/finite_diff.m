function D = finite_diff(sz,dim,bc)
% function to return finite differencing operator over arbitrary dimension
% by David Frey
%
% inputs:
% sz - size of array to compute finite difference over
% dim - dimension to compute finite difference over
% bc - boundary condition ('dirichlet' (default), 'neumann', or 'periodic')
%
% outputs:
% ob - finite differencing operator
%

    if nargin < 3 || isempty(bc)
        bc = 'dirichlet';
    end

    % set up fatrix operator
    D = fatrix2(...
        'idim', sz, ...
        'odim', sz, ...
        'forw', @(~,x) finite_diff_forw(x,dim,bc), ...
        'back', @(~,y) finite_diff_back(y,dim,bc) ...
        );

end

function y = finite_diff_forw(x,dim,bc)

    % get dimensions
    nd = ndims(x);

    % permute the object so that desired dimension is last
    x_tmp = permute(x,[setdiff(1:nd, dim, 'stable'),dim]);
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
    y_tmp = reshape(yv, size(x));
    y = permute(y_tmp,[setdiff(1:nd, dim, 'stable'),dim]);

end

function x = finite_diff_back(y,dim,bc)
    
    % get dimensions
    nd = ndims(y);

    % permute the object so that desired dimension is last
    y_tmp = permute(y,[setdiff(1:nd, dim, 'stable'),dim]);
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
    x_tmp = reshape(xv, size(y));
    x = permute(x_tmp,[setdiff(1:nd, dim, 'stable'),dim]);

end