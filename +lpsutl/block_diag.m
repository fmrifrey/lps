function ob = block_diag(Ablocks,usepar)
% construct block diagonal 'fatrix' operator with parallelization option
% by David Frey
%
% inputs:
% Ablocks - cell array of block operators
% usepar - option to parallelize calculation over blocks
%
% outputs:
% ob - block diagonal operator
%

    % set default parallelization
    if nargin < 2 || isempty(usepar)
        usepar = 0;
    end

    % get number of blocks
    nblocks = length(Ablocks);

    % create the fatrix object
    ob = fatrix2('idim', [Ablocks{1}.idim,nblocks], ...
        'odim', [Ablocks{1}.odim,nblocks], ...
        'forw', @(~,x)block_diag_forw(Ablocks,x,usepar), ...
        'back', @(~,y)block_diag_back(Ablocks,y,usepar) ...
        );

end

function y = block_diag_forw(Ablocks, x, usepar)

    % get number of blocks and set up cell arrays
    nblocks = length(Ablocks);
    xtmp= reshape(x,prod(Ablocks{1}.idim),nblocks);
    xblocks = cell(nblocks,1);
    for i = 1:nblocks
        xblocks{i} = reshape(xtmp(:,i),Ablocks{1}.idim);
    end
    yblocks = cell(nblocks,1);

    % loop through blocks and compute adjoint operation
    if usepar && nblocks > 1
        parfor i = 1:nblocks
            yblocks{i} = Ablocks{i} * xblocks{i};
        end
    else
        for i = 1:nblocks
            yblocks{i} = Ablocks{i} * xblocks{i};
        end
    end

    % concatenate block operations
    y = cat(length(Ablocks{1}.odim)+1,yblocks{:});

end

function x = block_diag_back(Ablocks, y, usepar)

    % get number of blocks and set up cell arrays
    nblocks = length(Ablocks);
    ytmp = reshape(y,prod(Ablocks{1}.odim),nblocks);
    yblocks = cell(nblocks,1);
    for i = 1:nblocks
        yblocks{i} = reshape(ytmp(:,i),Ablocks{1}.odim);
    end
    xblocks = cell(nblocks,1);

    % loop through blocks and compute adjoint operation
    if usepar && nblocks > 1
        parfor i = 1:nblocks
            xblocks{i} = Ablocks{i}' * yblocks{i};
        end
    else
        for i = 1:nblocks
            xblocks{i} = Ablocks{i}' * yblocks{i};
        end
    end

    % sum together block operations
    x = cat(length(Ablocks{1}.idim)+1,xblocks{:});

end