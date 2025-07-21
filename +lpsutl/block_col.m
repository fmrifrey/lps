function ob = block_col(Ablocks,usepar)
% construct block column 'fatrix' operator with parallelization option
% by David Frey
%
% inputs:
% Ablocks - cell array of block operators
% usepar - option to parallelize calculation over blocks
%
% outputs:
% ob - block column operator
%

    % set default parallelization
    if nargin < 2 || isempty(usepar)
        usepar = 0;
    end

    % get number of blocks
    nblocks = length(Ablocks);

    % create masks
    imask = Ablocks{1}.imask;
    omask = [];
    for i = 1:nblocks
        omask = cat(length(Ablocks{1}.odim)+1, omask, Ablocks{i}.omask);
    end
    omask = 1*omask > 0;

    % create the fatrix object
    ob = fatrix2( ...
        'idim', Ablocks{1}.idim, ...
        'odim', [Ablocks{1}.odim,nblocks], ...
        'imask', imask, ...
        'omask', omask, ...
        'forw', @(~,x)block_col_forw(Ablocks,x,usepar), ...
        'back', @(~,y)block_col_back(Ablocks,y,usepar) ...
        );

end

function y = block_col_forw(Ablocks, x, usepar)

    % get number of blocks and set up cell arrays
    nblocks = length(Ablocks);
    yblocks = cell(nblocks,1);

    % loop through blocks and compute adjoint operation
    if usepar && nblocks > 1
        parfor i = 1:nblocks
            yblocks{i} = Ablocks{i} * x;
            if isscalar(Ablocks{i}.idim) % handle 1D input case
                yblocks{i} = fatrix2_embed(Ablocks{i}.omask, ...
                    Ablocks{i}.odim, yblocks{i});
            end
        end
    else
        for i = 1:nblocks
            yblocks{i} = Ablocks{i} * x;
            if isscalar(Ablocks{i}.idim) % handle 1D input case
                yblocks{i} = fatrix2_embed(Ablocks{i}.omask, ...
                    Ablocks{i}.odim, yblocks{i});
            end
        end
    end

    % concatenate block operations
    y = cat(length(Ablocks{1}.odim)+1,yblocks{:});

end

function x = block_col_back(Ablocks, y, usepar)

    % get number of blocks and set up cell arrays
    nblocks = length(Ablocks);
    ytmp = reshape(y,prod(Ablocks{1}.odim),nblocks);
    yblocks = cell(nblocks,1);
    for i = 1:nblocks
        yblocks{i} = reshape(ytmp(:,i),[Ablocks{1}.odim,1]);
    end
    xblocks = cell(nblocks,1);

    % loop through blocks and compute adjoint operation
    if usepar && nblocks > 1
        parfor i = 1:nblocks
            xblocks{i} = Ablocks{i}' * yblocks{i};
            if isscalar(Ablocks{i}.odim) % handle 1D output case
                xblocks{i} = fatrix2_embed(Ablocks{i}.imask, ...
                    Ablocks{i}.idim, xblocks{i});
            end
        end
    else
        for i = 1:nblocks
            xblocks{i} = Ablocks{i}' * yblocks{i};
            if isscalar(Ablocks{i}.odim) % handle 1D output case
                xblocks{i} = fatrix2_embed(Ablocks{i}.imask, ...
                    Ablocks{i}.idim, xblocks{i});
            end
        end
    end

    % sum together block operations
    x = zeros(Ablocks{1}.idim);
    for i = 1:nblocks
        x = x + xblocks{i};
    end

end