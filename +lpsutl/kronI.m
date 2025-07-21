function ob = kronI(M, A, usepar)
% create repeated object over M times (i.e. I_M \otimes A)
% by David Frey
%
% inputs:
% M - size of identity matrix (number of blocks)
% A - block operator
% usepar - option to parallelize calculation over blocks
%
% outputs:
% ob - block diagonal operator
%

    % set default parallelization
    if nargin < 2 || isempty(usepar)
        usepar = 0;
    end

    % create the block diagonal operator
    Ablocks = repmat({A},M,1);
    ob = lpsutl.block_diag(Ablocks,usepar);

end