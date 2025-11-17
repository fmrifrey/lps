function A = A_volwise(Fs_in,Fs_out,Ws_in,Ws_out,smaps,usepar)
% creates a volume-wise foward operator for looping star data
% forward operator A = [[A{1}    0    ...    0      ]  * kronI(nvol,S)
%                       [0       A{2} ...    0      ]
%                       ...
%                       [0       ...         A{nvol}]]
% where A{i} = Ws_in{i}*Fs_in{i} + Ws_out{i}*Fs_out{i}, i = 1, ... nvol
% by David Frey
%
% inputs:
% Fs_in - cell array of volume-wise nufft operators for LpS spoke-in data
% Fs_out - cell array of volume-wise nufft operators for LpS spoke-out data
% Ws_in - cell array of weighting matrices (i.e. filters) for LpS spoke-in
%   data (use 1 for no weighting)
% Ws_out - cell array of weighting matrices (i.e. filters) for LpS spoke-out
%   data (use 1 for no weighting)
% smaps - sensitivity maps
% usepar - option to use parfor loop over volumes
%
% outputs:
% A - full constructed system operator
%

    idimv = Fs_in{1}.idim; % single vol input dimension
    odimv = Fs_in{1}.odim; % single vol output dimension
    nvol = length(Fs_in); % number of volumes

    % set defaults
    if nargin < 5 || isempty(smaps)
        smaps = ones(idimv);
    end
    nc = size(smaps,4); % number of coils

    if nargin < 6 || isempty(usepar)
        usepar = false;
    end

    % if weighting is a scalar, make it a diagonal
    if isscalar(Ws_in) && ~iscell(Ws_in)
        Ws_in = repmat({Gdiag(Ws_in*ones(odimv(:),1))},[nvol,1]);
    end
    if isscalar(Ws_out) && ~iscell(Ws_in)
        Ws_out = repmat({Gdiag(Ws_out*ones(odimv(:),1))},[nvol,1]);
    end
    
    % define forward operator
    function b = A_fwd(x)
        b = zeros([odimv,nc,nvol]);
        if usepar
            parfor v = 1:nvol
                Av = senseop(Ws_in{v}*Fs_in{v} + Ws_out{v}*Fs_out{v}, smaps);
                b(:,:,v) = Av * x(:,:,:,v);
            end
        else
            for v = 1:nvol
                Av = senseop(Ws_in{v}*Fs_in{v} + Ws_out{v}*Fs_out{v}, smaps);
                b(:,:,v) = Av * x(:,:,:,v);
            end
        end
    end

    % define adjoint operator
    function x = A_adj(b)
        x = zeros([idimv,nvol]);
        if usepar
            parfor v = 1:nvol
                Av = senseop(Ws_in{v}*Fs_in{v} + Ws_out{v}*Fs_out{v}, smaps);
                x(:,:,:,v) = Av' * b(:,:,v);
            end
        else
            for v = 1:nvol
                Av = senseop(Ws_in{v}*Fs_in{v} + Ws_out{v}*Fs_out{v}, smaps);
                x(:,:,:,v) = Av' * b(:,:,v);
            end
        end
    end
    
    % create fatrix object
    A = fatrix2('idim', [idimv,nvol], ...
        'odim', [odimv,nc,nvol], ...
        'forw', @(~,x) A_fwd(x), ...
        'back', @(~,b) A_adj(b));

end

function AS = senseop(A,smaps)
% adds sensitivity encoding to forward operator A
% by David Frey
%
% inputs:
% A - forward fatrix operator
% smaps - sensitivity maps
%
% outputs:
% AS - forward fatrix operator including sensitivity encoding
%

    % get number of coils
    nc = size(smaps,length(A.idim)+1);

    % get image mask
    msk = A.imask;

    % mask out the smaps
    smaps = reshape(smaps,[],nc);
    smaps_msk = smaps(msk(:),:);

    % create cell array of sensitivity operators
    AS_cell = cell(nc,1);
    for i = 1:nc
        AS_cell{i} = A*Gdiag(smaps_msk(:,i), 'mask', msk);
    end

    % concatenate into block column matrix
    AS = cat(1,AS_cell{:});

end