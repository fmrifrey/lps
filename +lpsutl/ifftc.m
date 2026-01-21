function out = ifftc(X,dim)
% function to compute ifft along any dim of a tensor
% by David Frey
%
% inputs:
% X - input tensor
% dim - dimensions to IFT over
%
% outputs:
% out - inverse fourier data
%

    % set default dimensions
    if nargin<2 || isempty(dim)
        dim = 1:ndims(X);
    elseif any(dim(:) > ndims(X)) || any(dim(:) < 1)
        error('dimensions out of range');
    end
    
    % define fourier transform with scaling and shifts
    ifftc1d = @(x,d) sqrt(size(x,d))*ifftshift(ifft(ifftshift(x,d),[],d),d);
    
    % fourier transform along each requested dimension
    out = X;
    for d = dim
        out = ifftc1d(out, d);
    end
    
end