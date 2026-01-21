function out = fftc(X,dim)
% function to compute fft along any dim of a tensor
% by David Frey
%
% inputs:
% X - input tensor
% dim - dimensions to FT over
%
% outputs:
% out - fourier data
%

    % set default dimensions
    if nargin<2 || isempty(dim)
        dim = 1:ndims(X);
    elseif any(dim(:) > ndims(X)) || any(dim(:) < 1)
        error('dimensions out of range');
    end
    
    % define fourier transform with scaling and shifts
    fftc1d = @(x,d) 1/sqrt(size(x,d))*fftshift(fft(fftshift(x,d),[],d),d);
    
    % fourier transform along each requested dimension
    out = X;
    for d = dim
        out = fftc1d(out, d);
    end
    
end