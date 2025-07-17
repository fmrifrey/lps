function out = fftc(X,N,dim)
% function to compute fft along any dim of a tensor
% by David Frey
%
% inputs:
% X - input tensor
% N - number of points in FT
% dim - dimensions to FT over
%
% outputs:
% out - fourier data
%

    % set default dimensions
    if nargin<3 || isempty(dim)
        dim = 1:ndims(X);
    elseif any(dim(:) > ndims(X)) || any(dim(:) < 1)
        error('dimensions out of range');
    end

    % set default N
    if nargin<2 || isempty(N)
        N = size(X,dim);
    elseif isscalar(N)
        N = N*ones(1,length(dim));
    end
    
    % define fourier transform with scaling and shifts
    fftc1d = @(x,n,d) 1/sqrt(size(x,d))*fftshift(fft(fftshift(x,d),n,d),d);
    
    % fourier transform along each requested dimension
    out = X;
    for i = 1:length(dim)
        out = fftc1d(out, N(i), dim(i));
    end
    
end