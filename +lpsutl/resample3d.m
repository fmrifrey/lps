function data_new = resample3d(data,sz_new)
% resamples the ND data along the first 3 dimensions
% by David Frey
%
% inputs:
% data - ND data matrix (at least 3D)
% sz_new - new 3D size
%
% outputs:
% data_new - resampled data
%

% get old size
sz = size(data);

% create sampling grids
[X, Y, Z] = ndgrid(linspace(-1, 1, sz(1)), ...
    linspace(-1, 1, sz(2)), ...
    linspace(-1, 1, sz(3)));
[Xq, Yq, Zq] = ndgrid(linspace(-1, 1, sz_new(1)), ...
    linspace(-1, 1, sz_new(2)), ...
    linspace(-1, 1, sz_new(3)));

% interpolate the data
if length(sz) > 3
    data_new = zeros([sz_new(1:3),sz(4:end)]);
    for i = 1:prod(sz(4:end))
        data_new(:,:,:,i) = interp3(Y, X, Z, data(:,:,:,i), Yq, Xq, Zq, 'cubic', 0);
    end
else
    data_new = interp3(Y, X, Z, data, Yq, Xq, Zq, 'cubic', 0);
end

end

