function s = loadh5struct(fname, groupName)
% loads in all data from h5 file as a structure
% by David Frey
%
% inputs:
% fname - name of h5 file to load in
% groupName - name of subgroup to load
%
% outputs:
% s - structure with same heirarchy as h5 file
%

    import lpsutl.*

    if nargin < 2
        groupName = '/';
    end

    info = h5info(fname, groupName);
    s = struct();

    % load datasets in this group
    for i = 1:length(info.Datasets)
        name = info.Datasets(i).Name;
        s.(name) = h5read(fname, [groupName '/' name]);
    end

    % recursively load subgroups
    for i = 1:length(info.Groups)
        fullGroupName = info.Groups(i).Name;
        [~, baseName] = fileparts(fullGroupName);
        groupField = matlab.lang.makeValidName(baseName);
        s.(groupField) = loadh5struct(fname, fullGroupName);
    end

end

