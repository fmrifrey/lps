function saveh5struct(fname,s,groupName)
% save heirarchical structure data an h5 file or subgroup
% by David Frey
%
% inputs:
% fname - name of h5 file to load in
% groupName - name of subgroup to load
% s - structure with same heirarchy as h5 file
%

    import lpsutl.*

    if nargin < 3 || isempty(groupName)
        groupName = '/';
    end

    % save datasets in this group
    if isstruct(s)
        fn = fieldnames(s);
        for i = 1:length(fn)
            % recursively save subgroups
            saveh5struct(fname,s.(fn{i}),[groupName '/' fn{i}]);
        end
    else
        switch class(s)
            case 'char'
                s = string(s);
            case 'logical'
                s = 1*s;
        end
        h5create(fname, groupName, size(s), 'Datatype', class(s));
        h5write(fname, groupName, s);
    end

end

