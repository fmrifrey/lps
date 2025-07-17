function idcs = snake_caipi_idcs(N, Ry, Rz, delta, Nacs)
% generate (row, col) indices in a snake order from the top left to bottom
% with Nacs fully sampled lines at center of kspace, and (Ry,Rz) CAIPI
% outside of ACS region
% by David Frey
%
% inputs:
% N - matrix size
% Ry - y acceleration (outside of ACS region)
% Rz - z acceleration (outside of ACS region)
% delta - CAIPI shift per line
% Nacs - number of ACS lines
%
% outputs:
% idcs - ordered y and z indicies
%

    % define center index
    center = floor(N/2);

    % initialize to center of grid
    idcs = [];

    % get CAIPI indicies
    CAIPI_idcs = [];
    for iy = 1:Ry:N
        z_shift = delta*floor((iy-1)/Ry);
        for iz = mod((0:Rz:N-1)  + z_shift, N) + 1
            CAIPI_idcs = [CAIPI_idcs; iy, iz];
        end
    end

    % intiialize z order
    iz_ordered = 1:N;
    flip_next = false;

    % loop through y locations
    for iy = 1:N

        % flip the z ordering
        if flip_next
            iz_ordered = flip(iz_ordered);
            flip_next = false;
        end

        for iz = iz_ordered
            if 2*abs(iy - center) <= Nacs && 2*abs(iz - center) <= Nacs
                % ACS index
                idcs = [idcs; iy, iz];
                flip_next = true;
            elseif ismember([iy,iz], CAIPI_idcs,'rows')
                % CAIPI index
                idcs = [idcs; iy, iz];
                flip_next = true;
            end
        end
        
    end

    % make sure no indicies are repeated
    idcs = unique(idcs,'rows','stable');

end
