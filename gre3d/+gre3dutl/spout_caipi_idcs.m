function idcs = spout_caipi_idcs(N, Ry, Rz, delta, Nacs)
% generate (row, col) indices in a spiral-out order from the center
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
    ny = center;
    nz = center;
    idcs = [ny, nz];

    % get CAIPI indicies
    CAIPI_idcs = [];
    for iy = 1:Ry:N
        z_shift = delta*floor((iy-1)/Ry);
        for iz = mod((0:Rz:N-1)  + z_shift, N) + 1
            CAIPI_idcs = [CAIPI_idcs; iy, iz];
        end
    end
    
    % define movement directions (right, down, left, up)
    directions = [0 1; 1 0; 0 -1; -1 0];

    % initialize
    step_size = 1; % initial step size
    breakloop = false;
    while true

        for d = 1:4  % iterate over directions
            actual_steps = step_size; % store current step size

            for i = 1:actual_steps
                ny = ny + directions(d, 1);
                nz = nz + directions(d, 2);
    
                if ny < 1 || ny > N || nz < 1 || nz > N
                    breakloop = true; % stop if bounds are reached
                    break
                elseif 2*abs(ny - center) <= Nacs && 2*abs(nz - center) <= Nacs
                    % ACS index
                    idcs = [idcs; ny, nz];
                elseif ismember([ny,nz],CAIPI_idcs,'rows')
                    % CAIPI index
                    idcs = [idcs; ny, nz];
                end
            end

            % check for loop break
            if breakloop
                break
            end

            % increase step size after completing vertical moves
            if mod(d, 2) == 0
                step_size = step_size + 1;
            end

        end

        % check for loop break
        if breakloop
            break
        end

    end

    % make sure no indicies are repeated
    idcs = unique(idcs,'rows','stable');

end

