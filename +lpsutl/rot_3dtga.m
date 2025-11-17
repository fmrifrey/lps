function R = rot_3dtga(prjn,intn)
% returns rotation matrix based on 3D TINY golden angle rotation sequence
% reference:
% Fyrdahl, A., Holst, K., Caidahl, K. et al. Generalization of
% three-dimensional golden-angle radial acquisition to reduce eddy current
% artifacts in bSSFP CMR imaging. Magn Reson Mater Phy 34, 109–118 (2021).
% https://doi.org/10.1007/s10334-020-00859-z
% by David Frey (djfrey@umich.edu)
%
% inputs:
% prjn - projection (thru-plane rotation) index
% intn - interleaf (in-plane rotation) index
%
% outputs:
% R - 3D golden angle rotation matrix
%

    % interleaf rotation: rotate by golden angle about z
    r_int = (intn-1)*pi*(3 - sqrt(5));
    R_int = eul2rotm([r_int,0,0],'ZYX');
    
    % projection rotation: rotate by 3D golden angles
    % ref: Generalization of three-dimensional golden-angle radial acquisition
    % to reduce eddy current artifacts in bSSFP CMR imaging (A, Fyrdahl et. al)
    phi1 = 0.4656; phi2 = 0.6823; % 3D golden ratios
    rp_prj = acos(mod((prjn-1) * phi1, 2)-1) + pi; % polar angle
    ra_prj = 2*pi*((prjn-1) * phi2); % azimuthal angle
    R_prj = eul2rotm([rp_prj,0,ra_prj],'ZYX');
    
    % multiply the rotations
    R = R_prj * R_int;

end

% define eul2rotm as to not require user to install robotics toolbox
function rotm = eul2rotm(eul,sequence)

    % default sequence
    if nargin > 2 || isempty(sequence)
        sequence = 'ZYX';
    else
        sequence = char(sequence);
    end

    % x rotation matrix
    eul_x = eul(char(lower(sequence))=='x');
    R.X = [cos(eul_x), -sin(eul_x), 0;
        sin(eul_x), cos(eul_x), 0;
        0, 0, 1];
    
    % y rotation matrix
    eul_y = eul(char(lower(sequence))=='y');
    R.Y = [cos(eul_y), 0, sin(eul_y);
        0, 1, 0;
        -sin(eul_y), 0, cos(eul_y)];
    
    % z rotation matrix
    eul_z = eul(char(lower(sequence))=='z');
    R.Z = [1, 0, 0;
        0, cos(eul_z), -sin(eul_z);
        0, sin(eul_z), cos(eul_z)];
    
    % multiply rotation matrices
    rotm = eye(3);
    for i = 3:-1:1
        rotm = R.(sequence(i)) * rotm;
    end

end

