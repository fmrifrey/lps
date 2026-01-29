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


function R = eul2rotm(eul,seq)
%This method is for internal use only. It may be removed in the future.

%EUL2ROTM Convert Euler angles to rotation matrix
%   R = EUL2ROTM(EUL, SEQ) converts 3D Euler angles into a rotation matrix.
%   The Euler angles are specified by the body-fixed (intrinsic) axis rotation
%   sequence, SEQ.
%
%   The following rotation sequences, SEQ, are supported: "ZYX", "ZYZ",
%   "XYZ", "ZXY", "ZXZ", "YXZ", "YXY", "YZX", "YZY", "XYX", "XZY", and
%   "XZX".
%
%   This is an internal function that does no input validation and is used
%   by user-facing functionality.
%
%   See also eul2rotm.

%   Copyright 2022-2023 The MathWorks, Inc.

%#codegen

    seq = convertStringsToChars(seq);

    R = zeros(3,3,size(eul,1),'like',eul);
    ct = cos(eul);
    st = sin(eul);

    % The parsed sequence will be in all upper-case letters and validated
    switch seq
      case 'ZYX'
        %     The rotation matrix R can be constructed as follows from
        %     input eul = [tz ty tx] and
        %     ct = cos(eul) = [cz cy cx] and
        %     st = sin(eul) = [sz sy sx]
        %
        %     R = [  cy*cz   sy*sx*cz-sz*cx    sy*cx*cz+sz*sx
        %            cy*sz   sy*sx*sz+cz*cx    sy*cx*sz-cz*sx
        %              -sy            cy*sx             cy*cx]
        %       = Rz(tz) * Ry(ty) * Rx(tx)

        cz = ct(:,1);

        cx = ct(:,3);
        cy = ct(:,2);
        sx = st(:,3);
        sy = st(:,2);
        sz = st(:,1);


        R(1,1,:) = cy.*cz;
        R(1,2,:) = sy.*sx.*cz - sz.*cx;
        R(1,3,:) = sy.*cx.*cz + sz.*sx;
        R(2,1,:) = cy.*sz;
        R(2,2,:) = sy.*sx.*sz + cz.*cx;
        R(2,3,:) = sy.*cx.*sz - cz.*sx;
        R(3,1,:) = -sy;
        R(3,2,:) = cy.*sx;
        R(3,3,:) = cy.*cx;

      case 'ZYZ'
        %     The rotation matrix R can be constructed as follows by
        %     input eul = [tz ty tz2] and
        %     ct = cos(eul) = [cz cy cz2] and
        %     st = sin(eul) = [sz sy sz2]
        %
        %     R = [  cz2*cy*cz-sz2*sz   -sz2*cy*cz-cz2*sz    sy*cz
        %            cz2*cy*sz+sz2*cz   -sz2*cy*sz+cz2*cz    sy*sz
        %                     -cz2*sy              sz2*sy       cy]
        %       = Rz(tz) * Ry(ty) * Rz(tz2)

        cz = ct(:,1);
        cy = ct(:,2);
        cz2 = ct(:,3);
        sz = st(:,1);
        sy = st(:,2);
        sz2 = st(:,3);

        R(1,1,:) = cz2.*cy.*cz - sz2.*sz;
        R(1,2,:) = -sz2.*cy.*cz - cz2.*sz;
        R(1,3,:) = sy.*cz;
        R(2,1,:) = cz2.*cy.*sz + sz2.*cz;
        R(2,2,:) = -sz2.*cy.*sz + cz2.*cz;
        R(2,3,:) = sy.*sz;
        R(3,1,:) = -cz2.*sy;
        R(3,2,:) = sz2.*sy;
        R(3,3,:) = cy;

      case 'XYZ'
        %     The rotation matrix R can be constructed as follows by
        %     input eul = [tx ty tz] and
        %     ct = cos(eul) = [cx cy cz] and
        %     st = sin(eul) = [sx sy sz]
        %
        %     R = [            cy*cz,           -cy*sz,     sy]
        %         [ cx*sz + cz*sx*sy, cx*cz - sx*sy*sz, -cy*sx]
        %         [ sx*sz - cx*cz*sy, cz*sx + cx*sy*sz,  cx*cy]
        %       = Rx(tx) * Ry(ty) * Rz(tz)

        cx = ct(:,1);
        cy = ct(:,2);
        cz = ct(:,3);
        sx = st(:,1);
        sy = st(:,2);
        sz = st(:,3);

        R(1,1,:) = cy.*cz;
        R(1,2,:) = -cy.*sz;
        R(1,3,:) = sy;
        R(2,1,:) = cx.*sz + cz.*sx.*sy;
        R(2,2,:) = cx.*cz - sx.*sy.*sz;
        R(2,3,:) = -cy.*sx;
        R(3,1,:) = sx.*sz - cx.*cz.*sy;
        R(3,2,:) = cz.*sx + cx.*sy.*sz;
        R(3,3,:) = cx.*cy;

      case 'ZXY'
        %     The rotation matrix R can be constructed as follows by
        %     input eul = [tz tx ty] and
        %     ct = cos(eul) = [cz cx cy] and
        %     st = sin(eul) = [sz sx sy]
        %
        %     R = [ cy*cz-sy*sx*sz, -sz*cx,  sy*cz+cy*sx*sz]
        %         [ cy*sz+sy*sx*cz,  cz*cx,  sy*sz-cy*sx*cz]
        %         [         -sy*cx,     sx,           cy*cx]
        %       = Rz(tz) * Rx(tx) * Ry(ty)

        cz = ct(:,1);
        cx = ct(:,2);
        cy = ct(:,3);
        sz = st(:,1);
        sx = st(:,2);
        sy = st(:,3);

        R(1,1,:) = cy.*cz - sy.*sx.*sz;
        R(1,2,:) = -sz.*cx;
        R(1,3,:) = sy.*cz + cy.*sx.*sz;
        R(2,1,:) = cy.*sz + sy.*sx.*cz;
        R(2,2,:) = cz.*cx;
        R(2,3,:) = sy.*sz - cy.*sx.*cz;
        R(3,1,:) = -sy.*cx;
        R(3,2,:) = sx;
        R(3,3,:) = cy.*cx;

      case 'ZXZ'
        %     The rotation matrix R can be constructed as follows by
        %     input eul = [tz tx tz2] and
        %     ct = cos(eul) = [cz cx cz2] and
        %     st = sin(eul) = [sz sx sz2]
        %
        %     R = [ cz2*cz-sz2*cx*sz, -sz2*cz-cz2*cx*sz,  sz*sx]
        %         [ cz2*sz+sz2*cx*cz, -sz2*sz+cz2*cx*cz, -cz*sx]
        %         [           sz2*sx,            cz2*sx,     cx]
        %       = Rz(tz) * Rx(tx) * Rz(tz2)

        cz = ct(:,1);
        cx = ct(:,2);
        cz2 = ct(:,3);
        sz = st(:,1);
        sx = st(:,2);
        sz2 = st(:,3);

        R(1,1,:) = cz2.*cz - sz2.*cx.*sz;
        R(1,2,:) = -sz2.*cz - cz2.*cx.*sz;
        R(1,3,:) = sz.*sx;
        R(2,1,:) = cz2.*sz + sz2.*cx.*cz;
        R(2,2,:) = -sz2.*sz + cz2.*cx.*cz;
        R(2,3,:) = -cz.*sx;
        R(3,1,:) = sz2.*sx;
        R(3,2,:) = cz2.*sx;
        R(3,3,:) = cx;

      case 'YXZ'
        %     The rotation matrix R can be constructed as follows by
        %     input eul = [ty tx tz] and
        %     ct = cos(eul) = [cy cx cz] and
        %     st = sin(eul) = [sy sx sz]
        %
        %     R = [  cy*cz+sy*sx*sz, -cy*sz+sy*sx*cz, sy*cx]
        %         [           sz*cx,           cz*cx,   -sx]
        %         [ -sy*cz+cy*sx*sz,  sy*sz+cy*sx*cz, cy*cx]
        %       = Ry(ty) * Rx(tx) * Rz(tz)

        cy = ct(:,1);
        cx = ct(:,2);
        cz = ct(:,3);
        sy = st(:,1);
        sx = st(:,2);
        sz = st(:,3);

        R(1,1,:) = cy.*cz + sy.*sx.*sz;
        R(1,2,:) = -cy.*sz + sy.*sx.*cz;
        R(1,3,:) = sy.*cx;
        R(2,1,:) = sz.*cx;
        R(2,2,:) = cz.*cx;
        R(2,3,:) = -sx;
        R(3,1,:) = -sy.*cz + cy.*sx.*sz;
        R(3,2,:) = sy.*sz + cy.*sx.*cz;
        R(3,3,:) = cy.*cx;

      case 'YXY'
        %     The rotation matrix R can be constructed as follows by
        %     input eul = [ty tx ty2] and
        %     ct = cos(eul) = [cy cx cy2] and
        %     st = sin(eul) = [sy sx sy2]
        %
        %     R = [  cy2*cy-sy2*cx*sy, sy*sx,  sy2*cy+cy2*cx*sy]
        %         [            sy2*sx,    cx,           -cy2*sx]
        %         [ -cy2*sy-sy2*cx*cy, cy*sx, -sy2*sy+cy2*cx*cy]
        %       = Ry(ty) * Rx(tx) * Ry(ty2)

        cy = ct(:,1);
        cx = ct(:,2);
        cy2 = ct(:,3);
        sy = st(:,1);
        sx = st(:,2);
        sy2 = st(:,3);

        R(1,1,:) = cy2.*cy - sy2.*cx.*sy;
        R(1,2,:) = sy.*sx;
        R(1,3,:) = sy2.*cy + cy2.*cx.*sy;
        R(2,1,:) = sy2.*sx;
        R(2,2,:) = cx;
        R(2,3,:) = -cy2.*sx;
        R(3,1,:) = -cy2.*sy - sy2.*cx.*cy;
        R(3,2,:) = cy.*sx;
        R(3,3,:) = -sy2.*sy + cy2.*cx.*cy;

      case 'YZX'
        %     The rotation matrix R can be constructed as follows by
        %     input eul = [ty tz tx] and
        %     ct = cos(eul) = [cy cz cx] and
        %     st = sin(eul) = [sy sz sx]
        %
        %     R = [  cy*cz, -sz*cx*cy+sy*sx,  cy*sx*sz+sy*cx]
        %         [     sz,           cz*cx,          -cz*sx]
        %         [ -sy*cz,  sy*cx*sz+cy*sx, -sy*sx*sz+cy*cx]
        %       = Ry(ty) * Rz(tz) * Rx(tx)

        cy = ct(:,1);
        cz = ct(:,2);
        cx = ct(:,3);
        sy = st(:,1);
        sz = st(:,2);
        sx = st(:,3);

        R(1,1,:) = cy.*cz;
        R(1,2,:) = -sz.*cx.*cy + sy.*sx;
        R(1,3,:) = cy.*sx.*sz + sy.*cx;
        R(2,1,:) = sz;
        R(2,2,:) = cz.*cx;
        R(2,3,:) = -cz.*sx;
        R(3,1,:) = -sy.*cz;
        R(3,2,:) = sy.*cx.*sz + cy.*sx;
        R(3,3,:) = -sy.*sx.*sz + cy.*cx;

      case 'YZY'
        %     The rotation matrix R can be constructed as follows by
        %     input eul = [ty tz ty2] and
        %     ct = cos(eul) = [cy cz cy2] and
        %     st = sin(eul) = [sy sz sy2]
        %
        %     R = [  cy2*cz*cy-sy2*sy, -cy*sz,  sy2*cz*cy+cy2*sy]
        %         [            cy2*sz,     cz,            sy2*sz]
        %         [ -cy2*cz*sy-sy2*cy,  sy*sz, -sy2*cz*sy+cy2*cy]
        %       = Ry(ty) * Rz(tz) * Ry(ty2)

        cy = ct(:,1);
        cz = ct(:,2);
        cy2 = ct(:,3);
        sy = st(:,1);
        sz = st(:,2);
        sy2 = st(:,3);

        R(1,1,:) = cy2.*cz.*cy - sy2.*sy;
        R(1,2,:) = -cy.*sz;
        R(1,3,:) = sy2.*cz.*cy + cy2.*sy;
        R(2,1,:) = cy2.*sz;
        R(2,2,:) = cz;
        R(2,3,:) = sy2.*sz;
        R(3,1,:) = -cy2.*cz.*sy - sy2.*cy;
        R(3,2,:) = sy.*sz;
        R(3,3,:) = -sy2.*cz.*sy + cy2.*cy;

      case 'XYX'
        %     The rotation matrix R can be constructed as follows by
        %     input eul = [tx ty tx2] and
        %     ct = cos(eul) = [cx cy cx2] and
        %     st = sin(eul) = [sx sy sx2]
        %
        %     R = [     cy,           sx2*sy,            cx2*sy]
        %         [  sy*sx, cx2*cx-sx2*cy*sx, -sx2*cx-cx2*cy*sx]
        %         [ -sy*cx, cx2*sx+sx2*cy*cx, -sx2*sx+cx2*cy*cx]
        %       = Rx(tx) * Ry(ty) * Rx(tx2)

        cx = ct(:,1);
        cy = ct(:,2);
        cx2 = ct(:,3);
        sx = st(:,1);
        sy = st(:,2);
        sx2 = st(:,3);

        R(1,1,:) = cy;
        R(1,2,:) = sx2.*sy;
        R(1,3,:) = cx2.*sy;
        R(2,1,:) = sy.*sx;
        R(2,2,:) = cx2.*cx - sx2.*cy.*sx;
        R(2,3,:) = -sx2.*cx - cx2.*cy.*sx;
        R(3,1,:) = -sy.*cx;
        R(3,2,:) = cx2.*sx + sx2.*cy.*cx;
        R(3,3,:) = -sx2.*sx + cx2.*cy.*cx;

      case 'XZY'
        %     The rotation matrix R can be constructed as follows by
        %     input eul = [tx tz ty] and
        %     ct = cos(eul) = [cx cz cy] and
        %     st = sin(eul) = [sx sz sy]
        %
        %     R = [          cy*cz,   -sz,          sy*cz]
        %         [ sz*cx*cy+sy*sx, cz*cx, sy*cx*sz-cy*sx]
        %         [ cy*sx*sz-sy*cx, cz*sx, sy*sx*sz+cy*cx]
        %       = Rx(tx) * Rz(tz) * Ry(ty)

        cx = ct(:,1);
        cz = ct(:,2);
        cy = ct(:,3);
        sx = st(:,1);
        sz = st(:,2);
        sy = st(:,3);

        R(1,1,:) = cy.*cz;
        R(1,2,:) = -sz;
        R(1,3,:) = sy.*cz;
        R(2,1,:) = sz.*cx.*cy + sy.*sx;
        R(2,2,:) = cz.*cx;
        R(2,3,:) = sy.*cx.*sz - cy.*sx;
        R(3,1,:) = cy.*sx.*sz - sy.*cx;
        R(3,2,:) = cz.*sx;
        R(3,3,:) = sy.*sx.*sz + cy.*cx;

      case 'XZX'
        %     The rotation matrix R can be constructed as follows by
        %     input eul = [tx tz tx2] and
        %     ct = cos(eul) = [cx cz cx2] and
        %     st = sin(eul) = [sx sz sx2]
        %
        %     R = [    cz,          -cx2*sz,            sx2*sz]
        %         [ sz*cx, cx2*cz*cx-sx2*sx, -sx2*cz*cx-cx2*sx]
        %         [ sz*sx, cx2*cz*sx+sx2*cx, -sx2*cz*sx+cx2*cx]
        %       = Rx(tx) * Rz(tz) * Rx(tx2)

        cx = ct(:,1);
        cz = ct(:,2);
        cx2 = ct(:,3);
        sx = st(:,1);
        sz = st(:,2);
        sx2 = st(:,3);

        R(1,1,:) = cz;
        R(1,2,:) = -cx2.*sz;
        R(1,3,:) = sx2.*sz;
        R(2,1,:) = sz.*cx;
        R(2,2,:) = cx2.*cz.*cx - sx2.*sx;
        R(2,3,:) = -sx2.*cz.*cx - cx2.*sx;
        R(3,1,:) = sz.*sx;
        R(3,2,:) = cx2.*cz.*sx + sx2.*cx;
        R(3,3,:) = -sx2.*cz.*sx + cx2.*cx;
    end

end
