function [Hs_in,Hs_out] = lps_setup_filters(Fs_in,Fs_out,cutoff,rolloff)
% sets up Fermi-shape filter for echo-in/out on kspace filter as a diagonal
% fatrix
% by David Frey (djfrey@umich.edu)
%
% inputs:
% Fs_in - cell array of spoke-in NUFFT operators for each kspace volume
% Fs_out - cell array of spoke-out NUFFT operators for each kspace volume
% cutoff - kspace cutoff (as fraction of kmax)
% rolloff - kspace filter rolloff (as fraction of kmax)
%
% outputs:
% Hs_in - spoke-in filter operator (diagonal fatrix)
% Hs_out - spoke-out filter operator (diagonal fatrix)
%

    % get number of volumes
    nvol = length(Fs_in);

    % create kspace filter
    kfilt = @(w) (w < pi*cutoff) .* 1./(1 + exp(2*pi*(w - pi*cutoff)/(pi*rolloff)));
    
    % loop through volumes and apply the filter
    Hs_in = cell(nvol,1);
    Hs_out = cell(nvol,1);
    for ivol = 1:nvol
        % get trajectory for current vol
        omega_in = Fs_in{ivol}.arg.arg{1};
        omega_out = Fs_out{ivol}.arg.arg{1};

        % calculate filter
        filt_in = kfilt(vecnorm(omega_in,2,2));
        filt_out = kfilt(vecnorm(omega_out,2,2));

        % diagonalize
        Hs_in{ivol} = Gdiag(filt_in);
        Hs_out{ivol} = Gdiag(filt_out);
    end

end

