function [g,g0,rf,t_ramp,k_in,k_out] = gen_lps_waveforms(varargin)
% generates gradient and rf waveforms for a single looping star TR, also
% returns trajectory
% by David Frey (djfrey@umich.edu)
%
% inputs:
% sys - pulseq mr system structure
% dwell - adc dwell time
% fov - field of view (cm)
% N - nominal 3D matrix size
% nspokes - number of spokes (or rf subpulses)
% nechoes - number of echoes (GREs, not including FID)
% tseg - time per spoke (us)
% trf - rf hard pulse width (us)
% fa - flip angle (deg)
% plotwavs - option to plot the waveforms
%
% outputs:
% g - gradient waveforms (Hz/m)
% g0 - gradient amplitudes at 0 and each echo time (Hz/m)
% rf - rf waveform (Hz)
% t_ramp - ramp time (s)
% k_in - kspace spoke-in trajectory (1/cm)
% k_out - kspace spoke-out trajectory (1/cm)
%

    % define default arguments
    arg.sys = mr.opts; % system structure
    arg.dwell = 4; % ADC sample rate (us)
    arg.fov = 20; % fov (cm)
    arg.N = 128; % nominal matrix size
    arg.nspokes = 23; % number of lps spokes
    arg.nechoes = 2; % number echoes
    arg.t_seg = 1120; % time/segment (us)
    arg.t_rf = 12; % time/rf pulse (us)
    arg.fa = 3; % rf flip angle (deg)
    arg.C = [0, 0, 0; 1, 0, 0; 0, 1, 1]; % fourier basis coefficient matrix
    arg.plotwavs = false; % option to plot the waveforms

    % parse inputs
    arg = vararg_pair(arg,varargin);

    % check that segment/rf widths are valid with given raster time
    dt_g = arg.sys.gradRasterTime*1e6; % (us)
    dt_adc = arg.dwell; % (us)
    dt_rf = arg.sys.rfRasterTime*1e6; % (us)
    if mod(arg.t_seg,dt_g) > 0 || mod(arg.t_seg,dt_rf) > 0
        error('tseg must be divisible by rf & gradient raster times');
    end
    if mod(arg.t_rf,dt_rf) > 0
        error('trf must be divisible by rf raster time');
    end

    % create the k-space and gradient basis functions
    nf = (size(arg.C, 1) - 1) / 2;
    B = @(t) fourier_series_basis(t, nf, arg.nspokes * arg.t_seg*1e-6 / nf);
    dB = @(t) fourier_series_basis_d1(t, nf, arg.nspokes * arg.t_seg*1e-6 / nf);

    % determine the scale factor based on minimum segment radius
    min_seg_dist = Inf;
    for s = 1:arg.nspokes*(arg.nechoes + 1)
        seg_dist = norm((B(s*arg.t_seg*1e-6) - B((s-1)*arg.t_seg*1e-6)) * arg.C, 2);
        min_seg_dist = min(min_seg_dist, seg_dist);
    end
    k_scale_factor = arg.N / (arg.fov*1e-2) / min_seg_dist;

    % create the gradient waveform function
    g_fun = @(t) k_scale_factor * dB(t) * arg.C; % (Hz/m)

    % construct looping star gradients
    nseg_g = round(arg.t_seg/dt_g);
    n_g = (0:(arg.nechoes+1)*arg.nspokes*nseg_g-1) + 0.5;
    g = g_fun(n_g*dt_g*1e-6);
    g0 = g_fun((0:arg.nechoes+1)'*arg.nspokes*arg.t_seg*1e-6);

    % get max gradient and amplitudes
    g_max = max(abs(g),[],'all'); % Hz/cm
    s_max = max(abs(diff(g)/(dt_g*1e-6)),[],'all'); % Hz/cm/s

    % get the gradients at ADC times to calculate k-space trajectory
    nseg_adc = round(arg.t_seg/dt_adc);
    n_adc = 0:(arg.nechoes+1)*arg.nspokes*nseg_adc-1;
    g_adc = g_fun(n_adc*dt_adc*1e-6);

    % calculate rf amplitude
    rf_amp = arg.fa / (360 * arg.t_rf*1e-6); % (Hz)
    assert(rf_amp <= arg.sys.maxB1, ...
        'rf amp exceeds limit with given parameters')

    % calculate ramp time
    t_ramp = g_max / s_max; % minimum to ensure slew is no greater (s)
    t_ramp = dt_g*1e-6 * ceil(t_ramp / (dt_g*1e-6)); % round to gradient raster (s)
    t_ramp = dt_rf*1e-6 * ceil(t_ramp / (dt_rf*1e-6)); % round to rf raster (s)
    
    % construct rf burst pulse
    nseg_rf = round(arg.t_seg/dt_rf);
    n_rf = 0:(arg.nspokes-1)*nseg_rf+round(arg.t_rf/dt_rf)-1;
    rf = rf_amp * (mod(n_rf,nseg_rf) < round(arg.t_rf/dt_rf)) .* (n_rf < arg.nspokes*nseg_rf);

    % calculate the full kspace trajectory for each spoke
    k_spokes = zeros([arg.nspokes*nseg_adc,3,arg.nspokes]);
    for v = 1:arg.nspokes
        idx_spoke = mod((v-1)*nseg_adc + (0:arg.nspokes*nseg_adc-1), arg.nspokes*nseg_adc)+1;
        k_spokes(:,:,v) = cumsum(g_adc(idx_spoke,:)*1e-2,1)*dt_adc*1e-6; % (1/cm)
    end
    k_spokes = permute(k_spokes,[1,3,2]); % [nseg*nspokes x nspokes x 3]

    % isolate in/out spokes
    k_out = reshape(k_spokes(1:nseg_adc,:,:),[],3); % spoke out
    ktmp = circshift(k_spokes,nseg_adc,1); % shift along fast time to get spoke in
    ktmp = circshift(ktmp,-1,2); % shift along spokes to align spokes in time
    k_in = reshape(ktmp(1:nseg_adc,:,:),[],3); % spoke in

    % plot the waveforms
    if arg.plotwavs
        yyaxis left
        plot(1e-3*dt_g*(0:size(g,1)-1),g);
        ylabel('gradient amp (Hz/m)')
        yyaxis right
        plot(1e-3*dt_rf*(0:length(rf)-1),rf)
        ylabel('rf amp (Hz)');
        xlabel('time (ms)')
        title('looping star waveforms');
    end

end

function B = fourier_series_basis(t,nf,dt)

    f = (1:nf) / (nf*dt); % fourier series frequencies
    B = ones(length(t), 2*nf+1); % DC in col 1

    for i = 1:nf
        omega = 2 * pi * f(i);
        B(:, 2*(i-1) + 2) = cos(omega * t);
        B(:, 2*(i-1) + 3) = sin(omega * t);
    end

end

function dB = fourier_series_basis_d1(t,nf,dt)

    f = (1:nf) / (nf*dt); % fourier series frequencies
    dB = zeros(length(t), 2*nf+1); % DC in col 1

    for i = 1:nf
        omega = 2 * pi * f(i);
        dB(:, 2*(i-1) + 2) = -omega * sin(omega * t);
        dB(:, 2*(i-1) + 3) = omega * cos(omega * t);
    end

end