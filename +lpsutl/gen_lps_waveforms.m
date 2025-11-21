function [gx,gy,g_amp,t_ramp,k_in,k_out] = gen_lps_waveforms(varargin)
% generates gradient and rf waveforms for a single looping star TR, also
% returns trajectory
% by David Frey (djfrey@umich.edu)
%
% inputs:
% sys - pulseq mr system structure
% fov - field of view (cm)
% N - nominal 3D matrix size
% nspokes - number of spokes (or rf subpulses)
% tseg - time per spoke (us)
%
% outputs:
% gx - x gradient waveforms for each spoke (Hz/m)
% gy - y gradient waveforms for each spoke (Hz/m)
% g_amp - gradient amplitude
% t_ramp - ramp time (s)
% k_in - kspace spoke-in trajectory (1/cm)
% k_out - kspace spoke-out trajectory (1/cm)
%

    % define default arguments
    arg.sys = mr.opts; % system structure
    arg.fov = 20; % fov (cm)
    arg.N = 128; % nominal matrix size
    arg.nspokes = 23; % number of lps spokes
    arg.t_seg = 1120; % time/segment (us)

    % parse inputs
    arg = vararg_pair(arg,varargin);

    % check that segment/rf widths are valid with given raster time
    dt_g = arg.sys.gradRasterTime*1e6; % (us)
    dt_adc = arg.sys.adcRasterTime*1e6; % (us)
    dt_rf = arg.sys.rfRasterTime*1e6; % (us)
    if mod(arg.t_seg,dt_g) > 0 || mod(arg.t_seg,dt_rf) > 0
        error('tseg must be divisible by rf & gradient raster times');
    end

    % calculate loop frequency
    f_loop = 1 / (arg.nspokes * arg.t_seg*1e-6); % (Hz)

    % calculate gradient amplitude & slew rate
    g_amp = pi*arg.N / (arg.t_seg*1e-6) / ...
        (arg.nspokes * arg.fov*1e-2 * 2 * sin(pi/arg.nspokes)); % (Hz/m)
    s_amp = 2*pi * f_loop * g_amp; % (Hz/m/s)
    assert(g_amp <= arg.sys.maxGrad, ...
        'gradient amp exceeds limit with given parameters')
    assert(s_amp <= arg.sys.maxSlew, ...
        'slew rate exceeds limit with given parameters')

    % construct looping star gradients (nseg x nspokes)
    nseg_g = round(arg.t_seg/dt_g);
    n_g = (0:arg.nspokes*nseg_g-1) + 0.5; % sampled at center of each raster
    gx = g_amp * cos(2*pi * f_loop * n_g*dt_g*1e-6); % (Hz/m)
    gy = g_amp * sin(2*pi * f_loop * n_g*dt_g*1e-6); % (Hz/m)
    gx = reshape(gx,nseg_g,arg.nspokes);
    gy = reshape(gy,nseg_g,arg.nspokes);

    % get the gradients at ADC times to calculate k-space trajectory
    nseg_adc = round(arg.t_seg/dt_adc);
    n_adc = 0:arg.nspokes*nseg_adc-1;
    gx_adc = interp1(dt_g*n_g,gx(:),dt_adc*n_adc,'nearest','extrap');
    gx_adc(1) = g_amp; gx_adc(end) = g_amp; % force 1 value at beginning and end
    gy_adc = interp1(dt_g*n_g,gy(:),dt_adc*n_adc,'nearest','extrap');
    gy_adc(1) = 0; gy_adc(end) = 0; % force 0 value at beginning and end

    % calculate ramp time
    t_ramp = g_amp / s_amp; % minimum to ensure slew is no greater (s)
    t_ramp = dt_g*1e-6 * ceil(t_ramp / (dt_g*1e-6)); % round to gradient raster (s)
    t_ramp = dt_rf*1e-6 * ceil(t_ramp / (dt_rf*1e-6)); % round to rf raster (s)

    % calculate the full kspace trajectory for each spoke
    k_spokes = zeros([arg.nspokes*nseg_adc,2,arg.nspokes]);
    for v = 1:arg.nspokes
        idx_spoke = mod((v-1)*nseg_adc + (0:arg.nspokes*nseg_adc-1), arg.nspokes*nseg_adc)+1;
        k_spokes(:,1,v) = cumsum(gx_adc(idx_spoke)*1e-2,1)*dt_adc*1e-6; % (1/cm)
        k_spokes(:,2,v) = cumsum(gy_adc(idx_spoke)*1e-2,1)*dt_adc*1e-6; % (1/cm)
    end
    k_spokes = permute(k_spokes,[1,3,2]); % [nseg*nspokes x nspokes x 2]

    % isolate in/out spokes
    k_out = reshape(k_spokes(1:nseg_adc,:,:),[],2); % spoke out
    ktmp = circshift(k_spokes,nseg_adc,1); % shift along fast time to get spoke in
    ktmp = circshift(ktmp,-1,2); % shift along spokes to align spokes in time
    k_in = reshape(ktmp(1:nseg_adc,:,:),[],2); % spoke in

end