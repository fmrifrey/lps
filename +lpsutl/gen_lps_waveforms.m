function [g,rf,t_ramp,k_in,k_out] = gen_lps_waveforms(varargin)
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
% trf - rf hard pulse width (us)
% fa - flip angle (deg)
% plotwavs - option to plot the waveforms
%
% outputs:
% g - gradient waveforms (Hz/m)
% rf - rf waveform (Hz)
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
    arg.t_rf = 12; % time/rf pulse (us)
    arg.fa = 3; % rf flip angle (deg)
    arg.plotwavs = false; % option to plot the waveforms

    % parse inputs
    arg = vararg_pair(arg,varargin);

    % check that segment/rf widths are valid with given raster time
    dt_g = round(arg.sys.gradRasterTime*1e6); % (us)
    dt_rf = round(arg.sys.rfRasterTime*1e6); % (us)
    if mod(arg.t_seg,dt_g) > 0 || mod(arg.t_seg,dt_rf) > 0
        error('tseg must be divisible by rf & gradient raster times');
    end
    if mod(arg.t_rf,dt_rf) > 0
        error('trf must be divisible by rf raster time');
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

    % construct looping star gradients
    nseg_g = round(arg.t_seg/dt_g);
    n_g = 0:2*arg.nspokes*nseg_g-1;
    gx = g_amp * cos(2*pi * f_loop * n_g*dt_g*1e-6); % (Hz/m)
    gy = g_amp * sin(2*pi * f_loop * n_g*dt_g*1e-6); % (Hz/m)

    % calculate rf amplitude
    rf_amp = arg.fa / (360 * arg.t_rf*1e-6); % (Hz)
    assert(rf_amp <= arg.sys.maxB1, ...
        'rf amp exceeds limit with given parameters')

    % calculate ramp time
    t_ramp = g_amp / s_amp; % minimum to ensure slew is no greater (s)
    t_ramp = dt_g*1e-6 * ceil(t_ramp / (dt_g*1e-6)); % round to gradient raster (s)
    t_ramp = dt_rf*1e-6 * ceil(t_ramp / (dt_rf*1e-6)); % round to rf raster (s)
    
    % construct rf burst pulse
    nseg_rf = round(arg.t_seg/dt_rf);
    n_rf = 0:(arg.nspokes-1)*nseg_rf+round(arg.t_rf/dt_rf)-1;
    rf = rf_amp * (mod(n_rf,nseg_rf) < round(arg.t_rf/dt_rf)) .* (n_rf < arg.nspokes*nseg_rf);

    % append the ramp-up to the gradients
    ramp_up = linspace(0,1,round(t_ramp*1e6/dt_g)); % (Hz/m)
    gx = [ramp_up*gx(1), gx, (1-ramp_up)*gx(end)];
    gy = [ramp_up*gy(1), gy, (1-ramp_up)*gy(end)];
    g = [gx(:),gy(:)];

    % calculate the full kspace trajectory for each spoke
    k_spokes = zeros([arg.nspokes*nseg_g,2,arg.nspokes]);
    for v = 1:arg.nspokes
        idx_spoke = round(t_ramp*1e6/dt_g) + (v-1)*nseg_g + (1:arg.nspokes*nseg_g);
        k_spokes(:,:,v) = cumsum(g(idx_spoke,:)*1e-2,1)*dt_g*1e-6; % (1/cm)
    end
    k_spokes = permute(k_spokes,[1,3,2]); % [nseg*nspokes x nspokes x 2]

    % isolate in/out spokes
    k_out = reshape(k_spokes(1:nseg_g,:,:),[],2); % spoke out
    ktmp = circshift(k_spokes,nseg_g,1); % shift along fast time to get spoke in
    ktmp = circshift(ktmp,-1,2); % shift along spokes to align spokes in time
    k_in = reshape(ktmp(1:nseg_g,:,:),[],2); % spoke in

    % plot the waveforms
    if arg.plotwavs
        yyaxis left
        plot(1e-3*dt_g*(0:size(g,1)-1),g);
        ylabel('gradient amp (Hz/m)')
        yyaxis right
        plot(t_ramp*1e3+1e-3*dt_rf*(0:length(rf)-1),rf)
        ylabel('rf amp (Hz)');
        xlabel('time (ms)')
        title('looping star waveforms');
    end

end