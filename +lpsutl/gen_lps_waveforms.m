function [g,rf,rf_del,k_in,k_out] = gen_lps_waveforms(varargin)
% generates gradient and rf waveforms for a single looping star TR, also
% returns trajectory
% by David Frey (djfrey@umich.edu)
%
% inputs:
% fov - field of view (cm)
% N - nominal 3D matrix size
% nspokes - number of spokes (or rf subpulses)
% tseg - time per spoke (us)
% trf - rf hard pulse width (us)
% fa - flip angle (flip)
% gmax - max gradient amplitude (G/cm)
% smax - max slew rate (G/cm/s)
% dt - raster time (us)
% plotwavs - option to plot the waveforms
%
% outputs:
% g - gradient waveforms (G/cm)
% rf - rf waveform (G)
% rf_del - number of samples to delay the rf waveform
% k_in - kspace spoke-in trajectory (1/cm)
% k_out - kspace spoke-out trajectory (1/cm)
%

    % define default arguments
    arg.fov = 20; % fov (cm)
    arg.N = 128; % nominal matrix size
    arg.nspokes = 23; % number of lps spokes
    arg.tseg = 1120; % time/segment (us)
    arg.trf = 12; % time/rf pulse (us)
    arg.fa = 4; % rf flip angle (deg)
    arg.gmax = 4; % max gradient amplitude (G/cm)
    arg.smax = 500; % max slew rate (G/cm/s)
    arg.dt = 4; % raster time (us)
    arg.plotwavs = false; % option to plot the waveforms

    gam = 4258; % GMR of H+ (Hz/G)

    % parse inputs
    arg = vararg_pair(arg,varargin);

    % check that segment/rf widths are valid with given raster time
    if mod(arg.tseg,arg.dt) > 0
        error('tseg must be divisible by raster time');
    end
    if mod(arg.trf,arg.dt) > 0
        error('trf must be divisible by raster time');
    end

    % calculate loop velocity
    om = 2*pi / (arg.nspokes * arg.tseg*1e-6);

    % calculate gradient amplitude & slew rate
    g_amp = pi*arg.N / gam / (arg.tseg*1e-6) / ...
        (arg.nspokes * arg.fov * 2 * sin(pi/arg.nspokes));
    s_amp = om * g_amp;
    assert(g_amp <= arg.gmax, ...
        'gradient amp exceeds limit with given parameters')
    assert(s_amp <= arg.smax, ...
        'slew rate exceeds limit with given parameters')

    % construct looping star gradients
    nseg = round(arg.tseg/arg.dt);
    n = 0:2*arg.nspokes*nseg-1;
    gx = g_amp * cos(om * n*arg.dt*1e-6);
    gy = g_amp * sin(om * n*arg.dt*1e-6);

    % calculate rf amplitude
    rf_amp = arg.fa / (360 * gam * arg.trf); % (G)

    % calculate ramp-up
    nramp = ceil(g_amp / arg.smax / (arg.dt*1e-6));
    ramp_up = linspace(0,1,nramp);
    nramp = length(ramp_up);
    
    % construct rf burst pulse
    nrf = round(arg.trf/arg.dt);
    n = 0:(arg.nspokes-1)*nseg+nrf-1;
    rf = rf_amp * (mod(n,nseg) < nrf) .* (n < arg.nspokes*nseg);
    rf_del = nramp;

    % append the ramp
    gx = [ramp_up*gx(1), gx, (1-ramp_up)*gx(end)];
    gy = [ramp_up*gy(1), gy, (1-ramp_up)*gy(end)];
    g = [gx(:),gy(:)];

    % calculate the full kspace trajectory for each spoke
    k_spokes = zeros([arg.nspokes*nseg,2,arg.nspokes]);
    for v = 1:arg.nspokes
        idx_spoke = nramp + (v-1)*nseg + (1:arg.nspokes*nseg);
        k_spokes(:,:,v) = gam*cumsum(g(idx_spoke,:),1)*arg.dt*1e-6;
    end
    k_spokes = permute(k_spokes,[1,3,2]); % [nseg*nspokes x nspokes x 2]

    % isolate in/out spokes
    k_out = reshape(k_spokes(1:nseg,:,:),[],2); % spoke out
    ktmp = circshift(k_spokes,nseg,1); % shift along fast time to get spoke in
    ktmp = circshift(ktmp,-1,2); % shift along spokes to align spokes in time
    k_in = reshape(ktmp(1:nseg,:,:),[],2); % spoke in

    % plot the waveforms
    if arg.plotwavs
        yyaxis left
        plot(1e-3*arg.dt*(0:size(g,1)-1),g);
        ylabel('gradient amp (G/cm)')
        yyaxis right
        plot(1e-3*arg.dt*(rf_del + (0:length(rf)-1)),rf)
        ylabel('rf amp (G)');
        xlabel('time (ms)')
        title('looping star waveforms');
    end

end