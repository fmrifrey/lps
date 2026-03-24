function lps_write_seq(varargin)
% create the sequence files for looping star in pulseq
% by David Frey (djfrey@umich.edu)
%
% inputs:
% sys - pulseq mr system structure
% dir - output directory name
% tr - repetition time (ms)
% fov - field of view (cm)
% N_nom - nominal 3D matrix size
% dummyshots - number of dummy shots (disdaqs) to play
% nrep - number of rotation sequence repetitions
% nint - number of interleaves (2D in-plane rotations)
% nprj - number of projections (3D thru-plane rotations)
% nspokes - number of looping star spokes
% t_seg - time per spoke (us)
% t_rf - rf hard pulse width (us)
% fa - flip angle (deg)
% plotseq - option to plot the sequence
% pislquant - number of TRs to use for prescan
% writepge - option to convert seq to pge file
%
% output files:
% lps.seq file - seq file for pulseq
% lps.pge file - seq file for pge2 interpreter
% seq_args.h5 - .h5 file containing copy of input arguments
%

    % set default arguments
    arg.sys = lpsutl.get_sys_defaults('ge'); % pulseq mr system structure
    arg.dir = pwd; % output directory name
    arg.tr = 100; % repetition time (ms)
    arg.fov = 20; % fov (cm)
    arg.N_nom = 128; % 3D matrix size
    arg.dummyshots = 20; % number of dummy shots
    arg.nrep = 1; % number of rotation sequence repetitions
    arg.nint = 1; % number of interleaves (2D in-plane rotations)
    arg.nprj = 16; % number of projections (3D thru-plane rotations)
    arg.nechoes = 1; % number of echoes to acquire (GREs only, no FID)
    arg.nspokes = 23; % number of lps spokes
    arg.psd_rf_wait = 100e-6; % RF–gradient delay (s), GE scanner-specific
    arg.psd_grd_wait = 100e-6;   % ADC–gradient delay (s), scanner-specific
    arg.coil = 'xrm'; % GE coil model
    arg.dwell = 4; % ADC sample rate (us)
    arg.t_seg = 1120; % time/segment (us)
    arg.t_rf = 16; % time/rf pulse (us)
    arg.fa = 4; % rf flip angle (deg)
    arg.C = [0, 0, 0; 1, 0, 0; 0, 1, 0]; % fourier basis coefficient matrix
    arg.plotseq = false; % option to plot the sequence
    arg.pislquant = 1; % number of TRs to use for prescan
    arg.writepge = true; % option to convert seq to pge file
    
    % parse arguments
    arg = vararg_pair(arg,varargin);
    
    % initialize sequence object
    seq = mr.Sequence(arg.sys);
    warning('OFF', 'mr:restoreShape');

    % create looping star waveforms
    [g_wav,g0,rf_wav,t_ramp] = lpsutl.gen_lps_waveforms( ...
        'sys', arg.sys, ... % pulseq mr system structure
        'fov', arg.fov, ... % fov (cm)
        'dwell', arg.dwell, ... % adc dwell time (us)
        'N', arg.N_nom, ... % nominal matrix size
        'nspokes', arg.nspokes, ... % number of lps spokes
        'nechoes', arg.nechoes, ... % number of echoes
        'C', arg.C, ... % fourier coefficient matrix
        't_seg', arg.t_seg, ... % time/segment
        't_rf', arg.t_rf, ... % rf pulse width
        'fa', arg.fa ... % rf flip angle (deg)
        );
    g_wav = reshape(g_wav, [], arg.nechoes+1, 3);

    % create gradient ramp waveform
    nramp = ceil(t_ramp/arg.sys.gradRasterTime);
    ramp_up = 1/nramp * ((0:nramp - 1).' + 0.5);
    
    % create gradient waveform objects for fID portion
    g_wav_fid = [ramp_up.*g0(1,:); squeeze(g_wav(:,1,:))];
    gx_fid = mr.makeArbitraryGrad('x', 0.99*g_wav_fid(:,1).', ...
        'system', arg.sys, ...
        'first', 0, ...
        'last', 0.99*g0(2,1));
    gy_fid = mr.makeArbitraryGrad('y', 0.99*g_wav_fid(:,2).', ...
        'system', arg.sys, ...
        'first', 0, ...
        'last', 0.99*g0(2,2));
    gz_fid = mr.makeArbitraryGrad('z', 0.99*g_wav_fid(:,3).', ...
        'system', arg.sys, ...
        'first', 0, ...
        'last', 0.99*g0(2,3));

    % create gradient waveform objects for GRE portion
    g_wav_gre = [reshape(g_wav(:,2:end,:),[],3); (1-ramp_up).*g0(end,:)];
    gx_gre = mr.makeArbitraryGrad('x', 0.99*g_wav_gre(:,1).', ...
        'system', arg.sys, ...
        'first', 0.99*g0(2,1), ...
        'last', 0);
    gy_gre = mr.makeArbitraryGrad('y', 0.99*g_wav_gre(:,2).', ...
        'system', arg.sys, ...
        'first', 0.99*g0(2,2), ...
        'last', 0);
    gz_gre = mr.makeArbitraryGrad('z', 0.99*g_wav_gre(:,3).', ...
        'system', arg.sys, ...
        'first', 0.99*g0(2,3), ...
        'last', 0);

    % only include non-zero gradients in gradient object cell array
    grads_fid = {gx_fid, gy_fid, gz_fid};
    grads_gre = {gx_gre, gy_gre, gz_gre};
    grads_fid = grads_fid(cellfun(@(g) max(abs(g.waveform)) > 0, grads_fid));
    grads_gre = grads_gre(cellfun(@(g) max(abs(g.waveform)) > 0, grads_gre));

    % create rf waveform
    rf_delay = t_ramp - arg.t_rf*1e-6/2;
    rf_delay = arg.sys.rfRasterTime*round(rf_delay/arg.sys.rfRasterTime);
    if rf_delay < arg.sys.rfDeadTime
        error('rf delay is less than deadtime')
    end
    rf = mr.makeArbitraryRf( ...
        rf_wav, arg.nspokes*arg.fa*pi/180, ...
        'delay', rf_delay, ...
        'use', 'excitation', ...
        'system', arg.sys);

    % create ADC
    nseg = round(arg.t_seg/arg.dwell);
    acq_len = arg.sys.adcSamplesDivisor * ceil(arg.nspokes*nseg*arg.nechoes / ...
        arg.sys.adcSamplesDivisor);
    adc = mr.makeAdc(acq_len, 'system', arg.sys, 'dwell', arg.dwell*1e-6);

    % create TR delay
    delay_time = arg.tr*1e-3;
    delay_time = delay_time - mr.calcDuration(gx_fid);
    delay_time = delay_time - mr.calcDuration(gx_gre);
    delay_time = arg.sys.blockDurationRaster*round(delay_time/arg.sys.blockDurationRaster);
    tr_delay = mr.makeDelay(delay_time);
    
    % define sequence blocks
    for i = (-arg.dummyshots-arg.pislquant+1):arg.nrep*arg.nprj*arg.nint
        isDummyTR = i <= -arg.pislquant;
        isReceiveGainCalibrationTR = i < 1 & i > -arg.pislquant;

        if i > 0
            % get interleaf, projection, and repetition index
            [iint,iprj,~] = ind2sub([arg.nint,arg.nprj,arg.nrep],i);
    
            % rotate the gradients based on a 3DTGA rotation sequence
            R = lpsutl.rot_3dtga(iprj, iint);
        else
            % rotate as if not a dummy shot
            R = lpsutl.rot_3dtga(i+arg.dummyshots+arg.pislquant+1, 1);
        end
        
        % R = eye(3);
        
        % make rotation object from rotation matrix
        rot = mr.makeRotation(R);

        % make TRID label
        lbl = mr.makeLabel('SET','TRID', 1 + isDummyTR + 2*isReceiveGainCalibrationTR);



        % write fID portion to sequence
        seq.addBlock(rf, grads_fid{:}, rot, lbl);

        % write GRE portion to sequence
        if isDummyTR
            seq.addBlock(grads_gre{:}, rot); % no adc
        else
            seq.addBlock(adc, grads_gre{:}, rot);
        end

        % add tr delay
        seq.addBlock(tr_delay);

    end
    
    % check whether the timing of the sequence is correct
    [ok, error_report] = seq.checkTiming;
    if (ok)
        fprintf('Timing check passed successfully\n');
    else
        fprintf('Timing check failed! Error listing follows:\n');
        fprintf([error_report{:}]);
        fprintf('\n');
    end
    
    % write out sequence and save args
    seq.setDefinition('FOV', arg.fov*1e-2*ones(1,3));
    seq.setDefinition('Name', 'lps');
    seq.write([arg.dir,'/lps.seq']);
    if arg.writepge
        ceq = seq2ceq([arg.dir,'/lps.seq']);
        sysGE = pge2.opts(arg.psd_rf_wait, arg.psd_grd_wait,...
            mr.convert(arg.sys.maxB1,'Hz','mT')*10, ... % max rf amplitude (G)
            mr.convert(arg.sys.maxGrad,'Hz/m','mT/m')*0.1, ... % max gradient amplitude (G/cm)
            mr.convert(arg.sys.maxSlew,'Hz/m/s','T/m/s')*0.1, ... % max slew rate (G/cm/ms)
            arg.coil, ... % coil model
            'GRAD_UPDATE_TIME', arg.sys.gradRasterTime, ...
            'RF_UPDATE_TIME', arg.sys.rfRasterTime, ...
            'adc_raster_time', arg.sys.adcRasterTime, ...
            'adc_ringdown_time', arg.sys.adcDeadTime, ...
            'rf_dead_time', arg.sys.rfDeadTime, ...
            'rf_ringdown_time', arg.sys.rfRingdownTime, ...
            'gamma', arg.sys.gamma*1e-4 ... % GMR (Hz/G)
            );
        pge2.writeceq(ceq, [arg.dir,'/lps.pge'], ...
            'sysGE', sysGE, 'pislquant', min(arg.pislquant,1));
    end
    if exist([arg.dir,'/seq_args.h5'], 'file')
        delete([arg.dir,'/seq_args.h5']);
    end
    lpsutl.saveh5struct([arg.dir,'/seq_args.h5'], arg);
    
    % the sequence is ready, so let's see what we got
    if arg.plotseq
        figure
        seq.plot();
    end

end
