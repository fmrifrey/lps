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
    arg.necho = 1; % number of gradient-recalled echoes to acquire
    arg.nspokes = 23; % number of lps spokes
    arg.psd_rf_wait = 100e-6; % RF–gradient delay (s), GE scanner-specific
    arg.psd_grd_wait = 100e-6;   % ADC–gradient delay (s), scanner-specific
    arg.coil = 'xrm'; % GE coil model
    arg.t_seg = 1120; % time/segment (us)
    arg.t_rf = 16; % time/rf pulse (us)
    arg.fa = 4; % rf flip angle (deg)
    arg.satellite = true; % option to use satellite trajectory
    arg.plotseq = false; % option to plot the sequence
    arg.pislquant = 1; % number of TRs to use for prescan
    arg.writepge = true; % option to convert seq to pge file
    
    % parse arguments
    arg = vararg_pair(arg,varargin);
    
    % initialize sequence object
    arg.sys.rfDeadTime = 0; % to start rf pulse at t=0 in fID block
    seq = mr.Sequence(arg.sys);
    warning('OFF', 'mr:restoreShape');

    % create looping star waveforms
    [g_wav,g0,rf_wav,t_ramp] = lpsutl.gen_lps_waveforms( ...
        'sys', arg.sys, ... % pulseq mr system structure
        'fov', arg.fov, ... % fov (cm)
        'N', arg.N_nom, ... % nominal matrix size
        'nspokes', arg.nspokes, ... % number of lps spokes
        't_seg', arg.t_seg, ... % time/segment
        't_rf', arg.t_rf, ... % rf pulse width
        'satellite', arg.satellite, ... % option to use satellite
        'fa', arg.fa ... % rf flip angle (deg)
        );
    
    % create gradient waveform objects for fID portion
    gx_fid = mr.makeArbitraryGrad('x', 0.99*g_wav(:,1).', ...
        'system', arg.sys, ...
        'first', 0.99*g0(1), ...
        'last', 0.99*g0(1));
    gy_fid = mr.makeArbitraryGrad('y', 0.99*g_wav(:,2).', ...
        'system', arg.sys, ...
        'first', 0.99*g0(2), ...
        'last', 0.99*g0(2));
    gz_fid = mr.makeArbitraryGrad('z', 0.99*g_wav(:,3).', ...
        'system', arg.sys, ...
        'first', 0.99*g0(3), ...
        'last', 0.99*g0(3));

    % create gradient waveform objects for GRE portion
    gx_gre = mr.makeArbitraryGrad('x', repmat(0.99*g_wav(:,1).', [1,arg.necho]), ...
        'system', arg.sys, ...
        'first', 0.99*g0(1), ...
        'last', 0.99*g0(1));
    gy_gre = mr.makeArbitraryGrad('y', repmat(0.99*g_wav(:,2).', [1,arg.necho]), ...
        'system', arg.sys, ...
        'first', 0.99*g0(2), ...
        'last', 0.99*g0(2));
    gz_gre = mr.makeArbitraryGrad('z', repmat(0.99*g_wav(:,3).', [1,arg.necho]), ...
        'system', arg.sys, ...
        'first', 0.99*g0(3), ...
        'last', 0.99*g0(3));

    % create gradient ramp objects
    nramp = ceil(t_ramp/arg.sys.gradRasterTime);
    ramp_wav = 1/nramp * ((0:nramp - 1) + 0.5);
    gx_ramp_up = mr.makeArbitraryGrad('x', ...
        0.99*g_wav(1,1)*ramp_wav, ...
        'system', arg.sys, ...
        'first', 0, ...
        'last', 0.99*g0(1));
    gx_ramp_down = mr.makeArbitraryGrad('x', ...
        0.99*g_wav(end,1)*(1-ramp_wav), ...
        'system', arg.sys, ...
        'first', 0.99*g0(1), ...
        'last', 0);
    gy_ramp_up = mr.makeArbitraryGrad('y', ...
        0.99*g_wav(1,2)*ramp_wav, ...
        'system', arg.sys, ...
        'first', 0, ...
        'last', 0.99*g0(2));
    gy_ramp_down = mr.makeArbitraryGrad('y', ...
        0.99*g_wav(end,2)*(1-ramp_wav), ...
        'system', arg.sys, ...
        'first', 0.99*g0(2), ...
        'last', 0);
    gz_ramp_up = mr.makeArbitraryGrad('z', ...
        0.99*g_wav(1,3)*ramp_wav, ...
        'system', arg.sys, ...
        'first', 0, ...
        'last', 0.99*g0(3));
    gz_ramp_down = mr.makeArbitraryGrad('z', ...
        0.99*g_wav(end,3)*(1-ramp_wav), ...
        'system', arg.sys, ...
        'first', 0.99*g0(3), ...
        'last', 0);

    % create rf waveform
    rf = mr.makeArbitraryRf( ...
        rf_wav, arg.nspokes*arg.fa*pi/180, ...
        'use', 'excitation', ...
        'system', arg.sys);

    % create ADC
    nseg = round(arg.t_seg*1e-6/arg.sys.adcRasterTime);
    acq_len = arg.sys.adcSamplesDivisor*ceil(arg.nspokes*nseg*arg.necho ...
        /arg.sys.adcSamplesDivisor);
    adc = mr.makeAdc(acq_len, ...
        'Duration', arg.sys.adcRasterTime*acq_len, ...
        'system', arg.sys);

    % create TR delay
    delay_time = arg.tr*1e-3;
    delay_time = delay_time - mr.calcDuration(gx_ramp_up);
    delay_time = delay_time - mr.calcDuration(gx_fid);
    delay_time = delay_time - mr.calcDuration(gx_gre);
    delay_time = delay_time - mr.calcDuration(gx_ramp_down);
    delay_time = arg.sys.blockDurationRaster*ceil(delay_time/arg.sys.blockDurationRaster);
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
        
        % make rotation object from rotation matrix
        rot = mr.makeRotation(R);

        if arg.satellite
            % write ramp-ups to sequence
            seq.addBlock(gx_ramp_up, gy_ramp_up, gz_ramp_up, rot, ...
                mr.makeLabel('SET','TRID', ...
                1 + isDummyTR + 2*isReceiveGainCalibrationTR));

            % write fID portion to sequence
            seq.addBlock(rf, gx_fid, gy_fid, gz_fid, rot);
    
            % write GRE portion to sequence
            if isDummyTR
                seq.addBlock(gx_gre, gy_gre, gz_gre, rot); % no adc
            else
                seq.addBlock(adc, gx_gre, gy_gre, gz_gre, rot);
            end
    
            % write ramp-downs to sequence
            seq.addBlock(gx_ramp_down, gy_ramp_down, gz_ramp_down, rot);

        else % don't include z gradient
            % write ramp-ups to sequence
            seq.addBlock(gx_ramp_up, gy_ramp_up, rot, ...
                mr.makeLabel('SET','TRID', ...
                1 + isDummyTR + 2*isReceiveGainCalibrationTR));

            % write fID portion to sequence
            seq.addBlock(rf, gx_fid, gy_fid, rot);

            % write GRE portion to sequence
            if isDummyTR
                seq.addBlock(gx_gre, gy_gre, rot); % no adc
            else
                seq.addBlock(adc, gx_gre, gy_gre, rot);
            end

            % write ramp-downs to sequence
            seq.addBlock(gx_ramp_down, gy_ramp_down, rot);
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
