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
% necho - number of echoes for multi-echo looping star
% nspokes - number of looping star spokes
% t_seg - time per spoke (us)
% t_rf - rf hard pulse width (us)
% fa - flip angle (deg)
% acq_fID - option to acquire fID segments
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
    arg.dummyshots = 1; % number of dummy shots
    arg.nrep = 1; % number of rotation sequence repetitions
    arg.nint = 1; % number of interleaves (2D in-plane rotations)
    arg.nprj = 16; % number of projections (3D thru-plane rotations)
    arg.necho = 1; % number of echoes for multi-echo lps
    arg.nspokes = 23; % number of lps spokes
    arg.t_seg = 1120; % time/segment (us)
    arg.t_rf = 16; % time/rf pulse (us)
    arg.fa = 4; % rf flip angle (deg)
    arg.acq_fIDs = false; % option to acquire fID signals
    arg.plotseq = false; % option to plot the sequence
    arg.pislquant = 1; % number of TRs to use for prescan
    arg.writepge = true; % option to convert seq to pge file
    
    % parse arguments
    arg = vararg_pair(arg,varargin);
    
    % initialize sequence object
    seq = mr.Sequence(arg.sys);
    warning('OFF', 'mr:restoreShape');

    % create looping star gradient waveforms
    [gx,gy,g_amp,t_ramp] = lpsutl.gen_lps_waveforms( ...
        'sys', arg.sys, ... % pulseq mr system structure
        'fov', arg.fov, ... % fov (cm)
        'N', arg.N_nom, ... % nominal matrix size
        'nspokes', arg.nspokes, ... % number of lps spokes
        't_seg', arg.t_seg ... % time/segment
        );

    % create fID rf and ADC objects
    rf_objs = cell(arg.nspokes,1);
    adc_fID_objs = cell(arg.nspokes,1);
    gx_fID_objs = cell(arg.nspokes,1);
    gy_fID_objs = cell(arg.nspokes,1);
    fID_acq_len = round((arg.t_seg*1e-6 - arg.t_rf*1e-6 - 2*arg.sys.rfRingdownTime)/arg.sys.adcRasterTime);
    fID_acq_len = arg.sys.adcSamplesDivisor*floor(fID_acq_len/arg.sys.adcSamplesDivisor);
    for spoke = 1:arg.nspokes
        % gradient objects
        gx_fID_objs{spoke} = mr.makeArbitraryGrad('x', ...
            round(0.99*gx(:,spoke)).', ...
            'system', arg.sys, ...
            'first', round(0.99*g_amp*cos(2*pi*(spoke-1)/arg.nspokes)), ...
            'last', round(0.99*g_amp*cos(2*pi*spoke/arg.nspokes)));
        gy_fID_objs{spoke} = mr.makeArbitraryGrad('y', ...
            round(0.99*gy(:,spoke)).', ...
            'system', arg.sys, ...
            'first', round(0.99*g_amp*sin(2*pi*(spoke-1)/arg.nspokes)), ...
            'last', round(0.99*g_amp*sin(2*pi*spoke/arg.nspokes)));

        % rf object
        rf_objs{spoke} = mr.makeBlockPulse( ...
            arg.fa*pi/180, ...
            'Duration', arg.t_rf*1e-6, ...
            'use', 'excitation', ...
            'Delay', 0, ...
            'system', arg.sys);

        % adc object
        adc_fID_objs{spoke} = mr.makeAdc( ...
            fID_acq_len, ...
            'Delay', arg.t_rf*1e-6 + arg.sys.rfRingdownTime, ...
            'Duration', arg.sys.adcRasterTime*fID_acq_len, ...
            'system', arg.sys);
    end

    % create gradient waveform objects for the GRE block
    gx_GRE_obj = mr.makeArbitraryGrad('x', ...
        round(0.99*repmat(gx(:).',[1,arg.necho])), ...
        'system', arg.sys, ...
        'first', round(0.99*g_amp), ...
        'last', round(0.99*g_amp));
    gy_GRE_obj = mr.makeArbitraryGrad('y', ...
        round(0.99*repmat(gy(:).',[1,arg.necho])), ...
        'system', arg.sys, ...
        'first', 0, ...
        'last', 0);

    % create GRE ADC objects
    GRE_acq_len = arg.necho*length(gx(:));
    adc_GRE_obj = mr.makeAdc( ...
        GRE_acq_len, ...
        'Delay', 0, ...
        'Duration', arg.sys.adcRasterTime*GRE_acq_len, ...
        'system', arg.sys);

    % create ramp objects
    nramp = ceil(t_ramp/arg.sys.gradRasterTime);
    ramp_wav = 1/nramp * ((0:nramp - 1) + 0.5);
    gx_ramp_up_obj = mr.makeArbitraryGrad('x', ...
        0.99*g_amp*ramp_wav, ...
        'system', arg.sys, ...
        'first', 0, ...
        'last', 0.99*g_amp);
    gx_ramp_down_obj = mr.makeArbitraryGrad('x', ...
        0.99*g_amp*(1-ramp_wav), ...
        'system', arg.sys, ...
        'first', 0.99*g_amp, ...
        'last', 0);

    % create TR delay
    tr_min = 2*t_ramp + (arg.necho+1)*(arg.t_seg*1e-6)*arg.nspokes;
    if arg.tr*1e-3 < tr_min
        error('specified tr (%.3f ms) < min tr (%.3f ms)', arg.tr, tr_min*1e3);
    end
    delay_time = arg.sys.blockDurationRaster*ceil((arg.tr*1e-3 - tr_min)/arg.sys.blockDurationRaster);
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
        rot_obj = mr.makeRotation(R);

        % make label object for current TR
        lbl_obj = mr.makeLabel('SET','TRID', ...
            1 + isDummyTR + 2*isReceiveGainCalibrationTR);

        % write ramp up to sequence
        seq.addBlock( ...
            gx_ramp_up_obj, ... % gradient ramp-up
            lbl_obj, ... % label object
            rot_obj ... % rotation object
            );

        % loop through spokes
        for spoke = 1:arg.nspokes
            % add fID block to sequence
            if isDummyTR || ~arg.acq_fIDs % no adc
                seq.addBlock( ...
                    rf_objs{spoke}, ... % rf pulses
                    gx_fID_objs{spoke}, gy_fID_objs{spoke}, ... % gradients
                    rot_obj ... % rotation object
                    );
            else
                seq.addBlock( ...
                    rf_objs{spoke}, ... % rf pulses
                    adc_fID_objs{spoke}, ... % fID adc objects
                    gx_fID_objs{spoke}, gy_fID_objs{spoke}, ... % gradients
                    rot_obj ... % rotation object
                    );
            end
        end

        % write GRE block to sequence
        if isDummyTR % no adc
            seq.addBlock( ...
                gx_GRE_obj, gy_GRE_obj, ... % gradients
                rot_obj ... % rotation object
                );
        else
            seq.addBlock( ...
                adc_GRE_obj, ... % echo adc object
                gx_GRE_obj, gy_GRE_obj, ... % gradients
                rot_obj ... % rotation object
                );
        end

        % write ramp down to sequence
        seq.addBlock( ...
            gx_ramp_down_obj, ... % ramp-down object
            rot_obj ... % rotation object
            );

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
        writeceq(ceq, [arg.dir,'/lps.pge'], 'pislquant', min(arg.pislquant,1));
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
