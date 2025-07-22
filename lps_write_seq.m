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
    arg.nspokes = 23; % number of lps spokes
    arg.t_seg = 1120; % time/segment (us)
    arg.t_rf = 16; % time/rf pulse (us)
    arg.fa = 4; % rf flip angle (deg)
    arg.plotseq = false; % option to plot the sequence
    arg.pislquant = 1; % number of TRs to use for prescan
    arg.writepge = true; % option to convert seq to pge file
    
    % parse arguments
    arg = vararg_pair(arg,varargin);
    
    % initialize sequence object
    seq = mr.Sequence(arg.sys);
    warning('OFF', 'mr:restoreShape');

    % create looping star waveforms
    [g_wav,rf_wav,rf_del] = lpsutl.gen_lps_waveforms( ...
        'sys', arg.sys, ... % pulseq mr system structure
        'fov', arg.fov, ... % fov (cm)
        'N', arg.N_nom, ... % nominal matrix size
        'nspokes', arg.nspokes, ... % number of lps spokes
        't_seg', arg.t_seg, ... % time/segment
        't_rf', arg.t_rf, ... % rf pulse width
        'fa', arg.fa ... % rf flip angle (deg)
        );
    G0 = padarray(g_wav.',[1,0],0,'post');

    % create rf (only take 1st half to play during block 1)
    rf = mr.makeArbitraryRf( ...
        rf_wav, arg.nspokes*arg.fa*pi/180, ...
        'delay', rf_del, ...
        'use', 'excitation', ...
        'system', arg.sys);

    % create ADC
    nseg = round(arg.t_seg*1e-6/arg.sys.gradRasterTime);
    acq_len = arg.sys.adcSamplesDivisor*ceil(arg.nspokes*nseg ...
        /arg.sys.adcSamplesDivisor);
    adc = mr.makeAdc(acq_len, ...
        'Duration', arg.sys.adcRasterTime*acq_len, ...
        'system', arg.sys);

    % create TR delay
    delay_time = arg.sys.blockDurationRaster*ceil((arg.tr*1e-3 - arg.sys.gradRasterTime*size(G0,2))/arg.sys.blockDurationRaster);
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
            iG = R * G0;
        else
            % rotate as if not a dummy shot
            R = lpsutl.rot_3dtga(i+arg.dummyshots+arg.pislquant+1, 1);
            iG = R * G0;
        end

        % write fID portion to sequence
        gx_fid = mr.makeArbitraryGrad('x', 0.99*iG(1,1:end/2), ...
            'system', arg.sys, ...
            'first', 0.99*iG(1,1), ...
            'last', 0.99*iG(1,end/2));
        gy_fid = mr.makeArbitraryGrad('y', 0.99*iG(2,1:end/2), ...
            'system', arg.sys, ...
            'first', 0.99*iG(2,1), ...
            'last', 0.99*iG(2,end/2));
        gz_fid = mr.makeArbitraryGrad('z', 0.99*iG(3,1:end/2), ...
            'system', arg.sys, ...
            'first', 0.99*iG(3,1), ...
            'last', 0.99*iG(3,end/2));
        seq.addBlock(rf, gx_fid, gy_fid, gz_fid, ...
            mr.makeLabel('SET','TRID', ...
            1 + isDummyTR + 2*isReceiveGainCalibrationTR));

        % write GRE portion to sequence
        gx_gre = mr.makeArbitraryGrad('x', 0.99*iG(1,end/2+1:end), ...
            'system', arg.sys, ...
            'first', 0.99*iG(1,end/2+1), ...
            'last', 0.99*iG(1,end));
        gy_gre = mr.makeArbitraryGrad('y', 0.99*iG(2,end/2+1:end), ...
            'system', arg.sys, ...
            'first', 0.99*iG(2,end/2+1), ...
            'last', 0.99*iG(2,end));
        gz_gre = mr.makeArbitraryGrad('z', 0.99*iG(3,end/2+1:end), ...
            'system', arg.sys, ...
            'first', 0.99*iG(3,end/2+1), ...
            'last', 0.99*iG(3,end));
        if isDummyTR
            seq.addBlock(gx_gre, gy_gre, gz_gre); % no adc
        else
            seq.addBlock(adc, gx_gre, gy_gre, gz_gre);
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
