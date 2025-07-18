function gre3d_write_seq(varargin)
% creates the pulseq files for 3d cartesian GRE readout
%
% by David Frey (djfrey@umich.edu)
%
% inputs:
% sys - pulseq mr system structure
% dir - output directory name
% dt - ADC sampling rate (s)
% te - TE (ms)
% tr - TR (ms)
% fov - fov (cm)
% slabfrac - excitation slab width (fraction of fov)
% fa - flip angle (deg)
% rfspoil - RF spoiling option
% fatsat - fat saturation option
% fatChemShift - fat chemical shift (ppm)
% N - 3D matrix size
% Nacs - number of ACS lines (non-accelerated at center of kspace)
% Ry - Ky acceleration factor
% Rz - Kz acceleration factor
% delta - CAIPI odd/even shift
% peorder - phase encode ordering ('snake' or 'spout')
% dummyshots - number of dummy shots (disdaqs) to play
% writepge - option to write pge file
% pislquant - number of prescan acquisitions
% plotseq - option to plot the sequence
%
% output files:
% gre3d.seq file - seq file for pulseq
% gre3d.pge file - pge file for pge2 interpreter
% seq_args.h5 - .h5 file containing copy of input arguments
%

    % set default arguments
    arg.sys = lpsutl.get_sys_defaults('ge'); % pulseq mr system structure
    arg.dir = pwd; % output destination for files
    arg.te = 'min'; % TE (ms)
    arg.tr = 30; % TR (ms)
    arg.fov = 16; % fov (cm)
    arg.slabfrac = 0.7; % excitation slab width (fraction of fov)
    arg.fa = 6; % flip angle (deg)
    arg.rfspoil = true; % RF spoiling option
    arg.fatsat = false; % fat saturation option
    arg.fatChemShift = 3.5; % fat chemical shift (ppm)
    arg.N = 128; % 3D matrix size
    arg.Nacs = 32; % width of fully sampled (ACS) region at center of kspace
    arg.Ry = 2; % Ky acceleration factor (outside ACS region)
    arg.Rz = 2; % Kz acceleration factor (outside ACS region)
    arg.delta = 1; % CAIPI odd/even shift
    arg.peorder = 'snake'; % pe ordering scheme
    arg.dt = 20e-6; % ADC sampling rate (s)
    arg.dummyshots = 100; % number of dummy shots (disdaqs) to play
    arg.writepge = true;
    arg.pislquant = 10;
    arg.plotseq = false; % option to plot the sequence

    % parse arguments
    arg = vararg_pair(arg,varargin);
    
    % get phase encode indicies
    if strcmpi(arg.peorder,'snake')
        pe_idcs = gre3dutl.snake_caipi_idcs(arg.N, arg.Ry, arg.Rz, arg.delta, arg.Nacs);
    elseif strcmpi(arg.peorder,'spout')
        pe_idcs = gre3dutl.spout_caipi_idcs(arg.N, arg.Ry, arg.Rz, arg.delta, arg.Nacs);
    else
        error('invalid option for peorder');
    end
    npe = length(pe_idcs);
    
    % Create a new sequence object
    seq = mr.Sequence(arg.sys);

    % create fat excitation pulse
    fatOffres = arg.sys.gamma*arg.sys.B0*arg.fatChemShift*1e-6;  % (Hz)
    fatBW = 200; % (Hz)
    fatPW = 12e-3; % (s)
    rf_fat = mr.makeGaussPulse(pi/2, ...
        'Duration', fatPW, ...
        'freqOffset', -fatOffres, ...
        'apodization', 0.42, ...
        'use', 'saturation', ...
        'timeBwProduct', fatBW*fatPW, ...
        'system', arg.sys);
    
    % create excitation
    [rf, gz] = mr.makeSincPulse(arg.fa*pi/180, ...
        'Duration', 3e-3, ...
        'SliceThickness', arg.slabfrac*arg.fov*1e-2, ...
        'apodization', 0.42, ...
        'use', 'excitation', ...
        'timeBwProduct', 4, ...
        'system', arg.sys);
    gzReph = mr.makeTrapezoid('z', ...
        'Area', -gz.area/2, ...
        'system', arg.sys);
    
    % create frequency encode gradient and ADC
    gx = mr.makeTrapezoid('x', ...
        'FlatArea', arg.N/(arg.fov*1e-2), ...
        'FlatTime', arg.N*arg.dt, ...
        'system', arg.sys);
    adc = mr.makeAdc(arg.N, ...
        'Duration', gx.flatTime, ...
        'Delay', gx.riseTime, ...
        'system', arg.sys);
    
    % create phase encode gradients
    gxPre = mr.makeTrapezoid('x', ...
        'Area', -gx.area/2, ...
        'system', arg.sys);
    gyPre = mr.makeTrapezoid('y', ...
        'Area', arg.N/(arg.fov*1e-2)/2, ...
        'Duration', mr.calcDuration(gxPre), ...
        'system', arg.sys);
    gzPre = mr.makeTrapezoid('z', ...
        'Area', arg.N/(arg.fov*1e-2)/2, ...
        'Duration', mr.calcDuration(gxPre), ...
        'system', arg.sys);
    
    % create spoilers
    gxSpoil = mr.makeTrapezoid('x', ...
        'Area', 4*arg.N/(arg.fov*1e-2), ...
        'system', arg.sys);
    gySpoil = mr.makeTrapezoid('y', ...
        'Area', 4*arg.N/(arg.fov*1e-2), ...
        'system', arg.sys);
    gzSpoil = mr.makeTrapezoid('z', ...
        'Area', 4*arg.N/(arg.fov*1e-2), ...
        'system', arg.sys);
    
    % determine echo time delay
    te_min = mr.calcDuration(rf)/2 + mr.calcDuration(gzReph) + ...
        mr.calcDuration(gzPre) + mr.calcDuration(gx)/2;
    te_min = arg.sys.gradRasterTime*ceil(te_min/arg.sys.gradRasterTime);
    if strcmpi(arg.te,'min')
        arg.te = te_min*1e3;
        te_delay = 0;
    elseif arg.te*1e-3 >= te_min
        te_delay = arg.te*1e-3 - te_min;
    else
        error('echo time must be >= %.3fms', te_min*1e3);
    end
    
    % determine repetition time delay
    tr_min = arg.fatsat*mr.calcDuration(rf_fat) + mr.calcDuration(gxSpoil) + ...
        mr.calcDuration(rf) + mr.calcDuration(gzReph) + ...
        mr.calcDuration(gzPre) + te_delay + mr.calcDuration(gx) + ...
        mr.calcDuration(gzPre);
    tr_min = arg.sys.gradRasterTime*ceil(tr_min / arg.sys.gradRasterTime);
    if strcmpi(arg.tr,'min')
        arg.tr = 1e3*tr_min;
        tr_delay = 0;
    elseif 1e-3*arg.tr >= tr_min
        tr_delay = 1e-3*arg.tr - tr_min;
    else
        error('repetition time must be >= %.3fms', tr_min*1e3);
    end
    
    % initialize rf phase
    rf_phase = 0;
    rf_inc = 0;
    
    % loop through shots
    for i = (-arg.dummyshots-arg.pislquant+1):npe
        isDummyTR = i <= -arg.pislquant;
        isReceiveGainCalibrationTR = i < 1 & i > -arg.pislquant;
        lbl = mr.makeLabel('SET', 'TRID', 1 + isDummyTR + 2*isReceiveGainCalibrationTR);

        % calculate RF and ADC phase
        rf.phaseOffset = rf_phase/180*pi;
        rf_fat.phaseOffset = rf_phase/180*pi;
        adc.phaseOffset = rf_phase/180*pi;
        rf_inc = mod(rf_inc + arg.rfspoil*117, 360.0);
        rf_phase = mod(rf_phase + rf_inc, 360.0);

        if arg.fatsat
            % add fat saturation
            seq.addBlock(rf_fat, lbl);
            seq.addBlock(gxSpoil, gySpoil, gzSpoil);

            % add excitation and refocuser
            seq.addBlock(rf, gz);
            seq.addBlock(gzReph);
        else
            % add spoiler
            seq.addBlock(gxSpoil, gySpoil, gzSpoil, lbl)

            % add excitation and refocuser
            seq.addBlock(rf, gz);
            seq.addBlock(gzReph);
        end
    
        % add te delay
        if te_delay > 0
            seq.addBlock(mr.makeDelay(te_delay));
        end
    
        % add phase encode gradients
        if i > 0 && ~isDummyTR
            iy = pe_idcs(i,1);
            iz = pe_idcs(i,2);
            pescy = -1 + 2*(iy-1)/arg.N;
            pescy = pescy + (pescy == 0)*eps;
            pescz = -1 + 2*(iz-1)/arg.N;
            pescz = pescz + (pescz == 0)*eps;
        else
            pescy = eps;
            pescz = eps;
        end
    
        % add phase encodes
        seq.addBlock(gxPre, ...
            mr.scaleGrad(gyPre, pescy), ...
            mr.scaleGrad(gzPre, pescz));
    
        % add readout
        if isDummyTR
            seq.addBlock(gx);
        else
            seq.addBlock(gx, adc);
        end
    
        % add rephasing and TR delay
        seq.addBlock(mr.scaleGrad(gyPre, -pescy), ...
            mr.scaleGrad(gzPre, -pescz));
        if tr_delay > 0
            seq.addBlock(mr.makeDelay(tr_delay));
        end
    
    end
    
    % check timing
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
    seq.setDefinition('Name', 'gre3d');
    seq.write([arg.dir,'/gre3d.seq']);
    if arg.writepge
        ceq = seq2ceq([arg.dir,'/gre3d.seq']);
        writeceq(ceq, [arg.dir,'/gre3d.pge'], 'pislquant', arg.pislquant);
    end
    if exist([arg.dir,'/seq_args.h5'], 'file')
        delete([arg.dir,'/seq_args.h5']);
    end
    lpsutl.saveh5struct([arg.dir,'/seq_args.h5'], arg);
    
    % plot the sequence
    if arg.plotseq
        seq.plot();
    end

end
