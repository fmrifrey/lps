function sys = get_sys_defaults(vendor)
% returns pulseq system structure with defaults for given vendor
% by David Frey (djfrey@umich.edu)
%
% inputs:
% vendor - 'ge' or 'siemens'
%
% outputs:
% sys - pulseq system structure
%

	switch lower(vendor)
		case 'ge'
			sys = mr.opts('MaxGrad', 50, 'GradUnit', 'mT/m', ...
				'MaxSlew', 170, 'SlewUnit', 'mT/m/ms', ...
				'rfDeadTime', 100e-6, ...
				'rfRingdownTime', 60e-6, ...
				'adcRasterTime', 4e-6, ...
				'gradRasterTime', 4e-6, ...
				'rfRasterTime', 2e-6, ...
				'blockDurationRaster', 4e-6, ...
				'B0', 3, ...
				'adcDeadTime', 0e-6);
		case 'siemens'
			sys = mr.opts; % pulseq already defaults to siemens
		otherwise
			error('invalid vendor');
	end

end
