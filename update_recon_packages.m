% gets the current compatible packages for running code from this repo
% by David Frey (djfrey@umich.edu)

%% MIRT
fprintf('updating MIRT... ')
system('[ -d "./MIRT" ] && rm -rf ./MIRT');
system('git clone git@github.com:JeffFessler/MIRT.git 2> /dev/null');
fprintf('done.\n')
run MIRT/setup.m

%% BART (for sensitivity map estimation with eSPIRIT)
fprintf('updating bart... ')
system('[ -d "./bart" ] && rm -rf ./bart');
system('git clone git@github.com:mrirecon/bart.git 2> /dev/null');
system('cd ./bart && make 2> /dev/null && cd ../')
fprintf('done.\n')
run bart/startup
