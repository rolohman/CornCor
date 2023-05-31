%params.slcdir='cropped_15420_1396_1600_320/SLC_vv/';
%params.slcdir='cropped_41000_11942_375_200/SLC_VV/';
params.slcdir='cropped_16162_1646_200_100/SLC_vv/';
params.outdir=[params.slcdir 'trial1/'];

%params.nx     = 1600;
%params.ny     = 320;
%params.nx=375;
%params.ny=200;
params.nx=200;
params.ny=100;
params.rx     = 10;
params.ry     = 10;
params.kernel = 'Gaussian';

params.OPTIONS = optimoptions('lsqnonlin','MaxIterations',20,'Display','none','FunctionTolerance',1e-3);



params.bpw     = 200; %weights for baseline dependence estimation
params.minbase = 50; %only use baselines > this amount

params.dely    = 10; %number of lines to process in each batch
params.k       = [1];
params.maxdt   = 45; %maximum time interval to use in calcs (days)
params.dt1     = 2;  %always make pairs starting with dates <= this date


   
