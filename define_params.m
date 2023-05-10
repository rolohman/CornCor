params.slcdir='cropped_41000_11942_375_200/SLC_VV/';
params.outdir=[params.slcdir 'trial3/'];

params.nx     = 375;
params.ny     = 200;
params.rx     = 10;
params.ry     = 10;
params.kernel = 'Gaussian';

params.OPTIONS = optimoptions('lsqnonlin','MaxIterations',50,'Display','none','FunctionTolerance',1e-2);



params.bpw     = 200; %weights for baseline dependence estimation
params.minbase = 50; %only use baselines > this amount

params.dely    = 10; %number of lines to process in each batch
params.k       = [1];
params.maxdt   = 10; %maximum time interval to use in calcs
params.dt1     = 3;  %always make pairs starting with dates <= this date


   
