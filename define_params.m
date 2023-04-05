params.slcdir='cropped_41000_11942_375_200/SLC_VV/';
params.outdir=[params.slcdir 'trial1/'];

params.nx     = 375;
params.ny     = 200;
params.rx     = 10;
params.ry     = 10;
params.kernel = 'Gaussian';

params.OPTIONS = optimoptions('lsqnonlin','MaxIterations',50,'Display','none');



params.bpw     = 200; %weights for baseline dependence estimation
params.minbase = 50; %only use baselines > this amount

params.dely    = 10;
params.k       = [1 2];


   