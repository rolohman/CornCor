params.outdir='cropped_38370_11488_4800_800/SLC_VV/trial1/';
params.slcdir='cropped_38370_11488_4800_800/SLC_VV/';

params.nx     = 4800;
params.ny     = 800;
params.rx     = 10;
params.ry     = 10;
params.kernel = 'Gaussian';

params.OPTIONS = optimoptions('lsqnonlin','MaxIterations',50,'Display','none');


params.bigbase = abpr'>50;
params.bpw     = 200; %weights for baseline estimation

params.dely = 2;

   