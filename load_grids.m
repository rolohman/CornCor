%function output=plot_pt(xpt,ypt,plotflag)
clear windn3 
%xpt=101
%ypt=105
define_params
%params.slcdir='cropped_41000_11942_375_200/SLC_VV/';
%params.outdir=[params.slcdir 'trial1/'];

nx=params.nx;
ny=params.ny;
rx=params.rx;
ry=params.ry;

init_dir(params)

[dates,slcnames]                  = list_slcs(params);
nd=length(dates);
ks=length(params.k);

filenames                         = make_filenames(params,dates,length(dates));
fids                              = open_files(filenames,'r');


bpslope   = fread(fids.bpslope,[nx,ny],'real*4');
c0        = fread(fids.c0,[nx,ny],'real*4');
mcor_orig = fread(fids.mcor_orig,[nx,ny],'real*4');
mcor_res  = fread(fids.mcor_res,[nx,ny],'real*4');

tphs_orig = fread(fids.tphs_orig,[nx,ny],'real*4');
tphs_res  = fread(fids.tphs_res,[nx,ny],'real*4');

for i=1:nd-1
    perm(i,:,:)=fread(fids.perm(i),[nx,ny],'real*4');
end
for i=1:nd
    hp(i,:,:)=fread(fids.hp(i),[nx,ny],'real*4');
    for j=1:ks
        mk(i,j,:,:)=fread(fids.mk(i,j),[nx,ny],'real*4');
    end
end
for j=1:ks
    tphs_resk(j,:,:)  = fread(fids.tphs_resk(j),[nx,ny],'real*4');
    slope(j,:,:)     = fread(fids.slope(j),[nx,ny],'real*4');
end

%while we only have one k
tphs_resk=squeeze(tphs_resk);
slope=squeeze(slope);
mk=squeeze(mk);

%[i,j]=find(tphs_res>0.95 & 

return

params.slcdir='cropped_41000_11942_375_200/SLC_VV/';
params.outdir=[params.slcdir 'trial2/'];

nx=params.nx;
ny=params.ny;
rx=params.rx;
ry=params.ry;

init_dir(params)

[dates,slcnames]                  = list_slcs(params);
nd=length(dates);
ks=length(params.k);

filenames                         = make_filenames(params,dates,length(dates));
fids                              = open_files(filenames,'r');


bpslope2   = fread(fids.bpslope,[nx,ny],'real*4');
c02        = fread(fids.c0,[nx,ny],'real*4');
mcor_orig2 = fread(fids.mcor_orig,[nx,ny],'real*4');
mcor_res2  = fread(fids.mcor_res,[nx,ny],'real*4');

tphs_orig2 = fread(fids.tphs_orig,[nx,ny],'real*4');
tphs_res2  = fread(fids.tphs_res,[nx,ny],'real*4');

for i=1:nd-1
    perm2(i,:,:)=fread(fids.perm(i),[nx,ny],'real*4');
end
for i=1:nd
    hp2(i,:,:)=fread(fids.hp(i),[nx,ny],'real*4');
    for j=1:ks
        mk2(i,j,:,:)=fread(fids.mk(i,j),[nx,ny],'real*4');
    end
end
for j=1:ks
    tphs_resk2(j,:,:)  = fread(fids.tphs_resk(j),[nx,ny],'real*4');
    slope2(j,:,:)     = fread(fids.slope(j),[nx,ny],'real*4');
end

%while we only have one k
tphs_resk2=squeeze(tphs_resk2);
slope2=squeeze(slope2);
mk2=squeeze(mk2);


params.slcdir='cropped_41000_11942_375_200/SLC_VV/';
params.outdir=[params.slcdir 'trial3/'];

nx=params.nx;
ny=params.ny;
rx=params.rx;
ry=params.ry;

init_dir(params)

[dates,slcnames]                  = list_slcs(params);
nd=length(dates);
ks=length(params.k);

filenames                         = make_filenames(params,dates,length(dates));
fids                              = open_files(filenames,'r');


bpslope3   = fread(fids.bpslope,[nx,ny],'real*4');
c03        = fread(fids.c0,[nx,ny],'real*4');
mcor_orig3 = fread(fids.mcor_orig,[nx,ny],'real*4');
mcor_res3  = fread(fids.mcor_res,[nx,ny],'real*4');

tphs_orig3 = fread(fids.tphs_orig,[nx,ny],'real*4');
tphs_res3  = fread(fids.tphs_res,[nx,ny],'real*4');

for i=1:nd-1
    perm3(i,:,:)=fread(fids.perm(i),[nx,ny],'real*4');
end
for i=1:nd
    hp3(i,:,:)=fread(fids.hp(i),[nx,ny],'real*4');
    for j=1:ks
        mk3(i,j,:,:)=fread(fids.mk(i,j),[nx,ny],'real*4');
    end
end
for j=1:ks
    tphs_resk3(j,:,:)  = fread(fids.tphs_resk(j),[nx,ny],'real*4');
    slope3(j,:,:)     = fread(fids.slope(j),[nx,ny],'real*4');
end

%while we only have one k
tphs_resk3=squeeze(tphs_resk3);
slope3=squeeze(slope3);
mk3=squeeze(mk3);

return

figure
subplot(2,2,1)
plot(mcor_orig(:),mcor_res(:),'.')
subplot(2,2,2)
plot(tphs_orig(:),tphs_res(:),'.')

subplot(2,2,3)
hist(slope(:),100)

subplot(2,2,4)
plot(slope(:),tphs_orig(:)-tphs_resk(:),'.')


