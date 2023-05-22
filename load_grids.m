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


%bpslope   = fread(fids.bpslope,[nx,ny],'real*4');
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
dres=tphs_resk-tphs_orig;

figure
subplot(2,2,1)
scatter(mcor_orig(:),mcor_res(:),5,c0(:),'filled')
grid on,box on,colorbar
xlabel('mcor orig'),ylabel('mcor res')
subplot(2,2,2)
scatter(tphs_orig(:),tphs_res(:),5,slope(:),'filled')
grid on,box on,colorbar
xlabel('tphs orig'),ylabel('tphs res')
title('slope')

subplot(2,2,3)
hist(slope(:),100)
xlabel('slopes')

subplot(2,2,4)
scatter(slope(:),dres(:),5,tphs_res(:),'filled')
grid on,box on,colorbar
xlabel('slope'),ylabel('d tres')
title('tphs res')

