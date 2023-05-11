%function output=plot_pt(xpt,ypt,plotflag)
clear windn3 
disp('hi')
%xpt=252;
%ypt=106;
%xpt=361;
%ypt=133;
%xpt=165;
%ypt=80;
%xpt=182;ypt=54;

define_params
params.slcdir='cropped_41000_11942_375_200/SLC_VV/';
params.outdir=[params.slcdir 'trial1/'];

nx=params.nx;
ny=params.ny;
rx=params.rx;
ry=params.ry;

init_dir(params)

[dates,slcnames]                  = list_slcs(params);
params.dely = 0;
params.dt1  = length(dates);

[dn,nd,intid,dt,ni,id1,id2,diags] = define_pairs(dates,params.maxdt,params.dt1);
[Gi0,Gr0]                         = make_G(ni,nd,id1,id2);
filenames                         = make_filenames(params,dates,nd);
%fids                              = open_files(filenames,'r');


load baselines.txt
bpr     = baselines(id1)-baselines(id2);
abpr    = abs(bpr);
bigbase = abpr'>params.minbase;

[windx,windy,wind,~,wind3] = make_kernel(params);

x1 = max(1,xpt-floor(length(windx)/2));
x2 = min(nx,xpt+floor(length(windx)/2));
y1 = max(1,ypt-floor(length(windy)/2));
y2 = min(ny,ypt+floor(length(windy)/2));
x0 = xpt-x1+1;
y0 = ypt-y1+1;
dx = x2-x1+1;
dy = y2-y1+1;
cpx    = load_slc_chunk(params,slcnames,x1,x2,y1,y2,'cpx');

windn         = conv2(ones(size(cpx,2),size(cpx,3)),wind,'same');
windn3(1,:,:) = windn;

amps   = cpx.*conj(cpx);
ampsum = sqrt(convn(amps,wind3,'same')./windn3);

disp('done loading')

[gamma,ints,cors,hp]=make_cor(cpx,intid,wind,windn,wind3,windn3,params);
disp('done making cor')

cors   = squeeze(cors(:,x0,y0));
hp     = squeeze(hp(:,x0,y0));
gamma  = squeeze(gamma(:,x0,y0));


triplets0       = trip_nm(gamma./cors,intid,1);
slopes          = squeeze(load_slc_chunk(params,{filenames.slope{1}},x1,x2,y1,y2,'r4'));
correction      = ones(nd,dx,dy);
for i=1:nd
    mk                = load_slc_chunk(params,{filenames.mk{i,1}},x1,x2,y1,y2,'r4');
    mk(mk>10)         = 10;
    mk(isnan(mk))     = 0;
    mk                = squeeze(mk);
    means             = atan(mk);
    %correction(i,:,:) = exp(1j*mk.*slopes).*conj(exp(1j*means));
   correction(i,:,:) = exp(1j*mk.*slopes);
    allmk(i)          = mk(x0,y0);
end

for i=1:nd
     hp0(i)                = load_slc_chunk(params,{filenames.hp{i}},xpt,xpt,ypt,ypt,'r4');
end
    
mags          = squeeze(ampsum(:,x0,y0));
diffmag       = mags(intid(:,2))-mags(intid(:,1));
diffmk        = allmk(intid(:,2))-allmk(intid(:,1));

slopes        = slopes(x0,y0)

figure
subplot(3,4,1)
triplot(angle(hp),dn,intid)
c1=caxis;
title('hp')

subplot(3,4,2)
triplot(cors,dn,intid)
c2=caxis;
title('cor')

subplot(3,4,3)
triplot(angle(triplets0),dn,intid)
c3=caxis;
title('phase closure')

subplot(3,4,4)
scatter(abs(diffmk(1,:))',[angle(hp.*exp(1j*atan(diffmk(1,:))'))].*sign(diffmk(1,:))',18,abs(diffmag),'filled')
hold on
plot(allmk,hp0+atan(allmk),'ko','markerfacecolor','k')
plot(allmk,angle(exp(1j*allmk.*slopes)),'r.')
xlabel('m')
ylabel('hp phase')
grid on
plot([0.4 0.4],[-pi pi],'k:')
axis([0 2.5 -pi pi])


[newgamma,newints,newcors,newhp]=make_cor(cpx.*conj(correction),intid,wind,windn,wind3,windn3,params);
newcors   = squeeze(newcors(:,x0,y0));
newhp     = squeeze(newhp(:,x0,y0));
newgamma  = squeeze(newgamma(:,x0,y0));

triplets1 = trip_nm(newgamma./newcors,intid,1);

subplot(3,4,5)
triplot(angle(newhp),dn,intid)
caxis(c1)
title('hp')
subplot(3,4,6)
triplot(abs(newgamma),dn,intid)
caxis(c2)
title('cor')
subplot(3,4,7)
triplot(angle(triplets1),dn,intid)
title('phase closure')
caxis(c3)
subplot(3,4,8)
plot(abs(diffmk'),angle(newhp).*sign(diffmk'),'.')
grid on
axis([0 2.5 -pi pi])
xlabel('m')
ylabel('hp phase')
