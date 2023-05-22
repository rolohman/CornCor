function output=plot_pt1(xpt,ypt,plotflag)
clear windn3 

define_params
nx=params.nx;
ny=params.ny;
rx=params.rx;
ry=params.ry;

init_dir(params)

[dates,slcnames]                  = list_slcs(params);
params.dely                       = 0;
params.dt1                        = length(dates);

[dn,nd,intid,dt,ni,id1,id2,~]     = define_pairs(dates,params.maxdt,params.dt1);
[Gi0,~]                           = make_G(ni,nd,id1,id2);
filenames                         = make_filenames(params,dates,nd);

% load baselines.txt
% bpr     = baselines(id1)-baselines(id2);
% abpr    = abs(bpr);
% bigbase = abpr'>params.minbase;

[windx,windy,wind,~,wind3] = make_kernel(params);

x1 = max(1,xpt-floor(length(windx)/2));x2 = min(nx,xpt+floor(length(windx)/2));
y1 = max(1,ypt-floor(length(windy)/2));y2 = min(ny,ypt+floor(length(windy)/2));
x0 = xpt-x1+1;
y0 = ypt-y1+1;

cpx    = load_slc_chunk(params,slcnames,x1,x2,y1,y2,'cpx');

windn         = conv2(ones(size(cpx,2),size(cpx,3)),wind,'same');
windn3(1,:,:) = windn;

amps   = cpx.*conj(cpx);
ampsum = sqrt(convn(amps,wind3,'same')./windn3);

disp('done loading')

[~,~,cors,hp] = make_cor(cpx,intid,wind,windn,wind3,windn3,params);
disp('done making cor')

cors   = squeeze(cors(:,x0,y0));
hp     = squeeze(hp(:,x0,y0));

slope           = squeeze(load_slc_chunk(params,{filenames.slope{1}},xpt,xpt,ypt,ypt,'r4'));
c0              = squeeze(load_slc_chunk(params,{filenames.c0},xpt,xpt,ypt,ypt,'r4'));

for i=1:nd
    mk                = load_slc_chunk(params,{filenames.mk{i,1}},x1,x2,y1,y2,'r4');
    mk(mk>10)         = 10;
    mk(isnan(mk))     = 0;
    mk                = squeeze(mk);
    allmk(i)          = mk(x0,y0);
end

for i=1:nd
     hp0(i)                = load_slc_chunk(params,{filenames.hp{i}},xpt,xpt,ypt,ypt,'r4');
end
for i=1:nd-1
    cp0(i)                = load_slc_chunk(params,{filenames.perm{i}},xpt,xpt,ypt,ypt,'r4');
end


mags          = squeeze(ampsum(:,x0,y0));
diffmag       = mags(intid(:,2))-mags(intid(:,1));
diffmk        = allmk(intid(:,2))-allmk(intid(:,1));

synthcorperm = c0.*exp(Gi0*log(cp0)');
synthcor = log(1./sqrt(diffmk'.^2+1))+log(c0)+Gi0*log(cp0)';
mincor   = 0.45;
res      = log(cors)-synthcor;
good     = res<0.2 & synthcor>log(mincor) & cors>log(mincor);
synthcor = exp(synthcor);
synthhp  = exp(1j*hp0(intid(:,2))).*conj(exp(1j*hp0(intid(:,1))));

figure('position',[1500 70 1276 911])
subplot(2,3,1)
triplot(angle(hp),dn,intid)
caxis([-pi pi])
title('hp')

subplot(2,3,2)
triplot(cors,dn,intid)
c2=caxis;
title('cor')

subplot(2,3,3)
triplot(synthcorperm,dn,intid)
caxis(c2)
title('c0+perm')

subplot(2,3,4)
triplot(angle(synthhp),dn,intid)
caxis([-pi pi])
title('reshp')

subplot(2,3,5)
triplot(synthcor(good),dn,intid(good,:))
caxis(c2)
title(['res=' num2str(norm(cors-synthcor))])

subplot(2,3,6)
scatter(abs(diffmk(1,:))',angle(hp).*sign(diffmk(1,:))',18,abs(diffmag),'filled')
hold on
plot(allmk,hp0,'ko','markerfacecolor','k')
plot(allmk,angle(exp(1j*allmk.*slope))-atan(allmk),'r.')
xlabel('m')
ylabel('hp phase')
grid on
%plot([0.4 0.4],[-pi pi],'k:')
%axis([0 2.5 -pi pi])
ax=axis;
axis([ax(1:2) -pi pi])
