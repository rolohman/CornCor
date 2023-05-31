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

%amps   = cpx.*conj(cpx);
%ampsum = sqrt(convn(amps,wind3,'same')./windn3);
ampsum=cpx(:,x0,y0).*conj(cpx(:,x0,y0));
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
mincor   = 0.45;
cp0(cp0<mincor)=NaN;
diags=dt==1;
first=intid(:,1)==1;
%mags          = squeeze(ampsum(:,x0,y0));
mags=squeeze(ampsum);
diffmag       = mags(intid(:,2))-mags(intid(:,1));
diffmk        = allmk(intid(:,2))-allmk(intid(:,1));

synthhps  = exp(1j*(diffmk'*slope-atan(diffmk')));
reshp=hp.*conj(synthhps);
synthcorperm = c0.*exp(Gi0*log(cp0)');
synthcor = log(1./sqrt(diffmk'.^2+1))+log(c0)+Gi0*log(cp0)';

res      = log(cors)-synthcor;
good     = res<0.2 & synthcor>log(mincor) & cors>log(mincor);
synthcor = exp(synthcor);
synthhp  = exp(1j*hp0(intid(:,2))).*conj(exp(1j*hp0(intid(:,1))));
reshp0=angle(exp(1j*(hp0-allmk*slope+atan(allmk))));
figure('position',[1500 70 1276 911])
subplot(3,4,1)
triplot(angle(hp),dn,intid)
caxis([-pi pi])
title('hp')

subplot(3,4,2)
triplot(cors,dn,intid)
caxis([0 1])
c2=caxis;
title('cor')

subplot(3,4,3)
triplot(synthcorperm,dn,intid)
caxis(c2)
title('c0+perm')

subplot(3,4,4)
triplot(diffmk,dn,intid)
title('diffmk')

subplot(3,4,5)
triplot(angle(synthhp),dn,intid)
caxis([-pi pi])
title('hp0 fit')

subplot(3,4,6)
triplot(angle(reshp),dn,intid)
caxis([-pi pi])
title('linefit res')

subplot(3,4,7)
triplot(synthcor(good),dn,intid(good,:))
caxis(c2)
title(['res=' num2str(norm(cors-synthcor))])

subplot(3,4,8)
plot(abs(diffmk(1,:))',angle(hp).*sign(diffmk(1,:))','k.','markersize',1)
hold on
scatter(allmk,hp0,36,mags,'filled','MarkerEdgeColor','k')
plot(allmk,angle(exp(1j*allmk.*slope))-atan(allmk),'r.')
xlabel('m')
ylabel('hp phase')
grid on
%plot([0.4 0.4],[-pi pi],'k:')
%axis([0 2.5 -pi pi])
ax=axis;
axis([ax(1:2) -pi pi])

subplot(3,4,9)
scatter(dn(2:end),angle(hp(diags)),30,cors(diags),'filled'),axis([min(dn) max(dn) -pi pi]),grid on

subplot(3,4,10)
scatter(dn(2:end),angle(reshp(diags)),30,cors(diags),'filled'),axis([min(dn) max(dn) -pi pi]),grid on


subplot(3,4,11)
scatter(dn(2:end),angle(hp(first)),30,cors(first),'filled'),axis([min(dn) max(dn) -pi pi]),grid on
hold on
plot(dn,hp0,'k.')

subplot(3,4,12)
scatter(dn(2:end),angle(reshp(first)),30,cors(first),'filled'),axis([min(dn) max(dn) -pi pi]),grid on
hold on
plot(dn,reshp0,'k.')

