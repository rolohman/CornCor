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

dt2=diff(dn(intid),[],2); %delta time in days
load baselines.txt
bpr     = baselines(id1)-baselines(id2);
abpr    = abs(bpr);
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

output.cpx=cpx;
output.abp=abpr;
output.intid=intid;
output.dn=dn;

clear cpx

cors   = squeeze(cors(:,x0,y0));
hp     = squeeze(hp(:,x0,y0));

ph     = transpose(invert_phsmat(hp,nd,intid));

diags=dt==1;
first=intid(:,1)==1;

rm=mean(ph(2:end).*conj(hp(first))); % just rotates ph relative to hp for plotting
ph=ph.*conj(rm);


%mags          = squeeze(ampsum(:,x0,y0));
mags=squeeze(ampsum);
diffmag       = mags(intid(:,2))-mags(intid(:,1));

figure('position',[38   231   648   617],'name',[num2str(xpt) ' ' num2str(ypt)])
subplot(2,2,1)
triplot(angle(hp),dn,intid)
caxis([-pi pi])
title('hp \phi')
datetick('x'),datetick('y'),axis tight

subplot(2,2,2)
triplot(cors,dn,intid)
caxis([0 1])
c2=caxis;
title('|\gamma|')
datetick('x'),datetick('y'),axis tight

subplot(4,1,3)
scatter(dt2,cors,10,angle(hp),'filled'),grid on
xlabel('days')
ylabel('|\gamma|')
title('color=hp \phi')

subplot(4,1,4)
scatter(dn(2:end),angle(hp(first)),30,cors(first),'filled'),axis([min(dn) max(dn) -pi pi]),grid on
hold on
plot(dn,angle(ph),'k.')
datetick
title('hp\phi vs. first date, color=|\gamma|')
ylabel('hp \phi')
axis tight
