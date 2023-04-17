%function output=plot_pt(xpt,ypt,plotflag)
xpt=101
ypt=101
define_params
oldparams=params;
nx=params.nx;
ny=params.ny;
rx=params.rx;
ry=params.ry;
params.dely=0;
init_dir(params)
adir           = [params.outdir 'models/'];
[dates,slcnames]                  = list_slcs(params);
[dn,nd,intid,dt,ni,id1,id2,diags] = define_pairs(dates);
[Gi0,Gr0]                         = make_G(ni,nd,id1,id2);
fids                              = open_files(params,dates,nd,'r');


load baselines.txt
bpr     = baselines(id1)-baselines(id2);
abpr    = abs(bpr);
bigbase = abpr'>params.minbase;

[windx,windy,wind,~,wind3] = make_kernel(params);
windn                      = conv2(ones(size(wind)),wind,'same');
windn3(1,:,:)              = windn;


cpx    = load_slc_chunk(params,slcnames,xpt-floor(length(windx)/2),xpt+floor(length(windx)/2),ypt-floor(length(windy)/2),ypt+floor(length(windy)/2),'cpx');
amps   = cpx.*conj(cpx);
ampsum = sqrt(convn(amps,wind3,'same')./windn3);

disp('done loading')

[gamma,ints,cors,hp]=make_cor(cpx,intid,wind,windn,wind3,windn3,params);
disp('done making cor')

cors   = squeeze(cors(:,rx*2+1,ry*2+1));
hp     = squeeze(hp(:,rx*2+1,ry*2+1));
gamma  = squeeze(gamma(:,rx*2+1,ry*2+1));


triplets0       = trip_nm(gamma./cors,intid,1);
slopes          = squeeze(load_slc_chunk(params,{[adir 'slopesk' num2str(params.k(1)) '.r4']},xpt-floor(length(windx)/2),xpt+floor(length(windx)/2),ypt-floor(length(windy)/2),ypt+floor(length(windy)/2),'r4'));
correction      = ones(nd,rx*4+1,ry*4+1);
for i=1:nd
    mk                = load_slc_chunk(params,{[adir 'mk/' dates(i).name '.mk' num2str(params.k(1))]},xpt-floor(length(windx)/2),xpt+floor(length(windx)/2),ypt-floor(length(windy)/2),ypt+floor(length(windy)/2),'r4');
    mk(mk>10)         = 10;
    mk(isnan(mk))     = 0;
    mk                = squeeze(mk);
    means             = atan(mk);
    correction(i,:,:) = exp(1j*mk.*slopes).*conj(exp(1j*means));
    allmk(i)          = mk(rx*2+1,ry*2+1);
end


mags          = squeeze(ampsum(:,rx*2+1,ry*2+1));
diffmag       = mags(intid(:,2))-mags(intid(:,1));
diffmk        = allmk(intid(:,2))-allmk(intid(:,1));

slopes        = slopes(rx*2+1,ry*2+1);

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
plot(allmk,allmk.*slopes,'r.')
xlabel('m')
ylabel('hp phase')
grid on
axis([0 2 -pi pi])


[newgamma,newints,newcors,newhp]=make_cor(cpx.*(correction),intid,wind,windn,wind3,windn3,params);
newcors   = squeeze(newcors(:,rx*2+1,ry*2+1));
newhp     = squeeze(newhp(:,rx*2+1,ry*2+1));
newgamma  = squeeze(newgamma(:,rx*2+1,ry*2+1));

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
axis([0 2 -pi pi])
xlabel('m')
ylabel('hp phase')
