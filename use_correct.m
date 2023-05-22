define_params
params.dely=1;
nx=params.nx;
ny=params.ny;
mincor = 0.45;
init_dir(params)
[dates,slcnames]                  = list_slcs(params);
[dn,nd,intid,dt,ni,id1,id2,diags] = define_pairs(dates,params.maxdt,params.dt1);
[Gi0,Gr0]                         = make_G(ni,nd,id1,id2);
filenames                         = make_filenames(params,dates,nd);
n=params.rx*params.ry;
shortid=find(dt==1);

skips=find(dt==2);

load baselines.txt
bpr     = baselines(id1)-baselines(id2);
abpr    = abs(bpr);
bigbase = abpr'>params.minbase;

[windx,windy,wind,windn,wind3,windn3] = make_kernel(params);


smallx=1:params.rx/2:nx;
smally=[params.ry*2 params.ry*2+[1:params.ry/2:params.dely params.dely+1]];
[sy,sx] = meshgrid(smally,smallx);
cs=linspace(0,1,100);
cutoffs=linspace(-0.01,0.1,20);
for j=50
%for j=1:ny
    j
    
    cpx          = load_slc_chunk(params,slcnames,1,nx,j-floor(length(windy)/2),j+floor(length(windy)/2)+params.dely,'cpx');
    
    mk             = load_slc_chunk(params,filenames.mk,1,nx,j-floor(length(windy)/2),j+floor(length(windy)/2)+params.dely,'r4');
    mk(mk>10)      = 10;
    mk(isnan(mk))  = 0;
    slopes         = load_slc_chunk(params,filenames.slope,1,nx,j-floor(length(windy)/2),j+floor(length(windy)/2)+params.dely,'r4');
    tphs_orig      = load_slc_chunk(params,{filenames.tphs_orig},1,nx,j-floor(length(windy)/2),j+floor(length(windy)/2)+params.dely,'r4');
    tphs_resk      = load_slc_chunk(params,filenames.tphs_resk,1,nx,j-floor(length(windy)/2),j+floor(length(windy)/2)+params.dely,'r4');
    
    dres = tphs_resk-tphs_orig;
     [~,ints,cors]         = make_cor(cpx,intid,wind,windn,wind3,windn3,params);
     h1 = hist(cors(:),cs);
 %figure
 %plot(cs,h1)
 %hold on
 
    disp('done loading')
    means          = atan(mk);
    
    for i=1:length(cutoffs)
    correction     = exp(1j*mk.*slopes).*conj(exp(1j*means));
    correction(repmat(dres,211,1,1)<cutoffs(i))=1;
    
    newhp          = cpx.*conj(correction); 
     
    [~,ints_new,cors_new] = make_cor(newhp,intid,wind,windn,wind3,windn3,params);
    disp('done making cor')

     h3 = hist(cors_new(:),cs);
 %    plot(cs,h3)
%     pause(0.01)
     allh3(i,:)=h3;
    end
end

save allh3 allh3

