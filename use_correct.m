define_params
params.dely=50;
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
x1=1;
x2=nx;

[windx,windy,wind,windn,wind3,windn3] = make_kernel(params);


smallx=1:params.rx/2:nx;
smally=[params.ry*2 params.ry*2+[1:params.ry/2:params.dely params.dely+1]];
[sy,sx] = meshgrid(smally,smallx);
cs=linspace(0,1,100);

Gt=[dn-dn(1)];
Ggt=inv(Gt'*Gt)*Gt';

for j=25
%for j=1:dely:ny
    j
       oldhp          = load_slc_chunk(params,filenames.hp,x1,x2,j-floor(length(windy)/2),j+floor(length(windy)/2)+params.dely,'r4');
 
       rates=reshape(Ggt*reshape(oldhp,nd,200*91),200,91);
       
%     cpx          = load_slc_chunk(params,slcnames,x1,x2,j-floor(length(windy)/2),j+floor(length(windy)/2)+params.dely,'cpx');
%     
%     mk             = load_slc_chunk(params,filenames.mk,x1,x2,j-floor(length(windy)/2),j+floor(length(windy)/2)+params.dely,'r4');
%     mk(mk>10)      = 10;
%     mk(isnan(mk))  = 0;
%     slopes         = load_slc_chunk(params,filenames.slope,x1,x2,j-floor(length(windy)/2),j+floor(length(windy)/2)+params.dely,'r4');
%     tphs_orig      = load_slc_chunk(params,{filenames.tphs_orig},x1,x2,j-floor(length(windy)/2),j+floor(length(windy)/2)+params.dely,'r4');
%     tphs_resk      = load_slc_chunk(params,filenames.tphs_resk,x1,x2,j-floor(length(windy)/2),j+floor(length(windy)/2)+params.dely,'r4');
%     dres = tphs_resk-tphs_orig;
%     [~,ints,cors,hp]         = make_cor(cpx,intid,wind,windn,wind3,windn3,params);
% 
%     tmp=repmat(dres,nd,1,1);
%     
%     disp('done loading')
%     means          = atan(mk);
%     
%     correction     = exp(1j*mk.*slopes);
%     correction(tmp<-0.01)=1;
%     
%     newhp          = cpx.*conj(correction);
%     
%     [~,ints_new,cors_new] = make_cor(newhp,intid,wind,windn,wind3,windn3,params);
%     disp('done making cor')
%     

end



