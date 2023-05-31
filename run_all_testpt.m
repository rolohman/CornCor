define_params
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
params.dely=0;
load baselines.txt
bpr     = baselines(id1)-baselines(id2);
abpr    = abs(bpr);
bigbase = abpr'>params.minbase;

[windx,windy,wind,windn,wind3,windn3] = make_kernel(params);
windn=windn(1:41,:);
windn3=windn3(:,1:41,:);

% 
% smallx=1:params.rx/2:nx;
% smally=[params.ry*2 params.ry*2+[1:params.ry/2:params.dely params.dely+1]];
% [sy,sx] = meshgrid(smally,smallx);

i0=1+floor(length(windx)/2);
j0=1+floor(length(windy)/2);

cpx = load_slc_chunk(params,slcnames,i-floor(length(windx)/2),i+floor(length(windx)/2),j-floor(length(windy)/2),j+floor(length(windy)/2),'cpx');
disp('done loading')

[gamma,ints,cors,hp]=make_cor(cpx,intid,wind,windn,wind3,windn3,params);
disp('done making cor')

%invert cor on coarse grid

data  = squeeze(cors(:,i0,j0))';
d     = log(data);
stds  = sqrt(-2*d);

[c0,cp,mk0,synth0] = invert_m_mat(d',nd,Gi0,intid);
mk0(mk0==0)=NaN;
mk0 = mk0+mymax(-mk0,3);
mk0(mk0<0)=0;
mk0(isnan(mk0))      = 0;

%d2                   = d'-c0-Gi0*cp;
%d2(isnan(synth0))    = NaN;

good    = isfinite(synth0);
start   = [c0;cp;mk0];
LB      = [-inf(nd,1);zeros(nd,1)];
UB      = [zeros(nd,1);inf(nd,1)];
[mod,~,resnew]     = lsqnonlin('corfit_all',start,LB,UB,params.OPTIONS,d(good),stds(good)',Gi0(good,:),Gr0(good,:),1);
[res1,synth] = corfit_all(mod,d,stds',Gi0,Gr0,1);

c0=mod(1);
cp=mod(2:nd);
mk1=mod(nd+1:end);

mk1         = mk1+ mymax(-mk1,5);
mk1(mk1<0)  = -mk1(mk1<0);


res   = d'-synth;
%good  = good & res<0.2 & synth>log(mincor) & d'>log(mincor);

mcor_orig = sqrt(mean(d(good).^2,'omitnan'));
mcor_res  = sqrt(mean(res(good).^2,'omitnan'));

meank1       = atan(mk1);
cdiag = d(diags)-c0;

%pull out all intervals with diag coherence of mincor or lower
permbad=cdiag<log(mincor);

%which ints are bad?
tmp = Gi0*permbad';
good = or(tmp==0 & isfinite(synth0),dt==1); %keep the ones on the diagonal for reference

%             synth = c0+Gi0*cp0;
%             res   = d'-synth;
%             good  = good & res<0.2 & synth>log(mincor) & d'>log(mincor);
%             [res,synth,dc0] = corfit(mk1,res(good),stds(good)',Gr0(good,:),1);
%             synth           = synth+c00+dc0+Gi0(good,:)*cp0;

%now do hp fit
aph    = squeeze(hp(:,i0,j0));

ph     = transpose(invert_phsmat(aph(good),nd,intid(good,:)));
ph(isnan(ph))=1;
aph2   = ph(intid(:,2)).*conj(ph(intid(:,1)));
res    = aph.*conj(aph2);
res    = res./abs(res);
tphs_res=abs(mean(res(good)));

allmk=diff(mk(intid),[],2);

%allhp=ph;
%alllp=[0;squeeze(gamma(shortid,i,k))];

ts     = linspace(-8,8,100);
gid    = and(mk1>std(mk1),isfinite(angle(ph)));

%gi2    = gi; %not storms
gi2=1:nd;

ph=ph./abs(ph);
ph(~isfinite(ph))=1;
if(sum(gid)>1)
    synths = mk1*ts;
    dat    = ph.*exp(1j*meank1);
    res    = dat.*conj(exp(1j*synths));
    fit    = sqrt(mean(angle(res(gid,:)).^2,1));
    besti  = find(fit==min(fit));
    besti  = round(mean(besti));
    slopes0k1=ts(besti);
    
elseif(sum(gid)==1)
    slopes0k1=angle(ph(gid)*exp(1j*meank1(gid)))/mk1(gid);
else
    slopes0k1=0;
end

[m1a,r1]  = lsqnonlin('phsfit_noint',slopes0k1,[],[],params.OPTIONS,mk1(gi2),ph(gi2),meank1(gi2),atan(mk1(gi2)));
[m1b,r2]  = lsqnonlin('phsfit_noint',0*slopes0k1,[],[],params.OPTIONS,mk1(gi2),ph(gi2),meank1(gi2),atan(mk1(gi2)));
if(r1>r2)
    slopesk1=m1b;
else
    slopesk1=m1a;
end
resk1 = phsfit_noint(slopesk1,mk1,ph,meank1,ones(size(mk1)));


tphs_orig     = abs(mean(ph));
tphs_resk      = abs(mean(exp(1j*resk1)));

