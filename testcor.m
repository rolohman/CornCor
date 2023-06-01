define_params
nx=params.nx;
ny=params.ny;

init_dir(params)
[dates,slcnames]                  = list_slcs(params);
params.dt1                        = length(dates);


[dn,nd,intid,dt,ni,id1,id2,diags] = define_pairs(dates,params.maxdt,params.dt1);
[Gi0,Gr0]                         = make_G(ni,nd,id1,id2);
n=params.rx*params.ry;
shortid=find(dt==1);

dt2=dn(intid(:,2))-dn(intid(:,1));

load baselines.txt
bpr     = baselines(id1)-baselines(id2);
abpr    = abs(bpr);
bigbase = abpr'>params.minbase;

%xpt=100;
%ypt=50;

cpx = load_slc_chunk(params,slcnames,xpt-5,xpt+5,ypt-5,ypt+5,'cpx');

amps   = cpx.*conj(cpx);
ampsum = mean(amps,[2 3]);

ampint=abs(diff(sqrt(ampsum(intid)),[],2));

ints=cpx(intid(:,2),:,:).*conj(cpx(intid(:,1),:,:));
gamma=mean(ints,[2 3])./sqrt(ampsum(intid(:,1)).*ampsum(intid(:,2)));

cors=abs(gamma);
d=log(cors);
cdiags=d(shortid);
stds=sqrt(-2*d);

c00=mymax(d,50);
mstart=[0;d(1:nd-1)-c00];
crm=flatten_back(mstart,50);
mstart=mstart-crm';
mstart(mstart>0)=0;
mstart=sqrt(1./(exp(mstart).^2)-1);
percstart=max(cors(dt2>1600));
as=sqrt(ampsum);
mstart=(as-min(as))/(max(as)-min(as));

res=d-c00;
for i=1:nd-1
    a=find(Gi0(:,i));
    cp0(i)=min(0,mymax(res(a),10));
end
cpstart=log(percstart+exp(cp0)*(1-percstart));

LB=[c00;cpstart';0];
UB=[max(d);zeros(nd-1,1);1];

start=[max(d)/2;cpstart'/5;0.6];
model=lsqnonlin('corfit_new2',start,LB,UB,[],d,stds,Gi0);
[res1,synth1,c0,cpp,scaleperc]=corfit_new2(model,cors,stds,Gi0);

LB=[c00;cpstart';zeros(nd,1);0];
UB=[max(d);zeros(nd-1,1);5*ones(nd,1);1];
start=[model(1:nd);mstart;percstart];
model=lsqnonlin('corfit_new',start,LB,UB,[],d,stds,Gi0,Gr0);
[res2,synth2]=corfit_new(model,cors,stds,Gi0,Gr0);

figure
subplot(2,3,1)
triplot(d,dn,intid),c=caxis;
subplot(2,3,2)
triplot(synth1,dn,intid),caxis(c)
subplot(2,3,3)
triplot(synth2,dn,intid),caxis(c)


subplot(4,3,7)
plot(model(2:nd))

subplot(4,3,10)
plot(model(nd+1:end-1))

subplot(2,3,5)
triplot(ampint,dn,intid);
subplot(2,3,6)
plot(sqrt(ampsum));


return



c2=cors-0.4;
c2(c2<0)=min(c2(c2>0));
d2=log(c2);

G=[ones(ni,1) Gi0];
G2=[ones(ni,1) log(0.6)+Gi0];


c2diags=d2(shortid);



return

lb=[-1;cdiags];
ub=zeros(nd,1);
lb2=[-1 ;c2diags];





mod1=lsqlin(G,d,[],[],[],[],lb,ub);
mod2=lsqlin(G2,d2,[],[],[],[],lb2,ub);

synth1=G*mod1;
synth2=G2*mod2;

figure
subplot(2,3,1)
triplot(d,dn,intid)
c=caxis;
subplot(2,3,2)
triplot(synth1,dn,intid)
caxis(c);
subplot(2,3,3)
triplot(synth2,dn,intid)
caxis(c)
subplot(2,3,4)
hist(d,100);
subplot(2,3,5)
hist(d-synth1,100);
subplot(2,3,6)
hist(d2-synth2,100);
return

stds=sqrt(-2*d);

[c00,cp0,mk0,synth0] = invert_m_mat(d',nd,Gi0,intid);
mk0(isnan(mk0))=0;

d2       = d'-c00-Gi0*cp0;
d2(isnan(synth0))=NaN;

good=isfinite(d2);
start   = [mk0];
LB      = [zeros(nd,1)];
UB      = [inf(nd,1)];
modk1   = lsqnonlin('corfit',start,LB,UB,params.OPTIONS,d2(good),stds(good)',Gr0(good,:),1);
[res,synth]      = corfit(modk1,d2(good),stds(good)',Gr0(good,:),1);
synth=synth+c00+Gi0(good,:)*cp0;

c=[-3 0];

