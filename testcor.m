load teststuff

data  = cors;
d     = log(data)';
define_params
d(d>0)=0;
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
figure
subplot(2,2,1)
triplot(d(good),dn,intid(good,:)),caxis(c)
subplot(2,2,2)
triplot(synth0,dn,intid),caxis(c)
subplot(2,2,3)
triplot(synth,dn,intid(good,:)),caxis(c)

figure
subplot(1,2,1)
plot(mk0)
hold on
plot(modk1)

subplot(1,2,2)
plot(cp0)
