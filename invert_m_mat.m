function [c0,cp,mk,synth,res] = invert_m_mat(cor,nd,Gi0,intid)
%expects log cor
d  = cor;
c0 = mymax(d(:),100);
d  = d-c0;


cpmin = zeros(1,nd-1);
w     = exp(100*d);
dw    = (d.*w);

for k=1:nd-1
    id       = and(isfinite(dw),Gi0(:,k)==1);
    cpmin(k) = sum(dw(id))/sum(w(id));
end
cpmin       = min(0,cpmin);
cp          = est_ct(d,Gi0,cpmin,intid);

dp = Gi0*cp';
d = d-dp;

id1=sub2ind([nd,nd],intid(:,1),intid(:,2));
id2=sub2ind([nd,nd],intid(:,2),intid(:,1));

crs=zeros(nd,nd);
crs(id1)=d;
crs(id2)=d;

crs(crs==0)      = NaN;
crs(crs>0)       = NaN;
crs = sqrt(1./exp(crs).^2-1);
mk0      = median(crs,1,'omitnan');

crsnew  = crs-abs(repmat(mk0,nd,1)-repmat(mk0',1,nd))+repmat(mk0,nd,1);
mk2      = median(crsnew,1,'omitnan')';

mk3      = mk2 + mymax(-mk2,50);
mk3(mk3<0) = 0;
%clear crs crsnew
mk=mk3;

crs2=abs(mk(intid(:,2))-mk(intid(:,1)));
synth=log(1./sqrt((1+crs2.^2)));

d=cor-synth;
c0=median(d);
synth=synth+c0+dp;
res=cor-synth;