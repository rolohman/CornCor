function [c0,cp,cr,synth,res] = invert_cormat(cor,nd,Gi0,intid)
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
cpmin          = min(0,cpmin);
cp1          = est_ct(d,Gi0,cpmin,intid);

dp = Gi0*cp1';
d = d-dp;


id1=sub2ind([nd,nd],intid(:,1),intid(:,2));
id2=sub2ind([nd,nd],intid(:,2),intid(:,1));

crs=zeros(nd,nd);
crs(id1)=d;
crs(id2)=d;

crs(crs==0)      = NaN;

cr      = median(crs,1,'omitnan');
pcr     = cr';

crsnew  = -abs(crs-(-abs(repmat(cr,nd,1)-repmat(pcr,1,nd))-repmat(pcr,1,nd)'));
cr      = median(crsnew,1,'omitnan');
clear crs crsnew

[crm,cp]       = flatten_back(cr,25,cp1,cpmin);
cr = cr-crm;


[crm,cp]       = flatten_front(cr,25,cp,cpmin);
cr=cr-crm;
cr           = cr-mymax(cr,50);

synth=c0+Gi0*cp'-abs(diff(cr(intid),[],2));
res=d-synth;
