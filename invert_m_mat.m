function [c0,cp,mk,synth,res] = invert_m_mat(logcor,stds,nd,Gi0,intid)
%expects log cor
d  = logcor;
mincor = 0.45;

wgts=1./stds.^2;
wgts(~isfinite(wgts))=0;

d(d>0)  = 0;
%tmpd(tmpd<-1) = -1;
dt            = intid(:,2)-intid(:,1);
diags         = dt==1;

[c0,cp,~,synth1]=est_cp(d,Gi0,diags,mincor);
d        = d-synth1;
%d(d>0.2) = NaN; %correction clearly bad
good     = isfinite(d);

id1  = sub2ind([nd,nd],intid(good,1),intid(good,2));
id2  = sub2ind([nd,nd],intid(good,2),intid(good,1));

crs         = nan(nd,nd);
crs(id1)    = d(good);
crs(id2)    = d(good);

crs(crs>0)  = NaN;
crs         = sqrt(1./exp(crs).^2-1);
mk0         = median(crs,1,'omitnan');

crsnew      = crs-abs(repmat(mk0,nd,1)-repmat(mk0',1,nd))+repmat(mk0,nd,1);
mk2         = median(crsnew,1,'omitnan')';

mk3         = mk2 + mymax(-mk2,50);
mk3(mk3<0)  = 0;
%clear crs crsnew
mk          = mk3;

crs2   = abs(mk(intid(:,2))-mk(intid(:,1)));
synth2 = log(1./sqrt((1+crs2.^2)));

%d     = logcor-synth2;
%c0=median(d);
synth = synth1+synth2;
%synth(~good)=NaN;
res   = logcor-synth;
%bad   = res>0.2;
%res(bad)=NaN;
%synth(bad)=NaN;
