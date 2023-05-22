function [dn,nd,intid,dt,ni,id1,id2,diags]=make_pairs(dates,maxdt,dt1)

dn           = [dates.dn]';
nd           = length(dates);
intid        = nchoosek(1:nd,2);
dt           = diff(intid,[],2);
dday         = diff(dn(intid),[],2);

use     = or(dday<=maxdt,intid(:,1)<=dt1);
intid1  = intid(use,:);
intid2  = intid(~use,:);
if(size(intid2,1)>1)
    p       = randperm(size(intid2,1),min(size(intid2,1),2000));
    intid   = [intid1;intid2(p,:)];
else
    intid=intid1;
end
dt      = diff(intid,[],2);


ni           = length(intid);
id1          = intid(:,1);
id2          = intid(:,2);
diags        = find(dt==1);