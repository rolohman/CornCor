function [dn,nd,intid,dt,ni,id1,id2,diags]=make_pairs(dates,maxdt,dt1)

dn           = [dates.dn]';
nd           = length(dates);
intid        = nchoosek(1:nd,2);
dt           = diff(intid,[],2);


use    = or(dt<=maxdt,intid(:,1)<=dt1);
intid  = intid(use,:);
dt     = dt(use);

ni           = length(intid);
id1          = intid(:,1);
id2          = intid(:,2);
diags        = find(dt==1);