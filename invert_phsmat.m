function  phsmod = invert_phsmat(phs,nd,intid)

id1=sub2ind([nd,nd],intid(:,2),intid(:,1));
id2=sub2ind([nd,nd],intid(:,1),intid(:,2));

crs=zeros(nd,nd);
crs(id1)=conj(phs);
crs(id2)=phs;

crs(crs==0)      = NaN;

cr               = mean(crs,2,'omitnan');
crsnew           = crs.*conj(repmat(cr,1,nd));
phsmod           = mean(crsnew,1,'omitnan');
%phsmod           = phsmod.*conj(phsmod(1));
