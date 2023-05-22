function triplets = trip_nm(allintf,intid,n)

%allintf needs to be just one pol for now, so ni,ny,nx
[ni,ny,nx] = size(allintf);

triplets   = nan(ni,ny,nx);

for i=1:ni
    if(diff(intid(i,:))<=n)
        %do nothing
    else
        shortid1        = intid(i,1):n:intid(i,2)-n;
        shortid2        = shortid1+n;
        %shortid1(end+1) = shortid2(end);
        %shortid2(end+1) = intid(i,2);
        
        %shortids        = unique(find(ismember(intid,[shortid1' shortid2'],'rows')));
        shorts=[shortid1' shortid2'];
        
        shortids=ismember(intid,shorts,'rows');
      
        
        triplets(i,:,:) = prod(allintf(shortids,:,:),1).*conj(allintf(i,:,:));
    end
end