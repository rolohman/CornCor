function triplot(data,dn,intid)
nd=length(dn);

id2=sub2ind([nd,nd],intid(:,2),intid(:,1));

jnk=nan(nd);

jnk(id2)=data;
pcolor(dn,dn,jnk'),shading flat,set(gca,'ydir','reverse');
