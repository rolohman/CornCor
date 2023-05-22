function write_output(fids,mcor_orig,mcor_res,allc0,allmk1,allcp,slopesk1,tphs_orig,tphs_res,tphs_resk,allhp,alllp)
nd=length(fids.hp);
ks=length(fids.slope);
fwrite(fids.c0,exp(allc0),'real*4');

fwrite(fids.mcor_orig,mcor_orig,'real*4');
fwrite(fids.mcor_res, mcor_res,'real*4');

fwrite(fids.tphs_orig,tphs_orig,'real*4');
fwrite(fids.tphs_res,tphs_res,'real*4');

for i=1:nd
    fwrite(fids.hp(i),angle(allhp(:,i)),'real*4');
    fwrite(fids.lp(i),angle(alllp(:,i)),'real*4');
    for j=1:ks
        fwrite(fids.mk(i,j),allmk1(:,i),'real*4');
    end
end
for i=1:nd-1
    fwrite(fids.perm(i),exp(allcp(:,i)),'real*4');
end
for j=1:ks
    fwrite(fids.tphs_resk(j),tphs_resk,'real*4');
    fwrite(fids.slope(j),slopesk1,'real*4');
end
