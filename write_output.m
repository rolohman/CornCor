function write_output(fids,allbps,mcor_orig,mcor_res,allc0,allcr,allmk1,allcp,slopesk1,tphs_orig,tphs_res,allhp)
nd=63;
    fwrite(fids.c0,exp(allc0),'real*4');
    fwrite(fids.bpslope,allbps,'real*4');
    
    fwrite(fids.mcor_orig,mcor_orig,'real*4');
    fwrite(fids.mcor_res, mcor_res,'real*4');
    
    fwrite(fids.tphs_orig,tphs_orig,'real*4');
    
    for i=1:nd
        fwrite(fids.rel(i),exp(allcr(:,i)),'real*4');
        fwrite(fids.mk(i,1),allmk1(:,i),'real*4');
        fwrite(fids.hp(i),angle(allhp(:,i)),'real*4');
    end
    for i=1:nd-1
        fwrite(fids.perm(i),exp(allcp(:,i)),'real*4');
    end
    j=1;
    %for j=1:length(params.k)
        fwrite(fids.tphs_res(j),tphs_res,'real*4');
        fwrite(fids.slope(j),slopesk1,'real*4');
    %end
