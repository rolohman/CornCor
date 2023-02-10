function write_output(fid,stuff)

    fwrite(fid0,exp(allc0(:,k)),'real*4');
    fwrite(fid00,exp(allc00(:,k)),'real*4');
    for i=1:nd
        fwrite(fidr(i),exp(squeeze(allcr(:,i,k))),'real*4');
        fwrite(fidk1(i),squeeze(allmk1(:,i,k)),'real*4');
        fwrite(fidk2(i),squeeze(allmk2(:,i,k)),'real*4');
        
    end
    for i=1:nd-1
        fwrite(fidp(i),exp(squeeze(allcp(:,i,k))),'real*4');
    end