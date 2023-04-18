function fids=open_files(filenames,rw)

nd=length(filenames.rel);
ks=length(filenames.slope);

fids.bpslope   = fopen(filenames.bpslope,rw);
fids.c0        = fopen(filenames.c0,rw);
fids.mcor_orig = fopen(filenames.mcor_orig,rw);
fids.mcor_res  = fopen(filenames.mcor_res,rw);

fids.tphs_orig = fopen(filenames.tphs_orig,rw);

for i=1:nd-1
    fids.perm(i)=fopen(filenames.perm{i},rw);
end
for i=1:nd
    fids.rel(i)=fopen(filenames.rel{i},rw);
    fids.hp(i)=fopen(filenames.hp{i},rw);
    for j=1:ks
        fids.mk(i,j)=fopen(filenames.mk{i,j},rw);
    end
end
for j=1:ks
    fids.tphs_res(j)  = fopen(filenames.tphs_res{j},rw);
    fids.slope(j)     = fopen(filenames.slope{j},rw);
end
