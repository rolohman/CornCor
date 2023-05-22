function fids=open_files(filenames,rw)

nd=length(filenames.hp);
ks=length(filenames.slope);

fids.c0        = fopen(filenames.c0,rw);
fids.mcor_orig = fopen(filenames.mcor_orig,rw);
fids.mcor_res  = fopen(filenames.mcor_res,rw);

fids.tphs_orig = fopen(filenames.tphs_orig,rw);
fids.tphs_res  = fopen(filenames.tphs_res,rw);

for i=1:nd-1
    fids.perm(i)=fopen(filenames.perm{i},rw);
end
for i=1:nd
   fids.hp(i)=fopen(filenames.hp{i},rw);
   fids.lp(i)=fopen(filenames.lp{i},rw);
    for j=1:ks
        fids.mk(i,j)=fopen(filenames.mk{i,j},rw);
    end
end
for j=1:ks
    fids.tphs_resk(j)  = fopen(filenames.tphs_resk{j},rw);
    fids.slope(j)     = fopen(filenames.slope{j},rw);
end
