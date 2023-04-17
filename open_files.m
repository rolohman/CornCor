function fids=open_files(params,dates,nd,rw)

adir           = [params.outdir 'models/'];
fids.bpslope   = fopen([adir 'corbpslope.r4'],rw);
fids.c0        = fopen([adir 'c0.r4'],rw);
fids.mcor_orig = fopen([adir 'mcor_orig.r4'],rw);
fids.mcor_res  = fopen([adir 'mcor_res.r4'],rw);

fids.tphs_orig = fopen([adir 'tphs_orig.r4'],rw);

for i=1:nd-1
    fids.perm(i)=fopen([adir 'perm/' dates(i).name '_' dates(i+1).name '.r4'],rw);
end
for i=1:nd
    fids.rel(i)=fopen([adir 'rel/' dates(i).name '.r4'],rw);
    fids.hp(i)=fopen([adir 'hp/' dates(i).name '.hp.r4'],rw);
    for j=1:length(params.k)
        fids.mk(i,j)=fopen([adir 'mk/' dates(i).name '.mk' num2str(params.k(j))],rw);
    end
end
for j=1:length(params.k)
    fids.tphs_res(j)  = fopen([adir 'tphs_resk' num2str(params.k(j)) '.r4'],rw);
    fids.slope(j)     = fopen([adir 'slopesk' num2str(params.k(j)) '.r4'],rw);
end
