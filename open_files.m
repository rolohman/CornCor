function fids=open_files(params,dates,nd)

adir           = [params.outdir 'models/'];
fids.bpslope   = fopen([adir 'corbpslope.r4'],'w');
fids.c0        = fopen([adir 'c0.r4'],'w');
fids.mcor_orig = fopen([adir 'mcor_orig.r4'],'w');
fids.mcor_res  = fopen([adir 'mcor_res.r4'],'w');

fids.tphs_orig = fopen([adir 'tphs_orig.r4'],'w');

for i=1:nd-1
    fids.perm(i)=fopen([adir 'perm/' dates(i).name '_' dates(i+1).name '.r4'],'w');
end
for i=1:nd
    fids.rel(i)=fopen([adir 'rel/' dates(i).name '.r4'],'w');
    fids.hp(i)=fopen([adir 'hp/' dates(i).name '.hp.r4'],'w');
    for j=1:length(params.k)
        fids.mk(i,j)=fopen([adir 'mk/' dates(i).name '.mk' num2str(params.k(j))],'w');
    end
end
for j=1:length(params.k)
    fids.tphs_res(j)  = fopen([adir 'tphs_resk' num2str(params.k(j)) '.r4'],'w');
    fids.slope(j)     = fopen([adir 'slopesk' num2str(params.k(j)) '.r4'],'w');
end
