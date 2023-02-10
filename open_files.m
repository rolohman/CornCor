function fids=open_files(params)

adir           = [params.outdir 'models'];
fids.bpslope   = fopen([adir 'corbpslope.r4'],'w');
fids.c0        = fopen([adir 'c0.r4'],'w');
fids.tcor_orig = fopen([adir 'tcor_orig.r4'],'w');
fids.tphs_orig = fopen([adir 'tphs_orig.r4'],'w');

for i=1:nd-1
    fids.perm(i)=fopen([adir 'perm/' dates(i).name '_' dates(i+1).name '.r4'],'w');
end
for i=1:nd
    fids.rel(i)=fopen([adir 'rel/' dates(i).name '.r4'],'w');
    for j=1:length(params.k)
        fids.mk(i,j)=fopen([adir 'mk/' dates(i).name '.mk' num2str(params.k(j))],'w');
    end
end
for j=1:length(params.k)
    fids.tcor_res(j)  = fopen([adir 'tcor_resk' num2str(params.k(j)) '.r4'],'w');
    fids.tphs_res(j)  = fopen([adir 'tphs_resk' num2str(params.k(j)) '.r4'],'w');
    fids.slope(i)     = fopen([adir 'slopesk' num2str(params.k(j)) '.r4'],'w');
end