function filenames=make_filenames(params,dates,nd)

adir           = [params.outdir 'models/'];
filenames.bpslope   = [adir 'corbpslope.r4'];
filenames.c0        = [adir 'c0.r4'];
filenames.mcor_orig = [adir 'mcor_orig.r4'];
filenames.mcor_res  = [adir 'mcor_res.r4'];

filenames.tphs_orig = [adir 'tphs_orig.r4'];
filenames.tphs_res  = [adir 'tphs_res.r4'];

for i=1:nd-1
    filenames.perm{i}=[adir 'perm/' dates(i).name '_' dates(i+1).name '.r4'];
end
for i=1:nd
    filenames.hp{i}=[adir 'hp/' dates(i).name '.hp.r4'];
    filenames.lp{i}=[adir 'lp/' dates(i).name '.lp.r4'];
    
    for j=1:length(params.k)
        filenames.mk{i,j}=[adir 'mk/' dates(i).name '.mk' num2str(params.k(j))];
    end
end
for j=1:length(params.k)
    filenames.tphs_resk{j}  = [adir 'tphs_resk' num2str(params.k(j)) '.r4'];
    filenames.slope{j}     = [adir 'slopesk' num2str(params.k(j)) '.r4'];
end
