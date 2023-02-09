function fids=open_files(params)

adir = [params.outdir 'models'];
for i=1:nd
    fids.fidr(i)=fopen([adir 'rel/' dates(i).name '.cor'],'w');
    fids.fidk1(i)=fopen([adir 'mk/' dates(i).name '.mk1'],'w');
    fids.fidk2(i)=fopen([adir 'mk/' dates(i).name '.mk2'],'w');
end
for i=1:nd-1
    fids.fidp(i)=fopen([adir 'perm/' dates(i).name '_' dates(i+1).name '.cor'],'w');   
end
fids.fid0=fopen([adir 'c0.cor'],'w');
fids.fid00=fopen([adir 'c00.cor'],'w');
