function init_dir(params)
outdir = params.outdir;

if(~exist(outdir,'dir'))
    disp([outdir ' doesn''t exist yet - creating.'])
    mkdir(outdir)
    adir     = [outdir 'models/'];
    if(~exist(adir,'dir'))
         mkdir(adir)
    end
    if(~exist([adir 'mk'],'dir'))
        mkdir([adir 'mk']);
    end
    if(~exist([adir 'rel'],'dir'))
        mkdir([adir 'rel']);
    end
    if(~exist([adir 'perm'],'dir'))
        mkdir([adir 'perm']);
    end
    save([outdir '/paramfile.mat'],'params');
else
    disp([outdir ' already exists, is that okay?'])    
end

