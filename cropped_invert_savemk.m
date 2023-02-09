init_dir(params)

[dates,slcnames]   = list_slcs(params);

%make all pairs
dn           = [dates.dn]';
nd           = length(dates);
intid        = nchoosek(1:nd,2);
dt           = diff(intid,[],2);
ni           = length(intid);
id1          = intid(:,1);
id2          = intid(:,2);
diags        = find(dt==1);

[Gi0,Gr0]    = make_G(ni,nd,id1,id2);

load baselines.txt
bpr     = baselines(id1)-baselines(id2);
abpr    = abs(bpr);
 
fids = open_files;

smallx=rx/2:rx:nx;
nsx=length(smallx);

[windx,windy,wind,windn,wind3,windn3] = make_kernel(params);

for j=1:ry:ny
    j
   
    cpx = load_slc_chunk(params,slcnames,1,nx,j-ry,j+ry+dely,'cpx');
    disp('done loading')

    [gamma,ints,cors,hp]=make_cor(cpx,ni,nx,dely,wind,windn3,windn);
    disp('done making cor')
    
    logcor   = log(cors);
    basewgts = exp(bpw*logcor);
  
    k=1;
    for i=1:length(smallx)            
            data  = squeeze(cors(:,smallx(i),k))';
            d     = log(data);
            g     = and(isfinite(d),bigbase);
            
            wgts  = diag(basewgts(g,smallx(i),k)); %wgts that allow robust max val.
            Gbase = [abpr(g).^2 ones(size(abpr(g)))];
            Ggbase= inv(Gbase'*wgts*Gbase)*Gbase'*wgts;
            
            basemod = Ggbase*d(g)';
            
            bps   = min(basemod(1),0);
          
            
            d       = d-bps*abpr'.^2;
            
            [c00,cp0,cr0]     = invert_cormat(d',nd,Gi0);
            
            cr0=cr0-mymax(cr0,50);
            allcr(i,:,k)=exp(cr0);
            cr0(cr0>0)=0;
            
            mk1_init=sqrt(1./exp(cr0).^2-1);
            mk2_init=2*sqrt(1./exp(cr0)-1);
            
            
            start   = [c00;cp0';mk1_init(2:end)'];
            LB      = [-inf(nd,1);zeros(nd-1,1)];
            UB      = [zeros(nd,1);inf(nd-1,1)];
            modk1 = lsqnonlin('corfit',start,LB,UB,OPTIONS,logcor(:,smallx(i),k),Gr0,Gi0,nd,1);
            start   = [c00;cp0';mk2_init(2:end)'];
            modk2 = lsqnonlin('corfit',start,LB,UB,OPTIONS,logcor(:,smallx(i),k),Gr0,Gi0,nd,2);
            allc00(i,k)=c00;
            allc0(i,k)=modk1(1);
            allcp(i,:,k)=modk1(2:nd);
            mk1 = [0;modk1(nd+1:end)];
            mk2 = [0;modk2(nd+1:end)];

            
            %mk1=mk1+mymax(-mk1,10);
            %mk2=mk2+mymax(-mk2,20);
            allmk1(i,:,k)=mk1;
            allmk2(i,:,k)=mk2;
         
    end  
    %write output files for this chunk

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
  
end


fclose('all');


