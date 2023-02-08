home=pwd;
params

slcdir='cropped_38370_11488_4800_800/SLC_VV/';

files=dir([slcdir '2*']);
dates=[];
for i=1:length(files)
    dates(i).name=files(i).name(1:8);
    dates(i).dn=datenum(dates(i).name,'yyyymmdd');
    slcnames{i}=[slcdir dates(i).name '/' dates(i).name '.slc.full'];
 end
dn = [dates.dn]';
nd     = length(dates);
%make all pairs

nx=4800;
ny=800;

rx           = 10;
ry           = 10;

windx        = exp(-(-rx*2:rx*2).^2/2/(rx/2)^2);
windy        = exp(-(-ry*2:ry*2).^2/2/(ry/2)^2);

windx        = windx/sum(windx);
windy        = windy/sum(windy);
wind         = windy'*windx;
     

%make all pairs
intid        = nchoosek(1:nd,2);
dt           = diff(intid,[],2);
ni           = length(intid);
id1          = intid(:,1);
id2          = intid(:,2);
diags        = find(dt==1);

load baselines.txt
bpr=baselines(id1)-baselines(id2);
abpr=abs(bpr);
bigbase = abpr'>50;

%nbin=20;
%[basebins,basebinedge]=discretize(abpr,nbin);
%basebinmid=basebinedge(1:end-1)+diff(basebinedge)/2;
bpw=200; %weights for baseline estimation

OPTIONS = optimoptions('lsqnonlin','MaxIterations',50,'Display','none');
    


adir     = [slcdir 'models/'];
if(~exist(adir,'dir'))
    disp(['creating directory ' adir]);
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

for i=1:nd
    fidr(i)=fopen([adir 'rel/' dates(i).name '.cor'],'w');
    fidk1(i)=fopen([adir 'mk/' dates(i).name '.mk1'],'w');
    fidk2(i)=fopen([adir 'mk/' dates(i).name '.mk2'],'w');
end
for i=1:nd-1
    fidp(i)=fopen([adir 'perm/' dates(i).name '_' dates(i+1).name '.cor'],'w');   
end
fid0=fopen([adir 'c0.cor'],'w');
fid00=fopen([adir 'c00.cor'],'w');

cidl     = tril(ones(nd),-1)==1;
errslope = 0.3;
mncor    = errslope/(1+errslope);
mncorl   = log(mncor);

Gi0=zeros(ni,nd-1); %intervals, perm cor loss
for i=1:ni
    Gi0(i,id1(i):id2(i)-1)=1;
end
Gr0=zeros(ni,nd); %rel cor on dates
for i=1:ni
    Gr0(i,id1(i))=1;
    Gr0(i,id2(i))=-1;
end
Gr=Gr0(:,2:end);

dely=2;
smallx=rx/2:rx:nx;

nsx=length(smallx);

windn        = conv2(ones(nx,dely+ry*2+1),wind,'same');
windn3(1,:,:)=windn;
wind3(1,:,:)=wind;
for j=1:ry:ny
    j
    x1=1;
    x2=nx;
    
    y1=j-ry;
    y2=j+ry+dely;
    newnx=x2-x1+1;
    cpx=load_chunk(slcnames,nx,ny,x1,x2,y1,y2,'cpx');
    disp('done loading')
  
    %cpx=cpx./abs(cpx);
    %cpx(isnan(cpx))=0;
    gamma=nan(ni,newnx,dely);
    ints   = nan(ni,newnx,dely);
    amps=cpx.*conj(cpx);
    ampsum=sqrt(convn(amps,wind3,'same')./windn3);
    cpx=cpx./abs(cpx);
    cpx(isnan(cpx))=1;
    for i=1:ni
        slc1=squeeze(cpx(intid(i,1),:,:));
        slc2=squeeze(cpx(intid(i,2),:,:));
       

        c   = slc2.*conj(slc1);
        asum = squeeze(ampsum(intid(i,1),:,:));
        bsum = squeeze(ampsum(intid(i,2),:,:));
        csum = conv2(c,wind,'same')./windn;
        cpx3 = csum./asum./bsum;
        %cpx3  = csum;
        cpx3(isnan(cpx3))=0;
        
        gamma(i,:,:)=cpx3(:,ry+[1:dely]);
        ints(i,:,:)=c(:,ry+[1:dely]);
    end

    disp('done making cor')
    
    cors   = abs(gamma);
    good   = cors>0;
    cors(cors>1)=1;
    cors(~good)=NaN;
    ngood    = squeeze(sum(good,1));
    goodid   = ngood>=10;
    hp       = ints.*conj(gamma);
    hp       = hp./abs(hp);
    hp(isnan(hp))=1;
    logcor   = log(cors);
    basewgts = exp(bpw*logcor);
    
    allc0=nan(nsx,dely);
    allc00=nan(nsx,dely);
    allcr=nan(nsx,nd,dely);
    allcp=nan(nsx,nd-1,dely);
    allmk1=nan(nsx,nd,dely);
    allmk2=nan(nsx,nd,dely);
 
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


