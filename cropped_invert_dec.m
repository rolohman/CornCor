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
gi           = sum(ismember(intid,[42 43]),2)==0; %for VV, don't use 2nd storm
gi           = ~ismember(1:nd,[42 43]); %for VV, don't use 2nd storm

OPTIONS       = optimoptions('lsqnonlin','MaxIterations',50,'Display','none');
adir          = [slcdir 'models/'];

slopenames = {'k1','k2'};

slopenum   = [1 2];
if(0)
    for i=1:length(slopenames)
        fidslope(i)    = fopen([adir 'slopes' slopenames{i} '.r4'],'w');
        fidslope0(i)    = fopen([adir 'slopes0' slopenames{i} '.r4'],'w');
        fidcor(i)      = fopen([adir 'cor' slopenames{i} '.r4'],'w');
        fidcorb(i)      = fopen([adir 'corb' slopenames{i} '.r4'],'w');

    end
end
if(1)
    fidr               = fopen([adir 'R.r4'],'w');

    %fido               = fopen([adir 'ocor.r4'],'w');
    fidrb               = fopen([adir 'Rb.r4'],'w');
    %fidob               = fopen([adir 'ocorb.r4'],'w');
end
dely    = 25;
smallx  = rx/2:rx:nx;
smally  = 1:ry:ny;
[sy,sx] = meshgrid(smally,smallx);

allmk1 = zeros(length(smallx),length(smally),nd);
for i=1:nd
    fid             = fopen([adir 'mk/' dates(i).name '.mk1'],'r');
    tmp             = fread(fid,[length(smallx),length(smally)],'real*4');
    tmp(tmp>10)     = 10;
    tmp(isnan(tmp)) = 0;
    allmk1(:,:,i)   = tmp;
    fclose(fid);
end
allmk2 = zeros(length(smallx),length(smally),nd);
for i=1:nd
    fid             = fopen([adir 'mk/' dates(i).name '.mk2'],'r');
    tmp             = fread(fid,[length(smallx),length(smally)],'real*4');
    tmp(tmp>10)     = 10;
    tmp(isnan(tmp)) = 0;
    allmk2(:,:,i)   = tmp;
    fclose(fid);
end

windn        = conv2(ones(nx,dely+ry*2+1),wind,'same');
windn3(1,:,:)=windn;
wind3(1,:,:)=wind;
for j=1:dely:ny
    j
    x1    = 1;
    x2    = nx;
    
    y1    = j-ry;
    y2    = j+ry+dely;
    newnx = x2-x1+1;
    [Y,X] = meshgrid(j-1+(1:dely),1:nx);
    
    cpx=load_chunk(slcnames,nx,ny,x1,x2,y1,y2,'cpx');
    disp('done loading')
    
    gamma  = nan(ni,newnx,dely);
    ints   = nan(ni,newnx,dely);
    amps   = cpx.*conj(cpx);
    ampsum = sqrt(convn(amps,wind3,'same')./windn3);
    for i=1:ni
        slc1 = squeeze(cpx(intid(i,1),:,:));
        slc2 = squeeze(cpx(intid(i,2),:,:));
        
        c    = slc2.*conj(slc1);
        
        asum = squeeze(ampsum(intid(i,1),:,:));
        bsum = squeeze(ampsum(intid(i,2),:,:));
        csum = conv2(c,wind,'same')./windn;
        cpx3 = csum./asum./bsum;
        cpx3(isnan(cpx3))=0;
        
        gamma(i,:,:) = cpx3(:,ry+[1:dely]);
        ints(i,:,:)  = c(:,ry+[1:dely]);
    end
    
    disp('done making gamma')
    
    hp       = ints.*conj(gamma);
    hp       = hp./abs(hp);
    hp(isnan(hp))=1;
  
    ocor        = nan(newnx,dely);
    cork1       = nan(newnx,dely);
    cork2       = nan(newnx,dely);
    slopesk1    = nan(newnx,dely);
    slopesk2    = nan(newnx,dely);

    amk1        = zeros(newnx,dely,nd);
    amk2        = zeros(newnx,dely,nd);
    GINT        = griddedInterpolant(sx,sy,squeeze(allmk1(:,:,1)));
    for i=2:nd
        GINT.Values = squeeze(allmk1(:,:,i));
        amk1(:,:,i) = GINT(X,Y);
        GINT.Values = squeeze(allmk2(:,:,i));
        amk2(:,:,i) = GINT(X,Y);
    end
    
    for i=1:newnx
        i
        
        for k=1:dely
            mk1         = squeeze(amk1(i,k,:));
            mk2         = squeeze(amk2(i,k,:));
            
            meank1      = atan(mk1);
            meank2      = atan(mk2./(1-0.25*mk2.^2));
            neg         = 1<0.25*mk2.^2;
            meank2(neg) = meank2(neg)+pi;
            
            %now do hp fit
            aph    = squeeze(hp(:,i,k));
            ph     = transpose(invert_phsmat(aph,nd));
            pa1    = angle(ph.*exp(1j*meank1)); %for calculation of R
            pa2    = angle(ph.*exp(1j*meank2)); %for calculation of R
            
            ts     = linspace(-8,8,100);
            gid    = and(mk1>0.4,isfinite(angle(ph)));
            gi2    = gi; %not storms
         
            if(sum(gid)>1)
                synths = mk1*ts;
                dat    = ph.*exp(1j*meank1);
                res    = dat.*conj(exp(1j*synths));
                fit    = sqrt(mean(angle(res(gid,:)).^2,1));
                besti  = find(fit==min(fit));
                besti  = round(mean(besti));
                slopes0k1(i,k)=ts(besti);
                
                synths = mk2*ts;
                dat    = ph.*exp(1j*meank2);
                res    = dat.*conj(exp(1j*synths));
                fit    = sqrt(mean(angle(res(gid,:)).^2,1));
                besti  = find(fit==min(fit));
                besti  = round(mean(besti));
                slopes0k2(i,k)=ts(besti);
            elseif(sum(gid)==1)
                slopes0k1(i,k)=angle(ph(gid)*exp(1j*meank1(gid)))/mk1(gid);
                slopes0k2(i,k)=angle(ph(gid)*exp(1j*meank2(gid)))/mk2(gid);
            else
                slopes0k1(i,k)=0;
                slopes0k2(i,k)=0;
            end
  
          [m1a,r1]  = lsqnonlin('phsfit_noint',slopes0k1(i,k),[],[],OPTIONS,mk1(gi2),ph(gi2),meank1(gi2),mk1(gi2));
          [m1b,r2]  = lsqnonlin('phsfit_noint',0*slopes0k1(i,k),[],[],OPTIONS,mk1(gi2),ph(gi2),meank1(gi2),mk1(gi2));
          if(r1>r2)
                  slopesk1(i,k)=m1b;
          else
                  slopesk1(i,k)=m1a;
          end
          resk1 = phsfit_noint(slopesk1(i,k),mk1,ph,meank1,ones(size(mk1)));

     %     [m2a,r1]  = lsqnonlin('phsfit_noint',slopes0k2(i,k),[],[],OPTIONS,mk2(gi2),ph(gi2),meank2(gi2),mk2(gi2));
     %     [m2b,r2]  = lsqnonlin('phsfit_noint',0*slopes0k2(i,k),[],[],OPTIONS,mk2(gi2),ph(gi2),meank2(gi2),mk2(gi2));
     %     if(r1>r2)
     %             slopesk2(i,k)=m2b;
     %     else
     %     
%	  	  slopesk2(i,k)=m2a;
%          end 
%	  resk2 = phsfit_noint(slopesk2(i,k),mk2,ph,meank2,ones(size(mk2)));

 %           ocor(i,k)       = abs(mean(ph));
 %           cork1(i,k)      = abs(mean(exp(1j*resk1)));
 %           cork2(i,k)      = abs(mean(exp(1j*resk2)));
  %          ocorb(i,k)      = abs(mean(ph(gid)));
  %          corbk1(i,k)     = abs(mean(exp(1j*resk1(gid))));
  %          corbk2(i,k)     = abs(mean(exp(1j*resk2(gid))));
n=sum(gid);
if(n>2)
	a1=mean(angle(ph));
	a2=mean(angle(resk1));
  R(i,k)=1-(sum((angle(ph)-a1).^2)/sum((angle(resk1)-a2).^2));
Rb(i,k)=(1-(1-R(i,k))*(n-1)/(n-2));
else
	R(i,k)=NaN;
	Rb(i,k)=NaN;
end
	end
    end

    
    %write output files for this chunk
    
    %fwrite(fido,ocor,'real*4');
    fwrite(fidr,R,'real*4');
    %fwrite(fidcor(1),cork1,'real*4');
    %fwrite(fidcor(2),cork2,'real*4');    
    %fwrite(fidob,ocorb,'real*4');
    fwrite(fidrb,Rb,'real*4');
    %5fwrite(fidcorb(1),corbk1,'real*4');
    %fwrite(fidcorb(2),corbk2,'real*4');
    %fwrite(fidslope(1),slopesk1,'real*4');
    %fwrite(fidslope(2),slopesk2,'real*4');
    %fwrite(fidslope0(1),slopes0k1,'real*4');
    %fwrite(fidslope0(2),slopes0k2,'real*4');
end


fclose('all');


