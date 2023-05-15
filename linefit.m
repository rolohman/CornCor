
clear windn3 
define_params
params.slcdir='cropped_41000_11942_375_200/SLC_VV/';
params.outdir=[params.slcdir 'trial2/'];

nx=params.nx;
ny=params.ny;
rx=params.rx;
ry=params.ry;

init_dir(params)

[dates,slcnames]                  = list_slcs(params);
nd=length(dates);
ks=length(params.k);

filenames                         = make_filenames(params,dates,length(dates));
fids                              = open_files(filenames,'r');


bpslope   = fread(fids.bpslope,[nx,ny],'real*4');
c0        = fread(fids.c0,[nx,ny],'real*4');

tphs_orig = fread(fids.tphs_orig,[nx,ny],'real*4');
tphs_res  = fread(fids.tphs_res,[nx,ny],'real*4');
for i=1:nd
    hp(i,:,:)=fread(fids.hp(i),[nx,ny],'real*4');
    for j=1:ks
        mk(i,j,:,:)=fread(fids.mk(i,j),[nx,ny],'real*4');
    end
end
for j=1:ks
    tphs_resk(j,:,:)  = fread(fids.tphs_resk(j),[nx,ny],'real*4');
    slope(j,:,:)     = fread(fids.slope(j),[nx,ny],'real*4');
end

%while we only have one k
tphs_resk=squeeze(tphs_resk);
slope=squeeze(slope);
mk=squeeze(mk);
ts     = linspace(-8,8,100);

for j=1:ny
    for i=1:nx
        ph=exp(1j*(squeeze(hp(:,i,j))));
        m=squeeze(mk(:,i,j));
        mf=atan(m);
        gid=m>std(m);
        synths = m*ts;
        ints(i,j)=angle(mean(ph(~gid)));
        %ph=ph.*conj(exp(1j*ints(i,j)));
        res    = ph.*conj(exp(1j*synths));
        res2=exp(1j*mf).*conj(exp(1j*synths));
        fit    = sqrt(mean(angle(res(gid,:)).^2,1));
        fit2=sqrt(mean(angle(res2(gid,:)).^2,1));
        besti  = find(fit==min(fit));
        besti  = round(mean(besti));
        
        newph(:,i,j)=angle(ph.*conj(exp(1j*synths(:,besti))));
        newslope(i,j)=ts(besti);
        besti=find(fit2==min(fit2));
        besti=round(mean(besti));
        newslope2(i,j)=ts(besti);
        newph2(:,i,j)=angle(exp(1j*mf).*conj(exp(1j*synths(:,besti))));
    end
end



m2=reshape(mk,nd,nx*ny);
p2=reshape(newph,nd,nx*ny);
p3=reshape(newph2,nd,nx*ny);

ms=linspace(0,4,40);
ps=linspace(-pi,pi,80);
p=ps(1:end-1)+diff(ps(1:2))/2;

good=find(tphs_resk(:)>0.6);
[h]=histcounts2(m2(good),p2(good),ms,ps);
h2=histcounts2(m2(good),p3(good),ms,ps);
clear bin
for i=1:size(h,1)
id=find(h(i,:)==max(h(i,:)));
bin(i)=round(mean(id));
id=find(m2(good)>ms(i) & m2(good)<=ms(i+1));
means(i)=mean(p2(good(id)));
meds(i)=median(p2(good(id)));
end

figure
subplot(1,2,1)
imagesc(ms(1:end-1),ps(1:end-1),log10(h)');
hold on
plot([0 4],[0 0],'k')
plot(ms(1:end-1),p(bin),'w.')
plot(ms(1:end-1),meds,'ro')
set(gca,'ydir','normal')
subplot(1,2,2)
imagesc(ms(1:end-1),ps(1:end-1),log10(h2)');
hold on
plot([0 4],[0 0],'k')
plot(ms(1:end-1),p(bin),'w.')
plot(ms(1:end-1),meds,'ro')
set(gca,'ydir','normal')

