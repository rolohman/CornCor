function data= load_chunk(filenames,nx,ny,x1,x2,y1,y2,type)

newdx=x2-x1+1;
newdy=y2-y1+1;
%type=1: cpx, type=2: rmg, type=3: r4

x0a=[];x0b=[];y0a=[];y0b=[];
if(y1<1)
    nystart=1-y1;
    y1=1;
    y0a=zeros(newdx,nystart);
end
if(y2>ny)
    nyend=y2-ny;
    y2=ny;
    y0b=zeros(newdx,nyend);
end
 if(x1<1)
    nxstart=1-x1;
    x1=1;
    x0a=zeros(nxstart,newdy);
end
if(x2>nx)
    nxend=x2-nx;
    x2=nx;
    x0b=zeros(nxend,newdy);
end   
    
n=length(filenames);
switch type
    case 'cpx' %cpx
        data=zeros(n,newdx,newdy);
        for i=1:n
            fid=fopen(filenames{i},'r');
            fseek(fid,(y1-1)*nx*8,-1);
            tmp=fread(fid,[nx*2,(y2-y1)+1],'real*4');
            fclose(fid);
            tmp=tmp(1:2:end,:)+1j*tmp(2:2:end,:);           
            tmp=tmp(x1:x2,:);
            data(i,:,:)=[x0a;y0a tmp y0b;x0b];
        end
    case 'rmg' %rmg
        data=zeros(n,newdx,newdy*2);
        for i=1:n
            fid=fopen(filenames{i},'r');
            fseek(fid,(y1-1)*nx*8,-1);
            tmp=fread(fid,[nx,(y2-y1)+1],'real*4');
            fclose(fid);
            tmp=tmp(x1:x2,:);
            data(i,:,:)=tmp;
        end
    case 'r4' %r4
        data=zeros(n,newdx,newdy);
        for i=1:n
            fid=fopen(filenames{i},'r');
            fseek(fid,(y1-1)*nx*4,-1);
            tmp=fread(fid,[nx,(y2-y1)+1],'real*4');
            fclose(fid);
            tmp=tmp(x1:x2,:);
            data(i,:,:)=tmp;
        end
end

