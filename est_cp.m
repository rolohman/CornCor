function [c0,cp,permbad,synth]=est_cp(d,Gi0,diags,mincor,perc)
[~,nint]=size(Gi0);
OPTIONS = optimset('Display','none');
dn=1:nint+1;
intid=nchoosek(dn,2);

c0    = max(d);
%cdiag = d(diags)-c0;
cdiag=d(diags)-c0;

res=d-c0;
for i=1:nint
    a=find(Gi0(:,i));
    cp0(i)=min(0,mymax(res(a),10));
end


[~,sortid]=sort(cp0);
cp=0*cdiag;

figure
triplot(res,dn,intid),colorbar
for i=1:nint
    a=find(Gi0(:,sortid(i)));
    cp(sortid(i))=min(0,mymax(res(a),10));
    cpp=Gi0(:,sortid(i))*cp(sortid(i));
    scaleperc = log(perc+exp(cpp)*(1-perc));
    res=res-scaleperc;
    triplot(res,dn,intid),colorbar
    pause
end
%pull out all intervals with diag coherence of mincor or lower
permbad=cdiag<log(mincor);
synth=Gi0*cp;

%which ints are bad?
%tmp = Gi0*permbad;
% %good = or(tmp==0,diags); %keep the ones on the diagonal for reference
% %ng   = sum(good);
% G0    = [ones(ni,1) Gi0];
% 
% 
% for i=1:ni
%     G(i,:)=G0(i,:)*wgts(i);
% end
% d2=d.*wgts;
% 
% G(~isfinite(G))=0;
% d2(~isfinite(d))=0;
% 
% 
% lb=[-1;cdiag-0.5];
% ub=zeros(nint+1,1);
% mod=lsqlin(G,d2,[],[],[],[],lb,ub,[],OPTIONS);
% 
% c0=mod(1);
% cp=mod(2:end);
% 
% synth=c0+Gi0*cp;
%synth(~good)=NaN;


