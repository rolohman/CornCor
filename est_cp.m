function [c0,cp,permbad,synth]=est_cp(d,Gi0,diags,mincor)
[~,nint]=size(Gi0);
OPTIONS = optimset('Display','none');

c0    = max(d);
%cdiag = d(diags)-c0;
cdiag=d(diags);

[~,sortid]=sort(cdiag);
cp=0*cdiag;

res=d;
for i=1:nint
    a=find(Gi0(:,sortid(i)));
    cp(sortid(i))=min(0,mymax(res(a),10));
    res=res-Gi0(:,sortid(i))*cp(sortid(i));
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


