function [c0,cp,permbad,synth]=est_cp(d,wgts,Gi0,diags,mincor)
[~,nint]=size(Gi0);
OPTIONS = optimset('Display','none');


c0    = max(d);
cdiag = d(diags)-c0;

%pull out all intervals with diag coherence of mincor or lower
permbad=cdiag<log(mincor);

%which ints are bad?
tmp = Gi0*permbad;
good = or(tmp==0,diags); %keep the ones on the diagonal for reference
ng   = sum(good);
G    = [ones(ng,1) Gi0(good,:)];

for i=1:ng
    G(i,:)=G(i,:)*wgts(i);
end
d2=d.*wgts;

G(~isfinite(G))=0;
d2(~isfinite(d))=0;


lb=[c0*2;cdiag-0.5];
ub=zeros(nint+1,1);
mod=lsqlin(G,d2(good),[],[],[],[],lb,ub,[],OPTIONS);

c0=mod(1);
cp=mod(2:end);

synth=c0+Gi0*cp;
synth(~good)=NaN;


