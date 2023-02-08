function res=corfit(model,cors,Gr0,Gi0,nd,k)
%expects logc0/cp cors;
c0=model(1);
cp=model(2:nd);
m=[0;model(nd+1:end)];

mdiff=Gr0*m;
if(k==1)
    gamma=1./sqrt(mdiff.^2+1);
else
    gamma=1./(mdiff.^2/4+1);
end
gamma(isinf(gamma))=0;
gamma=log(gamma);
cpp=Gi0*cp;
synth=c0+cpp+gamma;
res=cors-synth;


