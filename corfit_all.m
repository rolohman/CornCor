function [res,synth]=corfit(model,cors,stds,Gi0,Gr0,k)
%expects logc0/cp cors;
nd=size(Gr0,2);
%c0=model(1);
cp=model(1:nd-1);
m=[model(nd:end)];

c0=max(cors);

mdiff=Gr0*m;
if(k==1)
    gamma=1./sqrt(mdiff.^2+1);
else
    gamma=1./(mdiff.^2/4+1);
end
gamma(isinf(gamma))=0;
gamma=log(gamma);
cpp=Gi0*cp;
synth=c0+gamma+cpp;
res=exp(cors')-exp(synth);

res=res./(stds);

