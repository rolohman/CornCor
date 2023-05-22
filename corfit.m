function [res,synth,dc0]=corfit(model,cors,stds,Gr0,k)
%expects logc0/cp cors;
%c0=model(1);
%cp=model(2:nd);
%m=[0;model(nd+1:end)];
m=model;


mdiff=Gr0*m;
if(k==1)
    gamma=1./sqrt(mdiff.^2+1);
else
    gamma=1./(mdiff.^2/4+1);
end
gamma(isinf(gamma))=0;
gamma=log(gamma);
%cpp=Gi0*cp;
%synth=c0+gamma+cpp;
synth=gamma;
res=cors-synth;
dc0=median(res);
res=res-dc0;
%cors(cors>0)=0;
%stds=sqrt(-2*cors);
res=res./(stds);

