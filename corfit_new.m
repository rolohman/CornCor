function [res,synth,c0,gamma,cpp,scaleperc]=corfit_new(model,cors,stds,Gi0,Gr0)
%expects logs
nd=size(Gr0,2);
c0=model(1);
cp=model(1+[1:nd-1]);
m=model(nd+[1:nd]);
perc=model(nd*2+1);




mdiff=Gr0*m;
gamma=1./sqrt(mdiff.^2+1);



gamma=log(gamma);


cpp=Gi0*cp;
scaleperc = log(perc+exp(cpp)*(1-perc));


synth=c0+gamma+scaleperc;


res=exp(cors)-exp(synth);


res=res./(stds);

