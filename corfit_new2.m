function [res,synth,c0,cpp,scaleperc]=corfit_new2(model,cors,stds,Gi0)
%expects logs
nint=size(Gi0,2);
c0=model(1);
cp=model(1+[1:nint]);
perc=model(end);



cpp=Gi0*cp;
scaleperc = log(perc+exp(cpp)*(1-perc));


synth=c0+scaleperc;


res=exp(cors)-exp(synth);
if(sum(res>0.1)>100)
    res=1000*ones(size(res));
else

res=res./(stds);
end

