function [res,synth]=phsfit_noint(model,m,wphs,means,wgts)

slope=model(1);



synth=slope*m;
res0=wphs.*exp(1j*means).*conj(exp(1j*synth));
res=angle(res0).*wgts;



