function [res,synth]=phsfit(model,m,wphs,means,wgts)

slope=model(1);
inter=model(2);
%inter=0;

synth=slope*m+inter.*(sign(m));

res0=wphs.*exp(1j*means).*conj(exp(1j*synth));
res=angle(res0).*wgts;



