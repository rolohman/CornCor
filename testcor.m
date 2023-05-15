define_params
nx=params.nx;
ny=params.ny;

init_dir(params)
[dates,slcnames]                  = list_slcs(params);
[dn,nd,intid,dt,ni,id1,id2,diags] = define_pairs(dates,params.maxdt,params.dt1);
[Gi0,Gr0]                         = make_G(ni,nd,id1,id2);


data  = abs(gamma);
d     = log(data)';


mcor_orig=sqrt(mean(d.^2,'omitnan'));

[c00,cp0,mk0,synth0] = invert_m_mat(d',nd,Gi0,intid);

d       = d-[Gi0*cp0']';
start   = [c00;cp0';mk0(2:end)];
LB      = [-inf(nd,1);zeros(nd-1,1)];
UB      = [zeros(nd,1);inf(nd-1,1)];
modk1   = lsqnonlin('corfit',start,LB,UB,params.OPTIONS,d',Gr0,Gi0,nd,1);
[res,synth]      = corfit(modk1,d',Gr0,Gi0,nd,1);
mcor_res=sqrt(mean(res.^2,'omitnan'));


