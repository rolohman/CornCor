function [gamma,ints,cors,hp]=make_cor(cpx,intid,wind,windn,wind3,windn3,params)

ni     = size(intid,1);
nx     = params.nx;
dely   = params.dely;
ry     = params.ry;

gamma  = nan(ni,nx,dely);
ints   = nan(ni,nx,dely);
amps   = cpx.*conj(cpx);
ampsum = sqrt(convn(amps,wind3,'same')./windn3);

for i=1:ni
    slc1 = squeeze(cpx(intid(i,1),:,:));
    slc2 = squeeze(cpx(intid(i,2),:,:));
    
    c    = slc2.*conj(slc1);
    asum = squeeze(ampsum(intid(i,1),:,:));
    bsum = squeeze(ampsum(intid(i,2),:,:));
    csum = conv2(c,wind,'same')./windn;
    cpx3 = csum./asum./bsum;
    cpx3(isnan(cpx3)) = 0;
    
    gamma(i,:,:) = cpx3(:,ry+[1:dely]);
    ints(i,:,:)  = c(:,ry+[1:dely]);
end
cors          = abs(gamma);
good          = cors>0;
cors(cors>1)  = 1;
cors(~good)   = NaN;
hp            = ints.*conj(gamma);
hp            = hp./abs(hp);
hp(isnan(hp)) = 1;

ngood         = squeeze(sum(good,1));
goodid        = ngood>=10;

