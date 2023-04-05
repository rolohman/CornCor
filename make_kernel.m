function[windx,windy,wind,windn,wind3,windn3]=make_kernel(params)

rx     = params.rx;
ry     = params.ry;
kernel = params.kernel;
dely   = params.dely;
nx     = params.nx;

switch kernel
    case 'Gaussian'
        windx        = exp(-(-rx*2:rx*2).^2/2/(rx/2)^2);
        windy        = exp(-(-ry*2:ry*2).^2/2/(ry/2)^2);

    case 'Uniform'
        windx        = ones(1,rx);
        windy        = ones(1,ry);

    case 'Triangle'
        windx        = linspace(0,1,rx);
        windy        = linspace(0,1,ry);
        windx        = [windx fliplr(windx(1:end-1))];
        windy        = [windy fliplr(windy(1:end01))];

    otherwise
        disp('kernel type must be Gaussian, Uniform or Triangle')
        return
end
        
windx         = windx/sum(windx);
windy         = windy/sum(windy);
wind          = windy'*windx;


windn         = conv2(ones(nx,dely+ry*2+1),wind,'same');

windn3(1,:,:) = windn;
wind3(1,:,:)  = wind;
