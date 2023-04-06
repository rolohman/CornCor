define_params
nx=params.nx;
ny=params.ny;

init_dir(params)

[dates,slcnames]                  = list_slcs(params);
[dn,nd,intid,dt,ni,id1,id2,diags] = define_pairs(dates);
[Gi0,Gr0]                         = make_G(ni,nd,id1,id2);
fids                              = open_files(params,dates,nd);


load baselines.txt
bpr     = baselines(id1)-baselines(id2);
abpr    = abs(bpr);
bigbase = abpr'>params.minbase;

[windx,windy,wind,windn,wind3,windn3] = make_kernel(params);

for j=1:params.dely:ny
    j
   
    cpx = load_slc_chunk(params,slcnames,1,nx,j-params.ry,j+params.ry+params.dely,'cpx');
    disp('done loading')

    [gamma,ints,cors,hp]=make_cor(cpx,intid,wind,windn,wind3,windn3,params);
    disp('done making cor')
    
    %add scaling/saving 2d histogram
  
    for k=1:params.dely      
        for i=1:nx
            
            
            data  = squeeze(cors(:,i,k))';
            d     = log(data);
            
            g       = and(isfinite(d),bigbase);
            wgts    = diag(exp(params.bpw*d(g)));
            Gbase   = [abpr(g).^2 ones(size(abpr(g)))];
            Ggbase  = inv(Gbase'*wgts*Gbase)*Gbase'*wgts;       
            basemod = Ggbase*d(g)';            
            bps     = min(basemod(1),0);
            d       = d-bps*abpr'.^2;
            
            [c00,cp0,cr0] = invert_cormat(d',nd,Gi0);
            cr0(cr0>0)    = 0;
            synth         = c00+Gi0*cp0'-abs(Gr0*cr0');
            
            mk1_init = sqrt(1./exp(cr0).^2-1);
            mk2_init = 2*sqrt(1./exp(cr0)-1);
            
            return
            start   = [c00;cp0';mk1_init(2:end)'];
            LB      = [-inf(nd,1);zeros(nd-1,1)];
            UB      = [zeros(nd,1);inf(nd-1,1)];
            modk1   = lsqnonlin('corfit',start,LB,UB,params.OPTIONS,d',Gr0,Gi0,nd,1);
            start   = [c00;cp0';mk2_init(2:end)'];
            modk2   = lsqnonlin('corfit',start,LB,UB,params.OPTIONS,d',Gr0,Gi0,nd,2);
            allc00(i,k)=c00;
            allc0(i,k)=modk1(1);
            allcp(i,:,k)=modk1(2:nd);
            mk1 = [0;modk1(nd+1:end)];
            mk2 = [0;modk2(nd+1:end)];

  
            allmk1(i,:,k)=mk1;
            allmk2(i,:,k)=mk2;
            
            mk1         = squeeze(amk1(i,k,:));
            mk2         = squeeze(amk2(i,k,:));
            
            meank1      = atan(mk1);
            meank2      = atan(mk2./(1-0.25*mk2.^2));
            neg         = 1<0.25*mk2.^2;
            meank2(neg) = meank2(neg)+pi;
            

            
            %now do hp fit
            aph    = squeeze(hp(:,i,k));
            ph     = transpose(invert_phsmat(aph,nd));
            pa1    = angle(ph.*exp(1j*meank1)); %for calculation of R
            pa2    = angle(ph.*exp(1j*meank2)); %for calculation of R
            
            ts     = linspace(-8,8,100);
            gid    = and(mk1>0.4,isfinite(angle(ph)));
            gi2    = gi; %not storms
         
            if(sum(gid)>1)
                synths = mk1*ts;
                dat    = ph.*exp(1j*meank1);
                res    = dat.*conj(exp(1j*synths));
                fit    = sqrt(mean(angle(res(gid,:)).^2,1));
                besti  = find(fit==min(fit));
                besti  = round(mean(besti));
                slopes0k1(i,k)=ts(besti);
                
                synths = mk2*ts;
                dat    = ph.*exp(1j*meank2);
                res    = dat.*conj(exp(1j*synths));
                fit    = sqrt(mean(angle(res(gid,:)).^2,1));
                besti  = find(fit==min(fit));
                besti  = round(mean(besti));
                slopes0k2(i,k)=ts(besti);
            elseif(sum(gid)==1)
                slopes0k1(i,k)=angle(ph(gid)*exp(1j*meank1(gid)))/mk1(gid);
                slopes0k2(i,k)=angle(ph(gid)*exp(1j*meank2(gid)))/mk2(gid);
            else
                slopes0k1(i,k)=0;
                slopes0k2(i,k)=0;
            end
  
          [m1a,r1]  = lsqnonlin('phsfit_noint',slopes0k1(i,k),[],[],OPTIONS,mk1(gi2),ph(gi2),meank1(gi2),mk1(gi2));
          [m1b,r2]  = lsqnonlin('phsfit_noint',0*slopes0k1(i,k),[],[],OPTIONS,mk1(gi2),ph(gi2),meank1(gi2),mk1(gi2));
          if(r1>r2)
                  slopesk1(i,k)=m1b;
          else
                  slopesk1(i,k)=m1a;
          end
          resk1 = phsfit_noint(slopesk1(i,k),mk1,ph,meank1,ones(size(mk1)));

     %     [m2a,r1]  = lsqnonlin('phsfit_noint',slopes0k2(i,k),[],[],OPTIONS,mk2(gi2),ph(gi2),meank2(gi2),mk2(gi2));
     %     [m2b,r2]  = lsqnonlin('phsfit_noint',0*slopes0k2(i,k),[],[],OPTIONS,mk2(gi2),ph(gi2),meank2(gi2),mk2(gi2));
     %     if(r1>r2)
     %             slopesk2(i,k)=m2b;
     %     else
     %     
%	  	  slopesk2(i,k)=m2a;
%          end 
%	  resk2 = phsfit_noint(slopesk2(i,k),mk2,ph,meank2,ones(size(mk2)));

            ocor(i,k)       = abs(mean(ph));
            cork1(i,k)      = abs(mean(exp(1j*resk1)));
            cork2(i,k)      = abs(mean(exp(1j*resk2)));
        
            ocorb(i,k)      = abs(mean(ph(gid)));
            corbk1(i,k)     = abs(mean(exp(1j*resk1(gid))));
            corbk2(i,k)     = abs(mean(exp(1j*resk2(gid))));
n=sum(gid);
if(n>2)
	a1=mean(angle(ph));
	a2=mean(angle(resk1));
  R(i,k)=1-(sum((angle(ph)-a1).^2)/sum((angle(resk1)-a2).^2));
Rb(i,k)=(1-(1-R(i,k))*(n-1)/(n-2));
else
end
    end  
        end
    write_output(fids,stuff)
    
        end



fclose('all');


