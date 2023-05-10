define_params
nx=params.nx;
ny=params.ny;

init_dir(params)
[dates,slcnames]                  = list_slcs(params);
[dn,nd,intid,dt,ni,id1,id2,diags] = define_pairs(dates,params.maxdt,params.dt1);
[Gi0,Gr0]                         = make_G(ni,nd,id1,id2);
filenames                         = make_filenames(params,dates,nd);
fids                              = open_files(filenames,'w');

shortid=find(dt==1);

load baselines.txt
bpr     = baselines(id1)-baselines(id2);
abpr    = abs(bpr);
bigbase = abpr'>params.minbase;

[windx,windy,wind,windn,wind3,windn3] = make_kernel(params);
%j=133;
for j=1:params.dely:ny
    j
    
    cpx = load_slc_chunk(params,slcnames,1,nx,j-floor(length(windy)/2),j+floor(length(windy)/2)+params.dely,'cpx');
    disp('done loading')
    
    [gamma,ints,cors,hp]=make_cor(cpx,intid,wind,windn,wind3,windn3,params);
    disp('done making cor')
    
    %add scaling/saving 2d histogram

    for k=1:params.dely
        for i=1:nx
       %i=361
      
            data  = squeeze(cors(:,i,k))';
            d     = log(data);
            
            g       = and(isfinite(d),bigbase);
            wgts    = diag(exp(params.bpw*d(g)));
            Gbase   = [abpr(g).^2 ones(size(abpr(g)))];
            Ggbase  = inv(Gbase'*wgts*Gbase)*Gbase'*wgts;
            basemod = Ggbase*d(g)';
            bps     = min(basemod(1),0);
            allbps(i)=bps;
            d       = d-bps*abpr'.^2;
   
            mcor_orig(i)=sqrt(mean(d.^2,'omitnan'));
             
            %[c00,cp0,cr0] = invert_cormat(d',nd,Gi0,intid);
            %cr0(cr0>0)    = 0;
            %d=d-[Gi0*cp0']';

            %mk1_init = sqrt(1./exp(cr0).^2-1);
            %mk2_init = 2*sqrt(1./exp(cr0)-1);
             
            [c00,cp0,mk0] = invert_m_mat(d',nd,Gi0,intid);
            
            d       = d-[Gi0*cp0']';
            start   = [c00;mk0(2:end)];
            LB      = [-inf(1);zeros(nd-1,1)];
            UB      = [zeros(1,1);inf(nd-1,1)];
            modk1   = lsqnonlin('corfit',start,LB,UB,params.OPTIONS,d',Gr0,Gi0,nd,1);
            res      = corfit(modk1,d',Gr0,Gi0,nd,1);
            mcor_res(i)=sqrt(mean(res.^2,'omitnan'));
            
            mk1     = [0;modk1(2:end)];
            allc0(i)   = modk1(1);
            allcp(i,:) = cp0;
            allmk1(i,:) = mk1;
            %allcr(i,:)  = cr0;
            meank1      = atan(mk1);
            meank1 = -meank1;
            %now do hp fit
            aph    = squeeze(hp(:,i,k));
            
            ph     = transpose(invert_phsmat(aph,nd,intid));
            aph2   = ph(intid(:,2)).*conj(ph(intid(:,1)));
            res    = aph.*conj(aph2);
            res=res./abs(res);
            tphs_res(i)=abs(mean(res(:)));
            
            
            allhp(i,:)=ph;
            alllp(i,:)=[0;squeeze(gamma(shortid,i,k))];
            
            ts     = linspace(-8,8,100);
            gid    = and(mk1>std(mk1),isfinite(angle(ph)));
            
            %gi2    = gi; %not storms
            gi2=1:nd;
            
            ph=ph./abs(ph);
            ph(~isfinite(ph))=1;
            if(sum(gid)>1)
                synths = mk1*ts;
                dat    = ph.*exp(1j*meank1);
                res    = dat.*conj(exp(1j*synths));
                fit    = sqrt(mean(angle(res(gid,:)).^2,1));
                besti  = find(fit==min(fit));
                besti  = round(mean(besti));
                slopes0k1=ts(besti);
                
            elseif(sum(gid)==1)
                slopes0k1=angle(ph(gid)*exp(1j*meank1(gid)))/mk1(gid);
            else
                slopes0k1=0;
            end
            
            [m1a,r1]  = lsqnonlin('phsfit_noint',slopes0k1,[],[],params.OPTIONS,mk1(gi2),ph(gi2),meank1(gi2),atan(mk1(gi2)));
            [m1b,r2]  = lsqnonlin('phsfit_noint',0*slopes0k1,[],[],params.OPTIONS,mk1(gi2),ph(gi2),meank1(gi2),atan(mk1(gi2)));
            if(r1>r2)
                slopesk1(i)=m1b;
            else
                slopesk1(i)=m1a;
            end
            resk1 = phsfit_noint(slopesk1(i),mk1,ph,meank1,ones(size(mk1)));
            
            
            tphs_orig(i)     = abs(mean(ph));
            tphs_resk(i)      = abs(mean(exp(1j*resk1)));
            
            
        end
        write_output(fids,allbps,mcor_orig,mcor_res,allc0,allmk1,allcp,slopesk1,tphs_orig,tphs_res,tphs_resk,allhp,alllp);
    end
end


fclose('all');


