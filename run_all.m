define_params
nx=params.nx;
ny=params.ny;
mincor = 0.45;
init_dir(params)
[dates,slcnames]                  = list_slcs(params);
[dn,nd,intid,dt,ni,id1,id2,diags] = define_pairs(dates,params.maxdt,params.dt1);
[Gi0,Gr0]                         = make_G(ni,nd,id1,id2);
filenames                         = make_filenames(params,dates,nd);
fids                              = open_files(filenames,'w');
n=params.rx*params.ry;
shortid=find(dt==1);

load baselines.txt
bpr     = baselines(id1)-baselines(id2);
abpr    = abs(bpr);
bigbase = abpr'>params.minbase;

[windx,windy,wind,windn,wind3,windn3] = make_kernel(params);


smallx=1:params.rx/2:nx;
smally=[params.ry*2 params.ry*2+[1:params.ry/2:params.dely params.dely+1]];
[sy,sx] = meshgrid(smally,smallx);

for j=1:params.dely:ny
    j
    
    cpx = load_slc_chunk(params,slcnames,1,nx,j-floor(length(windy)/2),j+floor(length(windy)/2)+params.dely,'cpx');
    disp('done loading')
    
    [gamma,ints,cors,hp]=make_cor(cpx,intid,wind,windn,wind3,windn3,params);
    disp('done making cor')
    
    %invert cor on coarse grid
    for k=1:length(smally)
        
        for i=1:length(smallx)
            data  = squeeze(cors(:,smallx(i),smally(k)))';
            d     = log(data);
            stds  = sqrt(-2*d);
            
            [c00,cp0,mk0,synth0] = invert_m_mat(d',nd,Gi0,intid);
            mk0(mk0==0)=NaN;
            mk0 = mk0+mymax(-mk0,3);
            mk0(mk0<0)=0;         
            mk0(isnan(mk0))      = 0;
           
            d2                   = d'-c00-Gi0*cp0;
            d2(isnan(synth0))    = NaN;
            
            good    = isfinite(d2);
            start   = [mk0];
            LB      = [zeros(nd,1)];
            UB      = [inf(nd,1)];
            mk1     = lsqnonlin('corfit',start,LB,UB,params.OPTIONS,d2(good),stds(good)',Gr0(good,:),1);
            
            mk1         = mk1+ mymax(-mk1,3);
            mk1(mk1<0)  = 0;
                
            d2 = d-c00-Gi0*cp0;
            [~,synth] = corfit(mk1,d2,stds',Gr0,1);
            synth = synth+c00+Gi0*cp0;
            res   = d'-synth;
            good  = good & res<0.2 & synth>log(mincor) & d'>log(mincor);
            
            small_mcor_orig(i,k) = sqrt(mean(d(good).^2,'omitnan'));
            small_mcor_res(i,k)  = sqrt(mean(res(good).^2,'omitnan'));
            
            smallc0(i,k)=c00;
            smallcp(i,k,:)=cp0;
            smallmk1(i,k,:)=mk1;
        end
    end
    
    for k=params.ry*2+[1:params.dely]
        
        for i=1:nx
            data  = squeeze(cors(:,i,k))';
            d     = log(data);
            
            GINT        = griddedInterpolant(sx,sy,squeeze(smallmk1(:,:,1)));
            for l=1:nd
                GINT.Values = squeeze(smallmk1(:,:,l));
                mk1(l) = GINT(i,k);
            end
            for l=1:nd-1
                GINT.Values = squeeze(smallcp(:,:,l));
                cp(i,l) = GINT(i,k);
            end
            GINT.Values = smallc0;
            c0(i)   = GINT(i,k);
            GINT.Values = small_mcor_orig;
            mcor_orig(i)=GINT(i,k);
            GINT.Values = small_mcor_res;
            mcor_res(i)=GINT(i,k);
            allmk1(i,:)=mk1;
            meank1       = atan(mk1);
            
            
            
            
            cdiag = d(diags)-c0(i);
            
            %pull out all intervals with diag coherence of mincor or lower
            permbad=cdiag<log(mincor);
            
            %which ints are bad?
            tmp = Gi0*permbad';
            good = or(tmp==0 & d'>log(mincor),dt==1); %keep the ones on the diagonal for reference
            
            %             synth = c0+Gi0*cp0;
            %             res   = d'-synth;
            %             good  = good & res<0.2 & synth>log(mincor) & d'>log(mincor);
            %             [res,synth,dc0] = corfit(mk1,res(good),stds(good)',Gr0(good,:),1);
            %             synth           = synth+c00+dc0+Gi0(good,:)*cp0;
            
            %now do hp fit
            aph    = squeeze(hp(:,i,k));
            
            ph     = transpose(invert_phsmat(aph(good),nd,intid(good,:)));
            ph(isnan(ph))=1;
            aph2   = ph(intid(:,2)).*conj(ph(intid(:,1)));
            res    = aph.*conj(aph2);
            res    = res./abs(res);
            tphs_res(i)=abs(mean(res(good)));
            
            
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
        write_output(fids,mcor_orig,mcor_res,c0,allmk1,cp,slopesk1,tphs_orig,tphs_res,tphs_resk,allhp,alllp);
    end
end


fclose('all');


