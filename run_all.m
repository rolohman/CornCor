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

dt2=dn(intid(:,2))-dn(intid(:,1));
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
            stds  = sqrt(-2*d');
            
            ds=[1:12:100 logspace(log10(110),log10(max(dt2)),25)];
            bins=ds(1:end-1)+diff(ds)/2;
            for l=1:length(bins)
                id=dt2>=ds(l) & dt2<ds(l+1);
                maxt(l)=mymax(d(id),100);
                nbin(l)=sum(id);
            end
            Gbin=[ones(length(bins),1) sqrt(bins)' ];
            Ginv=inv(Gbin'*Gbin)*Gbin';
            mod=Ginv*maxt';
            c00=mod(1);
            pf=mod(2);
           
            synth00=c00+pf*sqrt(dt2);
            d2=d-synth00';
            
            [c0,cp0,mk0,synth0] = invert_m_mat(d2',stds,nd,Gi0,intid);
            mk0(mk0==0)=NaN;
            mk0 = mk0+mymax(-mk0,3);
            mk0(mk0<0)=0;         
            mk0(isnan(mk0))      = 0;
           
           % d2                   = d'-c00-Gi0*cp0;
            %d2(isnan(synth0))    = NaN;
            
            good    = isfinite(synth0);
            start   = [c0;d(diags)'-c0;mk0];
            start=start(2:end);
            LB      = [-inf(nd-1,1);zeros(nd,1)];
            UB      = [zeros(nd-1,1);inf(nd,1)];
            %start   = [mk0];
            %            LB      = [zeros(nd,1)];
            %            UB      = [inf(nd,1)];
            mod     = lsqnonlin('corfit_all',start,LB,UB,params.OPTIONS,d2(good),stds(good),Gi0(good,:),Gr0(good,:),1);
            [res1,synth]= corfit_all(mod,d2,stds,Gi0,Gr0,1);
            
            c0=c00+c0;
            cp=mod(1:nd-1);
            mk1=mod(nd:end);
            mk1         = mk1+ mymax(-mk1,5);
            mk1(mk1<0)  = -mk1(mk1<0);
           
            res   = d2'-synth;
            %good  = good & res<0.2 & synth>log(mincor) & d'>log(mincor);
            
            small_mcor_orig(i,k) = sqrt(mean(d(good).^2,'omitnan'));
            small_mcor_res(i,k)  = sqrt(mean(res(good).^2,'omitnan'));
            
            smallc0(i,k)=c0;
            smallcp(i,k,:)=cp;
            smallmk1(i,k,:)=mk1;
            smallpf(i,k,:)=pf;
        end
    end
    
    for k=params.ry*2+[1:params.dely]
        
        for i=1:nx
            i
            
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
            GINT.Values = smallpf;
            pf(i)=GINT(i,k);
            
            allmk1(i,:)=mk1;
            meank1       = atan(mk1);
            
            
            
 
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
        write_output(fids,mcor_orig,mcor_res,c0,pf,allmk1,cp,slopesk1,tphs_orig,tphs_res,tphs_resk,allhp,alllp);
    end
end


fclose('all');


