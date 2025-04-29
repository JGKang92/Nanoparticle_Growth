function [Jn, psi_n_val, rhoct] = Calculate_Jn_tilde_Ex_situ2(xfit, expdata, result, options)
    totmdata = result.totmdata;
    mean_theory = result.mean_theory;
    % p_cum_dev = expdata.pcum_dev_time;
    % pcumval = expdata.pcum_time_interp;
    T = expdata.tval;
    t = expdata.timedata;
    dndtmean = expdata.smooth_dndtmean;
    meannumberdata = expdata.smooth_meannumberdata;
    indices = options.indices;
    
    lnqfqb = makevec(T, xfit(indices{1})); %log(qf/qb)
    lnqeqb = makevec(T, xfit(indices{2})); %log(qe/qb)
    alpha_prime = makevec(T, xfit(indices{3})); %alpha + 4
    kar1inf = makevec(T, xfit(indices{4})); % kappa_a*rho_{1,\infty}
    kass_d1 = makevec(T, xfit(indices{5})); % kappa_a*sigma_s/D_1
    sigma_s = xfit(indices{6});
    supsat = makevec(T, xfit(indices{7})); % Supersaturation ratio, \rho_1(t)/rho_{1,\infty}
    shape_index = makevec(T, xfit(indices{8})); % shape index, \delta
    
    qfqb = exp(lnqfqb);
    qeqb = exp(lnqeqb);
    shape_factor = si2sf(shape_index);
    
    totm_meannumber_pp = csapi(T, totmdata./mean_theory);
    totm_meannumber_f = @(t)ppval(totm_meannumber_pp,t);
    d_totm_meannumber_pp = fnder(totm_meannumber_pp,1);
    d_totm_meannumber_f = @(t)fnval(d_totm_meannumber_pp,t);
    
    lent = length(t);
    Jn = cell(lent, 1);
    psi_n_val = cell(lent, 1);
    lnrhoct = zeros(lent, 1);
    for idx=1:lent
        tidx = find(T == t(idx));
        qfqb_ = qfqb(tidx);
        qeqb_ = qeqb(tidx);
        alp_ = alpha_prime(tidx);
        kar1inf_ = kar1inf(tidx);
        kass_d1_ = kass_d1(tidx);
        supsat_ = supsat(tidx);
        si_ = shape_index(tidx);
        sf_ = shape_factor(tidx);

        nval = expdata.nval_c{tidx};
        % Theory of pnt
        [Nf, Ne] = Calculate_Nf_Ne(nval, si_);
        [Nf1,Ne1] = Calculate_Nf_Ne(nval-1, si_);
        delta_Nf = Nf - Nf1;
        delta_Ne = Ne - Ne1;
        [kna, knd] = Calculate_kna_knd(nval, qfqb_, qeqb_, alp_, kar1inf_, kass_d1_, sigma_s, sf_, delta_Nf, delta_Ne);
        knd(nval == 1) = 0;
        pcum_val = expdata.pcum_val_c{tidx};
        pvalue = expdata.pvalue_c{tidx};
        pvalue_1 = expdata.pvalue_1_c{tidx};
        
        psi_n_val_ = d_totm_meannumber_f(t(idx)).*pcum_val + totm_meannumber_f(t(idx)).*expdata.ptimedev_c{tidx};
        psi_n_val{idx} = [nval(:), psi_n_val_(:)];
        
        Jn_ = (kna.*supsat_.*pvalue - knd.*pvalue_1).*(totmdata(idx)./mean_theory(idx));

        Jn{idx} = [nval(:), Jn_(:)];
        small_Jn_mean = sum(kna.*supsat_.*pvalue - knd.*pvalue_1);
        lnrhoct(idx) = (small_Jn_mean - dndtmean(idx)) ./ (meannumberdata(idx) - 1);
    end
    lnrhoct_pp = pchip(t, lnrhoct);
    lnrhoct_fun = @(x)ppval(lnrhoct_pp, x);
    rhoct = zeros(1,lent);
    for idx = 1:lent
        rhoct(idx) = exp(integral(@(x)lnrhoct_fun(x), t(1), t(idx)));
    end
end