function [Jn, psi_n_val, rhoct, ratio_comparison] = Calculate_Jn_tilde_In_situ(xfit, expdata, result, options)
    % ================================================== Part 1. parsing input arguments ==================================================
    totmdata = result.totmdata;
    mean_theory = result.mean_theory;
    p_cum_dev = expdata.p_cum_dev;
    pcumval = expdata.pcumval;
    timedata = expdata.timedata;
    dndtmean = expdata.smooth_dndtmean;
    meannumberdata = expdata.smooth_meannumberdata;
    indices = options.indices;
    sigma_s = options.sigma_s;

    lent = length(timedata);
    lnqfqb = makevec(timedata, xfit(indices{1})); %log(qf/qb)
    lnqeqb = makevec(timedata, xfit(indices{2})); %log(qe/qb)
    alpha_prime = makevec(timedata, xfit(indices{3})); %alpha + 4
    kar1inf = makevec(timedata, xfit(indices{4})); % kappa_a*rho_{1,\infty}
    kass_d1 = makevec(timedata, xfit(indices{5})); % kappa_a*sigma_s/D_1
    supsat = makevec(timedata, xfit(indices{7})); % Supersaturation ratio, \rho_1(t)/rho_{1,\infty}
    shape_index = makevec(timedata, xfit(indices{8})); % shape index, \delta

    qfqb = exp(lnqfqb);
    qeqb = exp(lnqeqb);
    shape_factor = si2sf(shape_index);
    
    % ================================================== Part 1. parsing input arguments ==================================================

    totm_meannumber_pp = csapi(timedata, totmdata./mean_theory);
    totm_meannumber_f = @(t)ppval(totm_meannumber_pp,t);
    d_totm_meannumber_pp = fnder(totm_meannumber_pp,1);
    d_totm_meannumber_f = @(t)fnval(d_totm_meannumber_pp,t);
    
    Jn = cell(lent, 1);
    psi_n_val = cell(lent, 1);
    lnrhoct = zeros(lent, 1);
    ratio_comparison = cell(lent, 1);
    for idx=1:lent
        qfqb_ = qfqb(idx);
        qeqb_ = qeqb(idx);
        alp_ = alpha_prime(idx);
        kar1inf_ = kar1inf(idx);
        kass_d1_ = kass_d1(idx);
        supsat_ = supsat(idx);
        si_ = shape_index(idx);
        sf_ = shape_factor(idx);

        p_time_dev = p_cum_dev{idx};
        nval = p_time_dev(:,1);
        % Theory of pnt
        [Nf, Ne] = Calculate_Nf_Ne(nval, si_);
        [Nf1,Ne1] = Calculate_Nf_Ne(nval-1, si_);
        delta_Nf = Nf - Nf1;
        delta_Ne = Ne - Ne1;
        [kna, knd] = Calculate_kna_knd(nval, qfqb_, qeqb_, alp_, kar1inf_, kass_d1_, sigma_s, sf_, delta_Nf, delta_Ne);
        if isfield(result, 'pnt_guess')
            pvalue = result.pnt_guess(idx, nval)';
            pvalue_1 = result.pnt_guess(idx, nval - 1)';
        else
            pcumval_idx = pcumval{idx};
            idx2 = find(pcumval_idx(:,1) == nval(1));
            pvalue = pcumval_idx(idx2:idx2+length(nval),2);
            pcum_val = pvalue(1:end-1);
            if idx2 ~= 1
                pvalue_1 = pcumval_idx(idx2-1:idx2+length(nval)-1,2);
            else
                pvalue_1 = pcumval_idx(idx2:idx2+length(nval)-1,2);
                pvalue_1 = [1; pvalue_1];
            end
            pvalue = -diff(pvalue);
            pvalue_1 = -diff(pvalue_1);
            psi_n_val{idx} = [nval, d_totm_meannumber_f(timedata(idx)).*pcum_val ...
                                + totm_meannumber_f(timedata(idx)).*p_time_dev(:,2)];
        end
        pvalue(pvalue<0) = 0;
        logval = pvalue(2:end) ./ pvalue(1:end-1);
        logval(isnan(logval)) = 1;
        logval(isinf(logval)) = 1;
        logval(logval==0) = 1;
        Jn{idx} = [nval, (kna.*supsat_.*pvalue - knd.*pvalue_1).*(totmdata(idx)./mean_theory(idx))];
        small_Jn_mean = sum(kna.*supsat_.*pvalue - knd.*pvalue_1);
        lnrhoct(idx) = (small_Jn_mean - dndtmean(idx)) ./ (meannumberdata(idx) - 1);
        ratio_comparison{idx} = [nval(1:end-1), ...
            exp((-delta_Nf(1:end-1).*log(qfqb_) - delta_Ne(1:end-1).*log(qeqb_) - alp_.*log(1+1./nval(1:end-1)))) - 1, ...
            exp((-delta_Nf(1:end-1).*log(qfqb_) - delta_Ne(1:end-1).*log(qeqb_) - alp_.*log(1+1./nval(1:end-1)) + log(logval))) - 1];
    end
    lnrhoct_pp = pchip(timedata, lnrhoct);
    lnrhoct_fun = @(t)ppval(lnrhoct_pp, t);
    rhoct = zeros(1,lent);
    for idx = 1:lent
        rhoct(idx) = exp(integral(@(t)lnrhoct_fun(t), timedata(1), timedata(idx)));
    end
end