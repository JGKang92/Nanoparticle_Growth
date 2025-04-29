function [pw_growth_rate] = Calculate_pnjn(xfit, expdata, result, options)
    
    indices = options.indices;
    sigma_s = options.sigma_s;

    timedata = expdata.timedata;

    lnqfqb = makevec(timedata, xfit(indices{1})); %ln(q_e/q_b)
    lnqeqb = makevec(timedata, xfit(indices{2})); %ln(q_f/q_b)
    alpha_prime = makevec(timedata, xfit(indices{3})); %\alpha + 4
    kar1inf = makevec(timedata, xfit(indices{4})); % \kappa_a * rho_{1, \infty}
    kass_d1 = makevec(timedata, xfit(indices{5})); % \kappa_a*\sigma_s/D_1
    supsat = makevec(timedata, xfit(indices{7})); %supersaturation ratio rho_1(t)/rho_{1,\infty}
    shape_index = makevec(timedata, xfit(indices{8})); %shape index(\delta)

    shape_factor = si2sf(shape_index);
    
    qfqb = exp(lnqfqb);
    qeqb = exp(lnqeqb);
    timedata = expdata.timedata;
    p_cum_dev = expdata.p_cum_dev;
    lent = length(timedata);
    pw_growth_rate = cell(lent,1);
    for idx=1:lent
        qfqb_ = qfqb(idx);
        qeqb_ = qeqb(idx);
        alpha_ = alpha_prime(idx);
        kar1inf_ = kar1inf(idx);
        kass_d1_ = kass_d1(idx);
        supsat_ = supsat(idx); 
        si_ = shape_index(idx); 
        sf_ = shape_factor(idx);

        ptimedev = p_cum_dev{idx};
        nval = ptimedev(:,1);
        [Nf, Ne] = Calculate_Nf_Ne(nval, si_);
        [Nf1,Ne1] = Calculate_Nf_Ne(nval-1, si_);
        delta_Nf = Nf - Nf1;
        delta_Ne = Ne - Ne1;
        val = (1 + 1./(nval-1)).^alpha_.*(qfqb_.^delta_Nf).*(qeqb_.^delta_Ne);
        % kna = 4.*pi.*sigma_s^2.*kar1inf_.*sf_.*(nval).^(2/3)./(1+kass_d1_.*(nval).^(1/3));
        knd = 4.*pi.*sigma_s^2.*kar1inf_.*sf_.*(nval-1).^(2/3)./(1+kass_d1_.*(nval-1).^(1/3))./val; 
        n_chemical_potential =  -delta_Nf.*log(qfqb_) - delta_Ne.*log(qeqb_) - alpha_.*log(1+1./nval);
        mono_chemical_potential = ones(length(nval),1).*log(supsat_);
        pvalue = result.pnt_guess(idx, :);
        pvalue = pvalue(nval)';
        pw_growth_rate{idx} = [nval, ...
            knd.*pvalue.*(exp(-n_chemical_potential + mono_chemical_potential) - 1)];
            % export_data = [nval, exp(-n_chemical_potential + mono_chemical_potential) - 1, ...
            %    knd.*pvalue.*(exp(-n_chemical_potential + mono_chemical_potential) - 1)];
    end
end