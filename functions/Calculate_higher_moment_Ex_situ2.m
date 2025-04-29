function [moment_raw_value, m1_theory] = Calculate_higher_moment_Ex_situ2(x, expdata, order, options)
    % Obtain the moment raw value
    nval = expdata.nval;
    pn = expdata.pntdata;
    nidx = isnan(nval);
    nval = nval(~nidx);
    pn = pn(:,~nidx);
    moment_raw_value = sum(nval.^order.*pn,2);
    % Calculate theoretical moment
    meannumberdata = expdata.smooth_meannumberdata;
    dndtmean = expdata.smooth_dndtmean;
    timedata = expdata.timedata;
    
    indices = options.indices;
    t = expdata.timedata;
    T = expdata.tval;

    lnqfqb = makevec(T, x(indices{1})); %log(qf/qb)
    lnqeqb = makevec(T, x(indices{2})); %log(qe/qb)
    alpha_prime = makevec(T, x(indices{3})); %alpha + 4
    kar1inf = makevec(T, x(indices{4})); % kappa_a*rho_{1,\infty}
    kass_d1 = makevec(T, x(indices{5})); % kappa_a*sigma_s/D_1
    sigma_s = x(indices{6});
    supsat = makevec(T, x(indices{7})); % Supersaturation ratio, \rho_1(t)/rho_{1,\infty}
    shape_index = makevec(T, x(indices{8})); % shape index, \delta

    qfqb = exp(lnqfqb);
    qeqb = exp(lnqeqb);
    shape_factor = si2sf(shape_index);
    
    lent = length(t);
    J_1 = zeros(lent,1);
    knet_1 = zeros(lent,1);
    tidxList = nan(lent, 1);
    for idx = 1:lent
        tidxList(idx) = find(T == t(idx));
    end
    for idx=1:lent
        tidx = tidxList(idx);
        qfqb_ = qfqb(tidx);
        qeqb_ = qeqb(tidx);
        alp_ = alpha_prime(tidx);
        kar1inf_ = kar1inf(tidx);
        kass_d1_ = kass_d1(tidx);
        supsat_ = supsat(tidx);
        si_ = shape_index(tidx);
        sf_ = shape_factor(tidx);

        [Nf, Ne] = Calculate_Nf_Ne(nval, si_);
        [Nf1,Ne1] = Calculate_Nf_Ne(nval-1, si_);
        delta_Nf = Nf - Nf1;
        delta_Ne = Ne - Ne1;
        [kna, knd] = Calculate_kna_knd(nval, qfqb_, qeqb_, alp_, kar1inf_, kass_d1_, sigma_s, sf_, delta_Nf, delta_Ne);
        % (n+1)^(1/3) mean
        kna_1 = (nval+1).^(order).*kna;
        knd_1 = nval.^(order).*knd;
        kna_mean = sum(kna.*pn(idx,:));
        kna_1_mean = sum(kna_1.*pn(idx,:));
        knd_1_mean = sum(knd_1.*pn(idx,:));

        % (n)^(1/3) mean
        kna_1p =  nval.^(order).*kna;
        knd_1p =  (nval-1).^(order).*knd;
        knd_mean = sum(knd.*pn(idx,:));
        kna_1p_mean = sum(kna_1p.*pn(idx,:));
        knd_1p_mean = sum(knd_1p.*pn(idx,:));
        J_1(idx) = kna_1_mean.*supsat_ - knd_1_mean - kna_1p_mean.*supsat_ + knd_1p_mean - kna_mean.*supsat_ + knd_mean ./ meannumberdata(idx);
        knet_1(idx) = kna_mean.*supsat_ - knd_mean;
    end
    dlogmeandata = dndtmean ./ meannumberdata;
    dlntotdata = (knet_1 - dlogmeandata(tidxList))./(meannumberdata(tidxList) - 1);
    dlntot_f = griddedInterpolant(timedata, dlntotdata, 'spline');
    dlntotf = @(t)dlntot_f(t);
    totmdata = zeros(lent,1);
    for idx=1:lent
        totmdata(idx) = exp(integral(dlntotf, timedata(1), timedata(idx)));
    end
    totmdata_g = griddedInterpolant(timedata, totmdata, 'spline');
    moment31_g = griddedInterpolant(timedata, J_1, 'spline');
    relmoment31 = zeros(lent,1);
    for idx=1:lent
        moment31f = @(t)moment31_g(t).*(totmdata_g(t)./totmdata_g(timedata(idx)));
        relmoment31(idx) = 1 + (moment_raw_value(1)./meannumberdata(1) - 1)./totmdata(idx) + ...
            integral(moment31f, timedata(1), timedata(idx));
    end
    m1_theory = relmoment31 .* meannumberdata(tidxList);
end
% changed by soar8nalra@gmail.com