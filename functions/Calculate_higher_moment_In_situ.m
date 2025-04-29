function [moment_raw_value, m1_theory] = Calculate_higher_moment_In_situ(x, expdata, order, options)
    % ================================================== Part 1. parsing input arguments ==================================================
    % Obtain the moment raw value
    T = expdata.timedata;

    indices = options.indices;
    numberdata = expdata.numberdata;
    moment_raw_value = mean(numberdata.^order,2, 'omitnan');
    % Calculate theoretical moment
    meannumberdata = expdata.smooth_meannumberdata;
    numberdata = expdata.smooth_numberdata;
    moments_t0 = mean(numberdata.^order, 2, 'omitnan');
    dndtmean = expdata.smooth_dndtmean;
    
    sigma_s = options.sigma_s;

    lnqfqb = makevec(T, x(indices{1})); %log(qf/qb)
    lnqeqb = makevec(T, x(indices{2})); %log(qe/qb)
    alpha_prime = makevec(T, x(indices{3})); %alpha + 4
    kar1inf = makevec(T, x(indices{4})); % kappa_a*rho_{1,\infty}
    kass_d1 = makevec(T, x(indices{5})); % kappa_a*sigma_s/D_1
    supsat = makevec(T, x(indices{7})); % Supersaturation ratio, \rho_1(t)/rho_{1,\infty}
    shape_index = makevec(T, x(indices{8})); % shape index, \delta

    qfqb = exp(lnqfqb);
    qeqb = exp(lnqeqb);
    shape_factor = si2sf(shape_index);
    % ================================================== Part 1. parsing input arguments ==================================================
    lent = length(T);
    J_1 = zeros(lent,1);
    knet_1 = zeros(lent,1);
    for idx=1:lent
        qfqb_ = qfqb(idx);
        qeqb_ = qeqb(idx);
        alp_ = alpha_prime(idx);
        kar1inf_ = kar1inf(idx);
        kass_d1_ = kass_d1(idx);
        supsat_ = supsat(idx);
        si_ = shape_index(idx);
        sf_ = shape_factor(idx);

        nval = numberdata(idx,:);
        nidx = ~isnan(nval);
        nval = nval(nidx);
        [Nf, Ne] = Calculate_Nf_Ne(nval, si_);
        [Nf1,Ne1] = Calculate_Nf_Ne(nval-1, si_);
        delta_Nf = Nf - Nf1;
        delta_Ne = Ne - Ne1;
        [kna, knd] = Calculate_kna_knd(nval, qfqb_, qeqb_, alp_, kar1inf_, kass_d1_, sigma_s, sf_, delta_Nf, delta_Ne);
        % (n+1)^(1/3) mean
        kna_1 = (nval+1).^(order).*kna;
        knd_1 = nval.^(order).*knd;
        % (n)^(1/3) mean
        kna_1p =  nval.^(order).*kna;
        knd_1p =  (nval-1).^(order).*knd;
        J_1(idx) = mean(kna_1.*supsat_ - knd_1 - kna_1p.*supsat_ + knd_1p - kna.*supsat_ + knd) ./ meannumberdata(idx);
        knet_1(idx) = mean(kna.*supsat_ - knd);
    end
    dlogmeandata = dndtmean ./ meannumberdata;
    dlntotdata = (knet_1 - dlogmeandata)./(meannumberdata - 1);
    dlntot_f = griddedInterpolant(T, dlntotdata, 'spline');
    
    T0 = T(1);
    totmdata = zeros(lent,1);
    for idx = 1:lent
        totmdata(idx) = exp(integral(@(t)dlntot_f(t), T0, T(idx)));
    end
    totmdata_g = griddedInterpolant(T, totmdata, 'spline');
    moment31_g = griddedInterpolant(T, J_1, 'spline');
    relmoment31 = zeros(lent,1);
    for idx = 1:lent
        tt = T(idx);
        moment31f = @(t)moment31_g(t).*(totmdata_g(t)./totmdata_g(tt));
        relmoment31(idx) = 1 + (moments_t0(1)./meannumberdata(1) - 1)./totmdata(idx) + ...
            integral(moment31f, T0, tt);
    end
    m1_theory = relmoment31 .* meannumberdata;
end
% changed by soar8nalra@gmail.com