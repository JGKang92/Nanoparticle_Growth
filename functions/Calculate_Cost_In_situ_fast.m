function [cost, costs, dtlnrhoct, prediction] = Calculate_Cost_In_situ_fast(x, expdata, options)
    % ================================================== Part 1. parsing input arguments ==================================================
    % See the Eqs. (S5-3), (S5-4), and (S5-5) for details
    t = expdata.timedata(:)';
    meannumberdata = expdata.smooth_meannumberdata;
    numberdata = expdata.smooth_numberdata;
    dndt = expdata.smooth_dndtmean;
    lent = numel(t);
    
    indices = options.indices;
    weights = options.weights;
    sigma_s = options.sigma_s;
    loss = options.loss;
    scale = options.scale;
    collect = options.collect;
    
    switch collect
        case "mean"
            collect_fun = @(loss) mean(loss(:).*weights(:));
        case "sum"
            collect_fun = @(loss) sum(loss(:).*weights(:));
        case "logmean"
            collect_fun = @(loss) log(mean(loss(:).*weights(:)));
        case "logsum"
            collect_fun = @(loss) log(sum(loss(:).*weights(:)));
        otherwise
            error("collecting function is not given correctly")
    end

    lnqfqb = makevec(t, x(indices{1})); %log(qf/qb)
    lnqeqb = makevec(t, x(indices{2})); %log(qe/qb)
    alpha_prime = makevec(t, x(indices{3})); %alpha + 4
    kar1inf = makevec(t, x(indices{4})); % kappa_a*rho_{1,\infty}
    kass_d1 = makevec(t, x(indices{5})); % kappa_a*sigma_s/D_1
    supsat = makevec(t, x(indices{7})); % Supersaturation ratio, \rho_1(t)/rho_{1,\infty}
    shape_index = makevec(t, x(indices{8})); % shape index, \delta
    qfqb = exp(lnqfqb);
    qeqb = exp(lnqeqb);
    shape_factor = si2sf(shape_index);

    % ================================================== Part 1. parsing input arguments ==================================================
    
    % ================================================== Part 2. Calcaulate \partial_t(ln(\rho_c(t))) ==================================================
    [Nf_vec, Ne_vec] = Calculate_Nf_Ne(numberdata, shape_index');
    [Nf1_vec, Ne1_vec] = Calculate_Nf_Ne(numberdata - 1, shape_index');
    deltaNf_vec = Nf_vec - Nf1_vec;
    deltaNe_vec = Ne_vec - Ne1_vec;
    vals_vec = (1 + 1 ./ (numberdata' - 1)) .^ alpha_prime .* (qfqb' .^ deltaNf_vec)' .* (qeqb' .^ deltaNe_vec)';
    kna_vec = 4 * pi * sigma_s^2 * kar1inf .* shape_factor .* numberdata' .^ (2/3) ./ (1 + kass_d1 .* numberdata' .^ (1/3));   
    knd_vec = 4 * pi * sigma_s^2 * kar1inf .* shape_factor .* (numberdata' - 1) .^ (2/3) ./ (1 + kass_d1 .* (numberdata' - 1) .^ (1/3)) ./ vals_vec;
    dtlnrhoct = (mean(kna_vec,1,'omitnan') .* supsat - mean(knd_vec,1,'omitnan') - dndt') ./ (meannumberdata' - 1); % \partial_t(ln(\rho_c(t))) Eq. (S5-5)
    % ================================================== Part 2. Calcaulate \partial_t(ln(\rho_c(t))) ==================================================

    % ================================================== Part 3. Calcaulate \partial_t(\psi_n^p(t))) experiment&theory ==================================================
    costs = zeros(1, lent);
    prediction = cell(1, lent);
    for idx = 1:lent
        qfqb_ = qfqb(idx);
        qeqb_ = qeqb(idx);
        alp_ = alpha_prime(idx);
        kar1inf_ = kar1inf(idx);
        kass_d1_ = kass_d1(idx);
        supsat_ = supsat(idx);
        si_ = shape_index(idx);
        sf_ = shape_factor(idx);

        nval = expdata.nval_c{idx};
        pvalue_1 = expdata.pvalue_1_c{idx}; %p(n-1, t)
        pvalue = expdata.pvalue_c{idx}; %p(n, t)
        pcum_val = expdata.pcum_val_c{idx}; % \psi_n^{p}(t) = \Sigma_{j=n}^{\infty}p(j, t)
        ptimedev = expdata.ptimedev_c{idx}; % \partial_t\psi_n^p(t)

        [Nf, Ne] = Calculate_Nf_Ne(nval, si_);
        [Nf1, Ne1] = Calculate_Nf_Ne(nval - 1, si_);
        delta_Nf = Nf - Nf1;
        delta_Ne = Ne - Ne1;
        val = (1 + 1./(nval - 1)).^alp_.*(qfqb_.^delta_Nf).*(qeqb_.^delta_Ne);
        kna = 4*pi*sigma_s^2.*kar1inf_.*sf_.*(nval - 1).^(2/3)./(1 + kass_d1_ .* (nval - 1) .^ (1/3));
        knd = kna ./ val;
        pred = kna .* supsat_ .* pvalue_1 - knd .* pvalue - pcum_val .* dtlnrhoct(idx); % \partial_t\psi_{n}^p(t) Equation (S5-4)
        costs(idx) = lossReg(pred, ptimedev, loss, scale);
        prediction{idx} = [nval, pred];
    end
    cost = collect_fun(costs);
    % ================================================== Part 3. Calcaulate \partial_t(\psi_n^p(t))) experiment&theory ==================================================
end
% changed by soar8nalra@gmail.com
% Checked right