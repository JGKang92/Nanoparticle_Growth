function [cost, costs, dtlnrhoct, prediction] = Calculate_Cost_Ex_situ_fast2(x, expdata, options)
    % ================================================== Part 1. parsing input arguments ==================================================
    % See the Eqs. (S5-3), (S5-4), and (S5-5) for details
    t = expdata.tval(:)';
    meannumberdata = expdata.smooth_meannumberdata;
    numberdata = expdata.numberdata;
    dndt = expdata.smooth_dndtmean;
    % lent = numel(t);
    
    indices = options.indices;
    weights = options.weights;
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
    sigma_s = unique(x(indices{6}));
    supsat = makevec(t, x(indices{7})); % Supersaturation ratio, \rho_1(t)/rho_{1,\infty}
    shape_index = makevec(t, x(indices{8})); % shape index, \delta

    qfqb = exp(lnqfqb);
    qeqb = exp(lnqeqb);
    shape_factor = si2sf(shape_index);
    lent_original = length(expdata.timedata);
    tidxList = nan(1, lent_original);
    for idx = 1:lent_original
        tidxList(idx) = find(t == expdata.timedata(idx));
    end
    qfqb_ = qfqb(tidxList);
    qeqb_ = qeqb(tidxList);
    alp_ = alpha_prime(tidxList);
    kar1inf_ = kar1inf(tidxList);
    kass_d1_ = kass_d1(tidxList);
    si_ = shape_index(tidxList);
    sf_ = shape_factor(tidxList);
    supsat_ = supsat(tidxList);
    % ================================================== Part 1. parsing input arguments ==================================================

    % ================================================== Part 2. Calcaulate \partial_t(ln(\rho_c(t))) ==================================================
    [Nf_vec, Ne_vec] = Calculate_Nf_Ne(numberdata, si_');
    [Nf1_vec, Ne1_vec] = Calculate_Nf_Ne(numberdata - 1, si_');
    deltaNf_vec = Nf_vec - Nf1_vec;
    deltaNe_vec = Ne_vec - Ne1_vec;
    vals_vec = (1 + 1 ./ (numberdata' - 1)) .^ alp_ .* (qfqb_' .^ deltaNf_vec)' .* (qeqb_' .^ deltaNe_vec)';
    kna_vec = 4 * pi * sigma_s^2 * kar1inf_ .* sf_ .* numberdata' .^ (2/3) ./ (1 + kass_d1_ .* numberdata' .^ (1/3));   
    knd_vec = 4 * pi * sigma_s^2 * kar1inf_ .* sf_ .* (numberdata' - 1) .^ (2/3) ./ (1 + kass_d1_ .* (numberdata' - 1) .^ (1/3)) ./ vals_vec;
     % \partial_t(ln(\rho_c(t))) Eq. (S5-5)
    dtlnrhoct = (mean(kna_vec,1,'omitnan') .* supsat_ - mean(knd_vec,1,'omitnan') - dndt(tidxList)') ./ (meannumberdata(tidxList)' - 1);
    if ~isvector(dtlnrhoct) || numel(dtlnrhoct) ~= lent_original
        error("dtlnrhoct should be vector of length %d", lent_original)
    end
    % ================================================== Part 2. Calcaulate \partial_t(ln(\rho_c(t))) ==================================================

    % ================================================== Part 3. Calcaulate \partial_t(\psi_n^p(t))) experiment&theory ==================================================
    costs = zeros(1, lent_original);
    prediction = cell(1, lent_original);
    for idx = 1:lent_original
        tidx = tidxList(idx);
        qfqb_ = qfqb(tidx);
        qeqb_ = qeqb(tidx);
        alp_ = alpha_prime(tidx);
        kar1inf_ = kar1inf(tidx);
        kass_d1_ = kass_d1(tidx);
        supsat_ = supsat(tidx);
        si_ = shape_index(tidx);
        sf_ = shape_factor(tidx);

        nval = expdata.nval_c{tidx}; % (lent, 1) cell
        pvalue_1 = expdata.pvalue_1_c{tidx}; %p(n-1, t)
        pvalue = expdata.pvalue_c{tidx}; %p(n, t)
        pcum_val = expdata.pcum_val_c{tidx}; % \psi_n^{p}(t) = \Sigma_{j=n}^{\infty}p(j, t)
        ptimedev = expdata.ptimedev_c{tidx}; % \partial_t\psi_n^p(t)
        
        TF = nval > expdata.minnum(idx) & nval < expdata.maxnum(idx);
        nval = nval(TF);
        pvalue_1 = pvalue_1(TF);
        pvalue = pvalue(TF);
        pcum_val = pcum_val(TF);
        ptimedev = ptimedev(TF);
        
        [Nf, Ne] = Calculate_Nf_Ne(nval, si_);
        [Nf1, Ne1] = Calculate_Nf_Ne(nval - 1, si_);
        delta_Nf = Nf - Nf1;
        delta_Ne = Ne - Ne1;
        val = (1 + 1./(nval - 1)).^alp_.*(qfqb_.^delta_Nf).*(qeqb_.^delta_Ne);

        %k_{n-1}^a*\rho_{1,\infty}
        kna = 4*pi*sigma_s^2.*kar1inf_.*sf_.*(nval - 1).^(2/3)./(1 + kass_d1_ .* (nval - 1) .^ (1/3));
        %k_{n}^d
        knd = kna ./ val;
        pred = kna .* supsat_ .* pvalue_1 - knd .* pvalue - pcum_val .* dtlnrhoct(idx); % \partial_t\psi_{n}^p(t) Equation (S5-4)
        costs(idx) = lossReg(pred, ptimedev, loss, scale);
        prediction{idx} = [nval(:), pred(:), ptimedev(:)];
    end
    cost = collect_fun(costs);
    % ================================================== Part 3. Calcaulate \partial_t(\psi_n^p(t))) experiment&theory ==================================================
end
% changed by soar8nalra@gmail.com