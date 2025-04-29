function [mean_theory, var_theory, totmdata] = Calculate_Mean_Var_Ex_situ2(x, expdata, options)
    % ================================================== Part 1. parsing input arguments ==================================================
    T = expdata.tval;
    indices = options.indices;
    
    lnqfqb = makevec(T, x(indices{1})); %ln(q_f/q_b)
    lnqeqb = makevec(T, x(indices{2})); %ln(q_e/q_b)
    alpha_prime = makevec(T, x(indices{3})); %\alpha + 4
    kar1inf = makevec(T, x(indices{4})); % \kappa_a * rho_{1, \infty}
    kass_d1 = makevec(T, x(indices{5})); % \kappa_a*\sigma_s/D_1
    sigma_s = x(indices{6});
    supsat = makevec(T, x(indices{7})); %supersaturation ratio rho_1(t)/rho_{1,\infty}
    shape_index = makevec(T, x(indices{8})); %shape index(\delta)
    if ~isscalar(sigma_s)
        error("\sigma_s should be scalar\n")
    end
    qfqb = exp(lnqfqb);
    qeqb = exp(lnqeqb);
    shape_factor = si2sf(shape_index);
    
    meannumberdata = expdata.smooth_meannumberdata;
    smdata = expdata.smooth_smdata;
    % numberdata = expdata.numberdata;
    dndtmean = expdata.smooth_dndtmean;
    % ================================================== Part 1. parsing input arguments ==================================================
    
    % ================================================== Part 2. Calculate <j_n(t)> and <n(t)j_n(t)>/<n(t)> ==================================================
    IP_method = options.IP_method;
    if IP_method == "csapi"
        IP_method = "spline";
    end
    lent = length(T);
    jn = nan(lent, 1);
    ntjnt = nan(lent, 1);

    for idx = 1:lent
        qfqb_ = qfqb(idx);
        qeqb_ = qeqb(idx);
        alpha_ = alpha_prime(idx);
        kar1inf_ = kar1inf(idx);
        kass_d1_ = kass_d1(idx);
        supsat_ = supsat(idx);
        
        si_ = shape_index(idx);
        sf_ = shape_factor(idx);

        % nval = numberdata(idx, :); % the observed number of monomers in nanoparticles at time index idx

        nval = expdata.nval;
        pn = expdata.pnt_interp(idx,:);
        
        TF = ~isnan(nval);
        nval = nval(TF);
        pn = pn(TF);
        
        % TF = nval>expdata.minnum(idx) & nval<expdata.maxnum(idx);
        % nval = nval(TF);

        [Nf, Ne] = Calculate_Nf_Ne(nval, si_);
        [Nf1, Ne1] = Calculate_Nf_Ne(nval - 1, si_);
        delta_Nf = Nf - Nf1;
        delta_Ne = Ne - Ne1;

        % rho_{1,\infty}*q_{n}star/(q_{1}star q_{n-1}star)
        val = (1 + 1./(nval - 1)).^alpha_.*(qfqb_.^delta_Nf).*(qeqb_.^delta_Ne);

        % rho_{1,\infty}*k_n^a Equation (M14)
        kna  = 4.*pi*sigma_s^2.*kar1inf_.*sf_.*nval.^(2/3)./(1+kass_d1_.*nval.^(1/3));
        nkna = nval.*kna;

        % k_n^d by Detailed balance condition: k_{n-1}a/k_nd = q(n)star/(q(n-1)star*q1star)
        knd  = 4.*pi*sigma_s^2.*kar1inf_.*sf_.*(nval - 1).^(2/3)./(1+kass_d1_.*(nval - 1).^(1/3))./val;
        nknd = (nval-1).*knd;
        kna_mean = sum(kna.*pn);
        knd_mean = sum(knd.*pn);
        nkna_mean = sum(nkna.*pn);
        nknd_mean = sum(nknd.*pn);
        % kna_mean = mean(kna);
        % knd_mean = mean(knd);
        % nkna_mean = mean(nkna);
        % nknd_mean = mean(nknd);
        % <j_n(t)> not <J_n(t)> equation (S5-6)
        jn(idx) = kna_mean*supsat_ - knd_mean;
        
        %<n(t)j_n(t)> equation (S5-10) divided by <n(t)>
        ntjnt(idx) = (nkna_mean*supsat_ - nknd_mean);
    end
    at = ntjnt./meannumberdata;
    at_f = griddedInterpolant(T, at, IP_method);
    atf = @(t)at_f(t);
    % ================================================== Part 2. Calculate <j_n(t)> and <n(t)j_n(t)>/<n(t)> ==================================================
    
    % ================================================== Part 3. Calculate Mean and Variance ==================================================
    dlogmeandata = dndtmean ./ meannumberdata;

    % the equation (S4-13)
    dlntotdata = (jn - dlogmeandata)./(meannumberdata - 1); 
    dlntot_f = griddedInterpolant(T, dlntotdata, IP_method);
    
    % knetdata = dlntotdata .* (1 - meannumberdata) + jn;
    % knet_f = griddedInterpolant(T, knetdata, IP_method);
    
    
    % mean_theory = zeros(lent, 1);
    T0 = T(1);
    % rho_{1,T}(t)/rho_{1,T}(t0)
    totmdata = zeros(lent, 1);
    for idx = 1:lent
        % rho_{1,T}(t)/rho_{1,T}(t0) the equation (S4-13)
        totmdata(idx) = exp(integral(@(t)dlntot_f(t), T0, T(idx))); %(S4-13)
    end
    % for idx = 1:lent
    %     mean_theory(idx) = meannumberdata(1).*exp(integral(@(t)knet_f(t), T0, T(idx)));
    % end

    mean_theory = meannumberdata;
    % rho_{1,T}(t0)/rho_{1,T}(t)
    totmdata_1 = 1./totmdata;
    totmdata_1f = griddedInterpolant(T, totmdata_1, IP_method);
    ratio_theory = zeros(lent, 1);
    for idx = 1:lent
        tt = T(idx);
        intf = @(t) (totmdata_1f(tt)./totmdata_1f(t)).*atf(t);
        % <n^2(t)>/<n(t)>, equation (S4-12) for q = 2
        ratio_theory(idx) = 1 + (smdata(1)/meannumberdata(1)-1).*totmdata_1f(tt) ...
            + 2.*integral(intf, T0, tt);
    end
    sm_theory = ratio_theory .* mean_theory;
    var_theory = sm_theory - mean_theory.^2;
    % ================================================== Part 3. Calculate Mean and Variance ==================================================
end