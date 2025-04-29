function [jn, meanjn] = Calculate_jn(x, expdata, options) % See Method Derivation of eq 2 for detail
    % ================================================== Part 1. parsing input arguments ==================================================
    T = expdata.timedata;
    
    indices = options.indices;
    
    lnqfqb = makevec(T, x(indices{1})); %ln(q_f/q_b)
    lnqeqb = makevec(T, x(indices{2})); %ln(q_e/q_b)
    alpha_prime = makevec(T, x(indices{3})); %\alpha + 4
    kar1inf = makevec(T, x(indices{4})); % \kappa_a * rho_{1, \infty}
    kass_d1 = makevec(T, x(indices{5})); % \kappa_a*\sigma_s/D_1
    sigma_s = options.sigma_s;
    supsat = makevec(T, x(indices{7})); %supersaturation ratio rho_1(t)/rho_{1,\infty}
    shape_index = makevec(T, x(indices{8})); %shape index(\delta)
    
    qfqb = exp(lnqfqb);
    qeqb = exp(lnqeqb);
    shape_factor = si2sf(shape_index);
    
    numberdata = expdata.smooth_numberdata;
    % ================================================== Part 1. parsing input arguments ==================================================
    
    % ================================================== Part 2. Calculate <j_n(t)> ==================================================
    lent = length(T);
    jn = cell(lent, 1);
    meanjn = cell(lent, 1);
    for idx = 1:lent
        qfqb_ = qfqb(idx);
        qeqb_ = qeqb(idx);
        alpha_ = alpha_prime(idx);
        kar1inf_ = kar1inf(idx);
        kass_d1_ = kass_d1(idx);
        supsat_ = supsat(idx); 
        si_ = shape_index(idx); 
        sf_ = shape_factor(idx);

        nval = min(numberdata(idx,:)):1:max(numberdata(idx,:));
        [Nf, Ne] = Calculate_Nf_Ne(nval, si_);
        [Nf1, Ne1] = Calculate_Nf_Ne(nval - 1, si_);
        delta_Nf = Nf - Nf1;
        delta_Ne = Ne - Ne1;
        
        % rho_{1,\infty}*q_{n}star/(q_{1}star q_{n-1}star)
        val = (1 + 1./(nval - 1)).^alpha_.*(qfqb_.^delta_Nf).*(qeqb_.^delta_Ne); 
        
        % rho_{1,\infty}*k_n^a Equation (M14)
        kna  = 4.*pi*sigma_s^2.*kar1inf_.*sf_.*nval.^(2/3)./(1+kass_d1_.*nval.^(1/3)); 
        
        % k_n^d by Detailed balance condition: k_{n-1}a/k_nd = q(n)star/(q(n-1)star*q1star)
        knd  = 4.*pi*sigma_s^2.*kar1inf_.*sf_.*(nval - 1).^(2/3)./(1+kass_d1_.*(nval - 1).^(1/3))./val; 
        
        % <j_n(t)> not <J_n(t)> equation (S5-6)
        jn_ = kna*supsat_ - knd;
        jntilde_ = jn_./(4*pi*sigma_s^2*sf_.*nval.^(2/3));
        jn{idx} = [nval', jn_', jntilde_'];
    end
    % ================================================== Part 2. Calculate <j_n(t)> ==================================================
    
end
% Checked right