function [mun_s, monomer_chemical_potential, free_energy] = Calculate_chemical_potential(x, expdata, options)
    indices = options.indices;

    if isfield(expdata, 'tval')
        t = expdata.tval;
    else
        t = expdata.timedata;
    end

    lnqfqb = makevec(t, x(indices{1}));
    lnqeqb = makevec(t, x(indices{2}));
    alpha_prime = makevec(t, x(indices{3}));
    supsat = makevec(t, x(indices{7})); %supersaturation ratio rho_1(t)/rho_{1,\infty}
    shape_index = makevec(t, x(indices{8})); %shape index(\delta)
    qfqb = exp(lnqfqb);
    qeqb = exp(lnqeqb);

    lent = length(t);
    nval = (1:1e4)';
    mun_s = zeros(1e4, lent);
    free_energy = zeros(numel(nval), lent);
    monomer_chemical_potential = zeros(1, lent);

    for idx = 1:lent
        qfqb_ = qfqb(idx);
        qeqb_ = qeqb(idx);
        alp_ = alpha_prime(idx);
        supsat_ = supsat(idx);
        si_ = shape_index(idx);

        [Nf, Ne] = Calculate_Nf_Ne(nval + 1, si_);
        [Nf1, Ne1] = Calculate_Nf_Ne(nval, si_);
        delta_Nf = Nf - Nf1;
        delta_Ne = Ne - Ne1;
        mun_s(:, idx) = -delta_Nf .* log(qfqb_) - delta_Ne .* log(qeqb_) - alp_ .* log(1 + 1 ./ nval);
        free_energy(:, idx) = -Nf .* log(qfqb_) - Ne .* log(qeqb_) - alp_ .* log(nval) - nval .* log(supsat_);
        monomer_chemical_potential(idx) = log(supsat_);
    end
end
% changed by soar8nalra@gmail.com