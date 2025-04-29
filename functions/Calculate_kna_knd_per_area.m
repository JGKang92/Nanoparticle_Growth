function [kna, knd] = Calculate_kna_knd_per_area(x, nval, expdata, options)
    t = expdata.timedata;
    sigma_s = options.sigma_s;
    indices = options.indices;
    lnqfqb = makevec(t, x(indices{1})); %log(qf/qb)
    lnqeqb = makevec(t, x(indices{2})); %log(qe/qb)
    alpha_prime = makevec(t, x(indices{3})); %alpha + 4
    kar1inf = makevec(t, x(indices{4})); % kappa_a*rho_{1,\infty}
    kass_d1 = makevec(t, x(indices{5})); % kappa_a*sigma_s/D_1
    
    shape_index = 1.5;
    qfqb = exp(lnqfqb);
    qeqb = exp(lnqeqb);
    shape_factor = si2sf(shape_index);
    
    [Nf, Ne] = Calculate_Nf_Ne(nval, shape_index);
    [Nf1,Ne1] = Calculate_Nf_Ne(nval-1, shape_index);
    delta_Nf = Nf - Nf1;
    delta_Ne = Ne - Ne1;
    
    [kna, knd] = Calculate_kna_knd(nval, ...
        mean(qfqb), mean(qeqb), mean(alpha_prime), mean(kar1inf), mean(kass_d1), sigma_s, ...
        shape_factor, delta_Nf, delta_Ne);
    area = 4*pi*sigma_s^2*shape_factor.*nval.^(2/3);
    kna = kna./area;
    knd = knd./area;
end
