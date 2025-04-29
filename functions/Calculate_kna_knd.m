function [kna, knd] = Calculate_kna_knd(nval, qfqb_, qeqb_, alpha_, kar1inf_, kass_d1_, sigma_s, sf_, delta_Nf, delta_Ne)
    val = (1 + 1 ./ (nval - 1)) .^ alpha_ .* (qfqb_ .^ delta_Nf) .* (qeqb_ .^ delta_Ne);
    factor = 4 .* pi .* sigma_s^2 .* kar1inf_ .* sf_;
    kna = factor.*nval .^ (2/3)./(1 + kass_d1_.*nval.^(1/3));
    knd = factor.*(nval - 1).^(2/3)./(1 + kass_d1_.*(nval - 1).^(1/3))./ val;
end