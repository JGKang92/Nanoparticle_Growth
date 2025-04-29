function mun_s = muns_function(nval, delta_ef, delta_ee, alp_prime, mn)
    [Nf, Ne] = Calculate_Nf_Ne(nval+1, mn);
    [Nf1,Ne1] = Calculate_Nf_Ne(nval, mn);
    delta_Nf = Nf - Nf1;
    delta_Ne = Ne - Ne1;
    mun_s = delta_Nf.*delta_ef + delta_Ne.*delta_ee - alp_prime.*log(1+1./nval);
end