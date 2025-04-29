function [raw_noise, eq_m12_relvar, eta2_inf] = Calculate_approximate_relvar(expdata, result, t0idx)
    % See Equation (M12) for detail
    mean_theory = result.mean_theory;
    var_theory = result.var_theory;
    if isfield(expdata, "varnumberdata")
        varnumberdata = expdata.varnumberdata;
        meannumberdata = expdata.meannumberdata;
        raw_noise = varnumberdata ./ meannumberdata.^2;
    else
        tidxList = find(ismember(expdata.tval, expdata.timedata));
        mean_theory = mean_theory(tidxList);
        var_theory = var_theory(tidxList);
        raw_noise = expdata.varr3data ./ expdata.meanr3data.^2;
    end
    options = optimoptions(@lsqnonlin, 'Display', 'none');
    
    f0value = var_theory(t0idx)./mean_theory(t0idx).^2 - 1./mean_theory(t0idx);
    % eta2_inf = 0.78;
    
    fitfun = @(eta2_inf) (1./(mean_theory(t0idx:end)) + f0value.*(mean_theory(t0idx)./mean_theory(t0idx:end)).^2 + ...
        eta2_inf.*(1 - (mean_theory(t0idx)./mean_theory(t0idx:end)).^2)) - raw_noise(t0idx:end);
    eta2_infinit = 0.7;
    eta2_inf = lsqnonlin(fitfun, eta2_infinit, [], [], options);
    eq_m12_relvar = 1./(mean_theory(t0idx:end)) + f0value.*(mean_theory(t0idx)./mean_theory(t0idx:end)).^2 + ...
        eta2_inf.*(1 - (mean_theory(t0idx)./mean_theory(t0idx:end)).^2);
end