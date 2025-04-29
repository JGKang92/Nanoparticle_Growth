function supsats = guess_supsat_Ex_situ2(x0, expdata, options)
    numInitialPoints = size(x0,1);
    if numInitialPoints == 1
        supsats = guess_supsat1(x0, expdata, options);
    else
        supsats = cell(numInitialPoints, 1);
        parfor ii = 1:numInitialPoints
            supsats{ii} = guess_supsat1(x0(ii,:), expdata, options);
        end
        supsats = cell2mat(supsats);
    end
end
function supsat = guess_supsat1(x0, expdata, options)
% This function guess supersaturation ratio using Eq. (M10b)
    indices = options.indices;
    T = expdata.tval;
    sigma_s = x0(options.indices{6});
    expdata = r2n_optim(expdata, sigma_s, options.IP_method, options.smooth_window);
    meannumberdata = expdata.smooth_meannumberdata;
    varnumberdata = expdata.smooth_varnumberdata;
    
    meannumber_pp = csapi(T, meannumberdata);
    f_dmeannumber = fnder(meannumber_pp, 1);
    dmeannumber = fnval(f_dmeannumber, T);
    varnumber_pp = csapi(T, varnumberdata);
    f_varnumber = fnder(varnumber_pp, 1);
    dvarnumber = fnval(f_varnumber, T);

    lnqfqb = makevec(T,x0(indices{1}));
    lnqeqb = makevec(T,x0(indices{2}));
    alpha_prime = makevec(T,x0(indices{3}));
    kar1inf = makevec(T,x0(indices{4}));
    kass_d1 = makevec(T,x0(indices{5}));
    shape_index = makevec(T,x0(indices{8}));
    
    qfqb = exp(lnqfqb);
    qeqb = exp(lnqeqb);
    shape_factor = si2sf(shape_index);
    
    lent = numel(T);
    deln_delkna_mean = nan(1, lent);
    deln_delknd_mean = nan(1, lent);

    for ii = 1:lent
        qfqb_ = qfqb(ii);
        qeqb_ = qeqb(ii);
        alp_ = alpha_prime(ii);
        kar1inf_ = kar1inf(ii);
        kass_d1_ = kass_d1(ii);
        si_ = shape_index(ii);
        sf_ = shape_factor(ii);
        nval = expdata.nval;
        nidx = ~isnan(nval);
        nval = nval(nidx);
        [Nf, Ne] = Calculate_Nf_Ne(nval, si_);
        [Nf1, Ne1] = Calculate_Nf_Ne(nval-1, si_);
        delta_Nf = Nf - Nf1;
        delta_Ne = Ne - Ne1;
        % val = (1 + 1 ./ (nval - 1)) .^ alp_ .* (qfqb_ .^ delta_Nf) .* (qeqb_ .^ delta_Ne);
        % kna = 4 * pi * sigma_s^2 * kar1inf_ .* sf_ .* nval .^ (2/3) ./ (1 + kass_d1_ .* nval .^ (1/3));
        % knd = 4 * pi * sigma_s^2 * kar1inf_ .* sf_ .* (nval - 1) .^ (2/3) ./ (1 + kass_d1_ .* (nval - 1) .^ (1/3)) ./ val;
        [kna, knd] = Calculate_kna_knd(nval, qfqb_, qeqb_, alp_, kar1inf_, kass_d1_, sigma_s, sf_, delta_Nf, delta_Ne);
        pn = expdata.pnt_interp(ii,:);
        deln = nval - expdata.smooth_meannumberdata(ii);
        % deln = nval - sum(nval.*pn);
        delkna = kna - sum(kna.*pn, 2);
        delknd = knd - sum(knd.*pn, 2);
        deln_delkna_mean(ii) = sum(deln.*delkna.*pn);
        deln_delknd_mean(ii) = sum(deln.*delknd.*pn);
    end
    % Eq. (M10b)
    supsat = ((dvarnumber(:) - dmeannumber(:))./2 + deln_delknd_mean(:))./deln_delkna_mean(:);
    
    supsat = supsat - min(supsat) + 1;
    % supsat(supsat<0) = 0;
    supsat = supsat(:)';
end
% Made by soar8nalra@gmail.com