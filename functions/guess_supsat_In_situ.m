function supsats = guess_supsat_In_situ(x0, expdata, options)
    numInitialPoints = size(x0,1);
    if numInitialPoints == 1
        supsats = guess_supsat1(x0, expdata, options);
    else
        supsats = cell(numInitialPoints, 1);
        for ii = 1:numInitialPoints
            supsats{ii} = guess_supsat1(x0(ii,:), expdata, options);
        end
        supsats = cell2mat(supsats);
    end
end
function supsat = guess_supsat1(x0, expdata, options)
% This function guess supersaturation ratio using Eq. (M10b)
    indices = options.indices;
    sigma_s = options.sigma_s;

    timedata = expdata.timedata;
    numberdata = expdata.smooth_numberdata;
    meannumberdata = expdata.smooth_meannumberdata;
    varnumberdata = expdata.smooth_varnumberdata;
    
    meannumber_pp = csapi(timedata, meannumberdata);
    f_dmeannumber = fnder(meannumber_pp, 1);
    dmeannumber = fnval(f_dmeannumber, timedata);
    varnumber_pp = csapi(timedata, varnumberdata);
    f_varnumber = fnder(varnumber_pp, 1);
    dvarnumber = fnval(f_varnumber, timedata);

    lnqfqb = makevec(timedata,x0(indices{1}));
    lnqeqb = makevec(timedata,x0(indices{2}));
    alpha_prime = makevec(timedata,x0(indices{3}));
    kar1inf = makevec(timedata,x0(indices{4}));
    kass_d1 = makevec(timedata,x0(indices{5}));
    shape_index = makevec(timedata,x0(indices{8}));
    
    qfqb = exp(lnqfqb);
    qeqb = exp(lnqeqb);
    shape_factor = si2sf(shape_index);
    
    lent = numel(expdata.timedata);
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
        nval = numberdata(ii,:);
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
        deln = nval - mean(nval);
        delkna = kna - mean(kna);
        delknd = knd - mean(knd);
        deln_delkna_mean(ii) = mean(deln.*delkna);
        deln_delknd_mean(ii) = mean(deln.*delknd);
    end
    % Eq. (M10b)
    supsat = ((dvarnumber(:) - dmeannumber(:))./2 + deln_delknd_mean(:))./deln_delkna_mean(:);
    
    supsat = supsat - min(supsat) + 1;
    % supsat(supsat<0) = 0;
    supsat = supsat(:)';
end
% Made by soar8nalra@gmail.com