function [exp_mundata, muns_mn, muns_log] = Extract_mun_In_situ_For_Figure5(xfit, expdata, options)
    % ================================================== Part 1. parsing input arguments ==================================================
    smooth_numberdata = expdata.smooth_numberdata;
    smooth_meannumberdata = expdata.smooth_meannumberdata;
    if isfield(expdata, "tval")
        T = expdata.tval;
    else
        T = expdata.timedata;
    end
    
    indices = options.indices;
    
    lent = length(T);
    
    lnqfqb = makevec(T, xfit(indices{1})); %ln(q_e/q_b)
    lnqeqb = makevec(T, xfit(indices{2})); %ln(q_f/q_b)
    alpha_prime = makevec(T, xfit(indices{3})); %\alpha + 4
    kar1inf = makevec(T, xfit(indices{4})); % \kappa_a * rho_{1, \infty}
    kass_d1 = makevec(T, xfit(indices{5})); % \kappa_a*\sigma_s/D_1
    if isfield(options,"sigma_s")
        sigma_s = options.sigma_s;
    else
        sigma_s = xfit(indices{6});
    end
    supsat = makevec(T, xfit(indices{7})); %supersaturation ratio rho_1(t)/rho_{1,\infty}
    shape_index = makevec(T, xfit(indices{8})); %shape index(\delta)
    
    qfqb = exp(lnqfqb);  
    qeqb = exp(lnqeqb);  
    shape_factor = si2sf(shape_index);
    
    if isfield(expdata, 'smooth_minnum')
        minval = min(expdata.smooth_minnum,[],'all','omitnan');
        maxval = max(expdata.smooth_maxnum,[],'all','omitnan');
    else
        minval = min(smooth_numberdata,[],'all','omitnan');
        maxval = max(smooth_numberdata,[],'all','omitnan');
    end
    % ================================================== Part 1. parsing input arguments ==================================================
    % exp_mudata_table = nan(maxval-minval+1, lent);
    exp_mudata_table1 = nan(maxval-minval+1, lent);
    correction = zeros(lent,1);
    meanf = csapi(T, smooth_meannumberdata);
    dmeanf = fnder(meanf,1);
    % Calculate correction arrry
    for tidx = 1:lent
        qfqb_ = qfqb(tidx);
        qeqb_ = qeqb(tidx);
        alp_ = alpha_prime(tidx);
        kar1inf_ = kar1inf(tidx);
        kass_d1_ = kass_d1(tidx);
        supsat_ = supsat(tidx); 
        si_ = shape_index(tidx); 
        sf_ = shape_factor(tidx);

        nval = smooth_numberdata(tidx,:);
        nval = nval(~isnan(nval));
        [Nf, Ne] = Calculate_Nf_Ne(nval, si_);
        [Nf1,Ne1] = Calculate_Nf_Ne(nval-1, si_);
        delta_Nf = Nf - Nf1;
        delta_Ne = Ne - Ne1;
        [kna, knd] = Calculate_kna_knd(nval, qfqb_, qeqb_, alp_, kar1inf_, kass_d1_, sigma_s, sf_, delta_Nf, delta_Ne);
        correction(tidx) = (mean(kna)*supsat_ - mean(knd) - fnval(dmeanf,T(tidx))) ...
                ./(smooth_meannumberdata(tidx)-1); % Equation (S5-5)
    end

    for idx = 1:lent
        qfqb_ = qfqb(idx);
        qeqb_ = qeqb(idx);
        alp_ = alpha_prime(idx);
        kar1inf_ = kar1inf(idx);
        kass_d1_ = kass_d1(idx);
        supsat_ = supsat(idx); 
        si_ = shape_index(idx); 
        sf_ = shape_factor(idx);
        
        nval = expdata.p_cum_dev{idx}(:,1);
        ptimedev = expdata.p_cum_dev{idx}(:,2);
        idx1 = find(expdata.pcumval{idx}(:,1) == nval(1));
        idxList = idx1-1:idx1 + length(nval);
        pvalue = expdata.pcumval{idx}(idxList(2:end),2);
        if idx1 == 1
            pvalue_1 = [0;expdata.pcumval{idx}(idxList(2:end-1),2)];
        else
            pvalue_1 = expdata.pcumval{idx}(idxList(1:end-1),2);
        end
        pcum_val = pvalue(1:end-1);
        pvalue = -diff(pvalue); %p(n, t)
        pvalue_1 = -diff(pvalue_1);

        % nval = expdata.nval_c{idx};
        % pvalue = expdata.pvalue_c{idx};
        % pvalue_1 = expdata.pvalue_1_c{idx};
        ratiodata = pvalue./pvalue_1;
        % pcum_val = expdata.pcum_val_c{idx};
        % p_time_dev_val = expdata.ptimedev_c{idx};
        
        % ratiodata = ones(size(pvalue));
        [Nf, Ne] = Calculate_Nf_Ne(nval, si_);
        [Nf1,Ne1] = Calculate_Nf_Ne(nval-1, si_);
        delta_Nf = Nf - Nf1;
        delta_Ne = Ne - Ne1;
        val = (1 + 1./(nval-1)).^alp_.*(qfqb_.^delta_Nf).*(qeqb_.^delta_Ne);  
        knd = 4.*pi.*sigma_s^2.*kar1inf_.*sf_.*(nval-1).^(2/3)./(1+kass_d1_.*(nval-1).^(1/3))./val; 
        exp_mudata = -0.5.*log((((ptimedev + correction(idx).*pcum_val)./(knd.*pvalue) + 1).*ratiodata).^2) ...
                + log(supsat_); 
        % exp_mudata_table(nval-minval+1,idx) = movmean(exp_mudata,1); 
        % exp_mudata(exp_mudata >= inf | exp_mudata <= -inf) = nan; 
        exp_mudata_table1(nval-minval+1,idx) = movmean(exp_mudata,100); 
    end
    mean_mudata = mean(exp_mudata_table1,2,'omitnan'); 
    std_mudata = std(exp_mudata_table1,1,2,'omitnan'); 
    nanidx = ~isnan(mean_mudata); 
    xval = minval:maxval;
    xval = xval(nanidx);
    mean_mudata = mean_mudata(nanidx); 
    std_mudata = std_mudata(nanidx); 
    std_idx = find(std_mudata >= 0); 
    mean_mudata = mean_mudata(std_idx); 
    std_mudata = std_mudata(std_idx); 
    xval = xval(std_idx);
    exp_mundata = [xval', mean_mudata, std_mudata];
    
    % Chemical potential from kna and knd assuming that shape index = mean of shape index
    nval = round(logspace(1,4,20));
    shape_index = mean(xfit(indices{8}));
    shapefactor = si2sf(shape_index);
    [Nf, Ne] = Calculate_Nf_Ne(nval, shape_index);
    [Nf1,Ne1] = Calculate_Nf_Ne(nval-1, shape_index);
    delta_Nf = Nf - Nf1;
    delta_Ne = Ne - Ne1;
    [kna, knd] = Calculate_kna_knd(nval, mean(qfqb), mean(qeqb), mean(alpha_prime), mean(kar1inf), mean(kass_d1), sigma_s, shapefactor, delta_Nf, delta_Ne);
    muns_log = [nval', log(knd./kna)'];

    nval = 20:1e4;
    [Nf, Ne] = Calculate_Nf_Ne(nval, shape_index);
    [Nf1,Ne1] = Calculate_Nf_Ne(nval-1, shape_index);
    delta_Nf = Nf - Nf1;
    delta_Ne = Ne - Ne1;
    muns_mn = [nval; -delta_Nf.*log(mean(qfqb)) - delta_Ne.*log(mean(qeqb)) - mean(alpha_prime).*log(1+1./nval)]';
    % if unique(shape_index) == -1
    %     muns_mn = [nval; -(2/3).*log(mean(qfqb))./nval.^(1/3) - (1/3).*log(mean(qeqb))./nval.^(2/3) - mean(alpha_prime)./nval]';
    % end
end