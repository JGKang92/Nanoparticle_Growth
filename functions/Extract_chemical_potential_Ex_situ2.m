function [exp_mundata, muns_mn, muns_log] = Extract_chemical_potential_Ex_situ2(xfit, expdata, options)
    % ================================================== Part 1. parsing input arguments ==================================================
    % smooth_numberdata = expdata.smooth_numberdata;
    smooth_meannumberdata = expdata.smooth_meannumberdata;
    T = expdata.tval;
    t = expdata.timedata;
    
    indices = options.indices;
    
    lent = length(t);
    
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
    
    minval = min(expdata.minnum,[],'all','omitnan');
    maxval = max(expdata.maxnum,[],'all','omitnan');
    % ================================================== Part 1. parsing input arguments ==================================================
    
    meanf = csapi(T, smooth_meannumberdata);
    dmeanf = fnder(meanf,1);
    
    % Calculate dtlnrhoct arrry
    dtlnrhoct = zeros(lent,1);
    for idx = 1:lent
        tidx = find(t(idx) == T);
        qfqb_ = qfqb(tidx);
        qeqb_ = qeqb(tidx);
        alp_ = alpha_prime(tidx);
        kar1inf_ = kar1inf(tidx);
        kass_d1_ = kass_d1(tidx);
        supsat_ = supsat(tidx); 
        si_ = shape_index(tidx); 
        sf_ = shape_factor(tidx);

        % nval = smooth_numberdata(idx,:);
        nval = expdata.nval;
        pn = expdata.pntdata(idx,:);
        nval = nval(~isnan(nval));
        [Nf, Ne] = Calculate_Nf_Ne(nval, si_);
        [Nf1,Ne1] = Calculate_Nf_Ne(nval-1, si_);
        delta_Nf = Nf - Nf1;
        delta_Ne = Ne - Ne1;
        [kna, knd] = Calculate_kna_knd(nval, qfqb_, qeqb_, alp_, kar1inf_, kass_d1_, sigma_s, sf_, delta_Nf, delta_Ne);
        kna_mean = sum(kna.*pn);
        knd_mean = sum(knd.*pn);
        dtlnrhoct(idx) = (kna_mean*supsat_ - knd_mean - fnval(dmeanf,T(tidx))) ...
                ./(smooth_meannumberdata(tidx)-1); % Equation (S5-5)
    end
    n = numel(expdata.nval_c{1});
    exp_mudata_table = nan(n, lent);
    exp_mudata_table1 = nan(n, lent);
    for idx = 1:lent
        tidx = find(T == t(idx));
        qfqb_ = qfqb(tidx);
        qeqb_ = qeqb(tidx);
        alp_ = alpha_prime(tidx);
        kar1inf_ = kar1inf(tidx);
        kass_d1_ = kass_d1(tidx);
        supsat_ = supsat(tidx); 
        si_ = shape_index(tidx); 
        sf_ = shape_factor(tidx);

        nval = expdata.nval_c{tidx};
        p_time_dev_val = expdata.ptimedev_c{tidx};
        pcum_val = expdata.pcum_val_c{tidx};
        pvalue = expdata.pvalue_c{tidx};
        pvalue_1 = expdata.pvalue_1_c{tidx};
        ratiodata = pvalue./pvalue_1;
        % ratiodata = ones(size(pvalue));
        [Nf, Ne] = Calculate_Nf_Ne(nval, si_);
        [Nf1,Ne1] = Calculate_Nf_Ne(nval-1, si_);
        delta_Nf = Nf - Nf1;
        delta_Ne = Ne - Ne1;
        val = (1 + 1./(nval-1)).^alp_.*(qfqb_.^delta_Nf).*(qeqb_.^delta_Ne);  
        knd = 4.*pi.*sigma_s^2.*kar1inf_.*sf_.*(nval-1).^(2/3)./(1+kass_d1_.*(nval-1).^(1/3))./val; 
        exp_mudata = -0.5.*log((((p_time_dev_val + dtlnrhoct(idx).*pcum_val)./(knd.*pvalue) + 1).*ratiodata).^2) ...
            + log(supsat_); 
        exp_mudata_table(:,idx) = movmean(exp_mudata,1); 
        exp_mudata(exp_mudata >= inf | exp_mudata <= -inf) = nan; 
        exp_mudata_table1(:,idx) = exp_mudata; 
        if size(exp_mudata_table1,1) ~= n || size(exp_mudata_table1,2) ~= lent
            error("The size of exp_mutdata_table1 should be (%d %d)",n, lent)
        end
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
    % npoints = 60;
    % xidx_1 = floor(logspace(log10(xval(1)),log10(xval(end)),npoints));
    % xidx = zeros(1,npoints);
    %for idx=1:npoints
    %    [~, xidx(idx)] = min(abs(xval - xidx_1(idx)));
    %end
    % exp_mundata = [xval(xidx)', mean_mudata(xidx), std_mudata(xidx)];
    exp_mundata = [xval', mean_mudata, std_mudata];
    % Chemical potential from kna and knd assuming that shape index = 1.5
    nval = round(logspace(1,4,20));
    shape_index = mean(xfit(indices{8}));
    shapefactor = si2sf(shape_index);
    [Nf, Ne] = Calculate_Nf_Ne(nval, shape_index);
    [Nf1,Ne1] = Calculate_Nf_Ne(nval-1, shape_index);
    delta_Nf = Nf - Nf1;
    delta_Ne = Ne - Ne1;
    [kna, knd] = Calculate_kna_knd(nval, mean(qfqb), mean(qeqb), mean(alpha_prime), mean(kar1inf), mean(kass_d1), sigma_s, shapefactor, delta_Nf, delta_Ne);
    % val = (1 + 1./(nval-1)).^mean(alpha_prime).*(mean(qfqb).^delta_Nf).*(mean(qeqb).^delta_Ne);   
    % kna = 4.*pi.*sigma_s^2.*mean(kar1inf).*shapefactor.*nval.^(2/3)./(1+mean(kass_d1).*nval.^(1/3));
    % knd = 4.*pi.*sigma_s^2.*mean(kar1inf).*shapefactor.*(nval-1).^(2/3)./(1+mean(kass_d1).*(nval).^(1/3))./val;
    muns_log = [nval', log(knd./kna)'];

    nval = 20:10000;
    [Nf, Ne] = Calculate_Nf_Ne(nval, shape_index);
    [Nf1,Ne1] = Calculate_Nf_Ne(nval-1, shape_index);
    delta_Nf = Nf - Nf1;
    delta_Ne = Ne - Ne1;
    muns_mn = [nval', (-delta_Nf.*log(mean(qfqb)) - delta_Ne.*log(mean(qeqb)) - mean(alpha_prime).*log(1+1./nval))'];
end