function expdata_new = r2n_optim(expdata, sigma_s, IP_method, smooth_window)
    expdata_new = expdata;
    t = expdata.timedata;
    T = expdata.tval;
    lent = numel(T);
    r3data = expdata.r3data;
    numberdata = round(r3data./sigma_s.^3);
    maxnum = max(numberdata,[],2,'omitnan');
    minnum = min(numberdata,[],2,'omitnan');
    % Obtain most probable value and histogram data
    lent_original = size(r3data,1);
    
    % Calculate 1 - cdf
    nval = min(minnum):max(maxnum);
    len_n = length(nval);
    pntdata = nan(lent_original, len_n);
    survivaldata = nan(lent_original, len_n);
    for idx=1:lent_original
        num_data = numberdata(idx,:);
        num_data = num_data(~isnan(num_data));
        [value, edges] = histcounts(num_data, ...
            Normalization = 'cdf', ...
            BinLimits = [nval(1), nval(end)]);
        value_ = [0, value];
        switch IP_method
            case "pchip"
                cfun = @(x)pchip(edges, value_, x); 
            case "makima"
                cfun = @(x)makima(edges, value_, x);
            case "csapi"
                cfun = @(x)csapi(edges, value_, x);
        end
        
        cumval = cfun(nval);
        dcumval = gradient(cumval);
        dcumval(dcumval<0) = 0;
        cumval_new = cumval(1) + cumsum(dcumval);
        cumval_new(cumval_new<0) = 0;
        cumval_new(cumval_new>1) = 1;
        survivaldata(idx,:) = 1-cumval_new;
        
        pval = gradient(cumval_new);
        pntdata(idx,:) = pval./sum(pval);
    end
    % meannumberdata = sum(nval.*pntdata,2);
    % varnumberdata = sum(nval.^2.*pntdata,2) - meannumberdata.^2;
    
    % Time point interpolation of survival prob: survival_interp
    % and obtain time derivative of survival prob: pcum_dev_time
    survival_interp = nan(lent, len_n);
    pcum_dev_time = nan(lent, len_n);
    for idx = 1:len_n
        p = survivaldata(:,idx);
        switch IP_method
            case "pchip"
                pcum_time_f = pchip(t, p);
            case "makima"
                pcum_time_f = makima(t, p);
            case "csapi"
                pcum_time_f = csapi(t, p);
        end
        p = ppval(pcum_time_f,T);
        p(p<0) = 0;
        p(p>1) = 1;
        survival_interp(:,idx) = p;
        d_pcum_time = fnder(pcum_time_f, 1);
        pcum_dev_time(:,idx) = fnval(d_pcum_time, T);
    end
    % find and replace point where survival probability increases
    for ii = 1:lent
        p = survival_interp(:,ii);
        grad = gradient(p);
        TF = grad>0;
        grad(TF) = 0;
        survival_interp(:,ii) = p(1) + cumsum(grad);
    end
    % Time point interpolation of pnt
    pnt_interp = nan(lent, len_n);
    for idx = 1:len_n
        p = pntdata(:, idx);
        switch IP_method
            case "pchip"
                pcum_time_f = pchip(t, p);
            case "makima"
                pcum_time_f = makima(t, p);
            case "csapi"
                pcum_time_f = csapi(t, p);
        end
        p = ppval(pcum_time_f, T);
        p(p<0) = 0;
        pnt_interp(:,idx) = p;
    end
    pnt_interp = pnt_interp./sum(pnt_interp,2);
    smooth_meannumberdata = sum(nval.*pnt_interp,2);
    smooth_varnumberdata = sum(nval.^2.*pnt_interp, 2) - smooth_meannumberdata.^2;
    smooth_relstd = sqrt(smooth_varnumberdata)./smooth_meannumberdata;
    if ~isempty(smooth_window)
        smooth_meannumberdata = smoothdata(smooth_meannumberdata, 'gaussian',smooth_window(1));
        smooth_relstd = smoothdata(smooth_relstd, 'gaussian', smooth_window(2));
        smooth_varnumberdata = (smooth_relstd.*smooth_meannumberdata).^2;
    end
    smooth_smdata = smooth_varnumberdata + smooth_meannumberdata.^2;
    p_mean = csapi(T, smooth_meannumberdata);
    f = fnder(p_mean, 1);
    smooth_dndtmean = ppval(f, T);

    ptimedev_c = cell(lent, 1);
    pvalue_c = cell(lent, 1);
    pvalue_1_c = cell(lent, 1);
    nval_c = cell(lent, 1);
    pcum_val_c = cell(lent, 1);
    for ii = 1:lent
        nval_ = nval(1:end-1);
        ptimedev = pcum_dev_time(ii,1:end-1); %\partial_t\psi(n, t)
        pvalue = pnt_interp(ii,2:end); %p(n,t)
        pvalue_1 = pnt_interp(ii,1:end-1); %p(n-1,t)
        pcum_val = survival_interp(ii,1:end-1); %\psi(n, t)
        
        ptimedev_c{ii} = ptimedev;
        nval_c{ii} = nval_;
        pvalue_c{ii} = pvalue;
        pvalue_1_c{ii} = pvalue_1;
        pcum_val_c{ii} = pcum_val;
    end
    
    expdata_new.pnt_interp = pnt_interp;
    expdata_new.nval = nval;
    expdata_new.pntdata = pntdata;
    expdata_new.minnum = minnum; % lent*1 vector
    expdata_new.maxnum = maxnum; % lent*1 vector
    expdata_new.numberdata = numberdata; 
    expdata_new.nval_c = nval_c; %lent*1 cell
    expdata_new.ptimedev_c = ptimedev_c;%lent*1 cell
    expdata_new.pvalue_c = pvalue_c;%lent*1 cell
    expdata_new.pvalue_1_c = pvalue_1_c;%lent*1 cell
    expdata_new.pcum_val_c = pcum_val_c;%lent*1 cell
    
    expdata_new.p_mean = p_mean;
    expdata_new.smooth_meannumberdata = smooth_meannumberdata;
    expdata_new.smooth_smdata = smooth_smdata;
    expdata_new.smooth_dndtmean = smooth_dndtmean;
    expdata_new.smooth_varnumberdata = smooth_varnumberdata;

end