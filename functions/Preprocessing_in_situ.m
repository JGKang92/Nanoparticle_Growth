function expdata = Preprocessing_in_situ(timedata, radidata, rho_s, dt, smooth_windows, filename, options)
arguments
    timedata
    radidata
    rho_s
    dt
    smooth_windows
    filename
    options.trim = 200;
    options.smooth_trajectory = true;
    options.window = 1;
end
    sigmas = (3/(4*pi*rho_s))^(1/3);
    radidata(radidata == 0) = nan;
    numberdata = round((radidata./sigmas).^3);
    numberdata(numberdata == 0) = nan;
    numberdata = remove_trajectories(numberdata, nannum = 3);
    [lent, ntrj] = size(numberdata);
    
    meannumberdata = mean(numberdata, 2, 'omitnan');
    varnumberdata = var(numberdata, 1, 2, 'omitnan');
    maxnum = max(numberdata,[],2,'omitnan');
    minnum = min(numberdata,[],2,'omitnan');

    %% trajectory smoothing using FFT filter method
    if options.smooth_trajectory
        smooth_numberdata = nan(lent, ntrj);
        smooth_dndtdata = nan(lent, ntrj);
        warning('off', 'signal:internal:filteringfcns:SignalLengthForIIR')
        for idx=1:ntrj
            fprintf("%3d/%3d\n",idx,ntrj)
            nidx = ~isnan(numberdata(:,idx)); % find the time indices where idx-th trajectory of numberdata is not nan
            rnum = sum(nidx); % number of time points that are not nan
            ndata = numberdata(nidx,idx); 
            ftrdata = [ndata(end:-1:1); ndata; ndata(end:-1:1)];
            tdata = timetable(seconds((1:length(ftrdata))'), ftrdata);
            smooth_tdata = lowpass(tdata, 0.01);
            sm_data = smooth_tdata.ftrdata(rnum+1:2*rnum);
            smooth_numberdata(nidx, idx) = sm_data;
            smpp = csapi(timedata(nidx),sm_data);
            smfprime = fnder(smpp,1);
            smooth_dndtdata(nidx,idx) = fnval(smfprime, timedata(nidx));
            for ii = 1:8
                fprintf("\b")
            end
        end
        smooth_numberdata = round(smooth_numberdata);
    else
        smooth_numberdata = numberdata;
        smooth_dndtdata = nan(lent, ntrj);
        for ii = 1:ntrj
            nidx = ~isnan(numberdata(:,ii));
            smpp = csapi(timedata(nidx),smooth_numberdata(nidx,ii));
            smfprime = fnder(smpp,1);
            smooth_dndtdata(:,ii) = fnval(smfprime, timedata);
        end
    end
    smooth_maxnum = max(smooth_numberdata,[],2,'omitnan');
    smooth_minnum = min(smooth_numberdata,[],2,'omitnan');
    smooth_meannumberdata = mean(smooth_numberdata, 2, 'omitnan');
    smooth_varnumberdata = var(smooth_numberdata, 1, 2, 'omitnan');
    smooth_smdata = smooth_varnumberdata + smooth_meannumberdata.^2;
    smooth_dndtmean = mean(smooth_dndtdata, 2, 'omitnan');
    %% Obtain cdf P(n <= z, t) and pdf P(n = z, t)
    window = options.window;
    nbin = 30;
    pntdata = cell(lent, 1);
    edges = zeros(nbin + 1, lent);
    smooth_value = zeros(nbin, lent);
    for idx = 1:lent
        num_data = smooth_numberdata(idx,:);
        num_data = num_data(~isnan(num_data));
        % [value1, edge] = histcounts(num_data, nbin, 'Normalization', 'cdf'); % P(n <= z)
        % [value2,  ~  ] = histcounts(num_data, nbin, 'Normalization', 'pdf'); % P(n = z)
        % value(:,idx) = 1 - value1 + value2; % P(n >= z)
        [value, edges(:,idx)] = histcounts(num_data, nbin, 'Normalization', 'pdf');
        smooth_value(:,idx) = movmean(value,5); % moving averaged pdf values
        pntdata{idx} = [edges(1:end-1,idx) + (edges(2,idx) - edges(1,idx)) / 2, movmean(value',window)];
    end
    %% Obtain P(n >= z, t) for z = smooth_minnum:smooth_maxnum by interpolation
    % pcumval = cell(lent,1);
    % for idx=1:lent
    %     edge = edges(1:end-1,idx) + (edges(2,idx) - edges(1,idx))/2;
    %     pfun = @(x)pchip(edge,smooth_value(:,idx),x);
    %     nval = smooth_minnum(idx):smooth_maxnum(idx);
    %     pval = pfun(nval);
    %     pcumval{idx} = [nval', pval'];
    % end
    %% Cumulative function
    pcumval = cell(lent,1);
    for idx=1:lent
        pfun = @(x)pchip(edges(1:end-1,idx),smooth_value(:,idx),x); %interpolation of moving averaged pdf values
        nval = smooth_minnum(idx):smooth_maxnum(idx);
        norm_const = sum(pfun(nval)); % normalization constant
        pval = pfun(nval) ./ norm_const; % normalize
        pval = pval(end:-1:1);
        pcumval_temp = cumsum(pval);
        pcumval{idx} = [nval', pcumval_temp(end:-1:1)'];
    end
    smooth_globalminnum = min(smooth_minnum, [], 'omitnan');
    smooth_globalmaxnum = max(smooth_maxnum, [], 'omitnan');
    pcum_time = nan(length(smooth_globalminnum:smooth_globalmaxnum),lent);
    idx_adj = smooth_globalminnum - 1;
    for idx=1:lent
        nval = smooth_minnum(idx):smooth_maxnum(idx);
        pcum_time(nval - idx_adj,idx) = pcumval{idx}(:,2);
    end
    %% Obtain the time derivative of pcum_time
    p_cum_dev = cell(lent,1);
    % derivative at first time index
    nval1 = smooth_minnum(1):smooth_maxnum(1);
    nval2 = smooth_minnum(2):smooth_maxnum(2);
    nval_common = intersect(nval1,nval2) - idx_adj;
    nval_common = nval_common(1:end-1);
    p_cum_dev{1} = [nval_common' + idx_adj, (pcum_time(nval_common,2) - pcum_time(nval_common,1)) / dt];
    % derivative at last time index
    nval1 = smooth_minnum(end-1):smooth_maxnum(end-1);
    nval2 = smooth_minnum(end):smooth_maxnum(end);
    nval_common = intersect(nval1,nval2) - idx_adj;
    nval_common = nval_common(1:end-1);
    p_cum_dev{end} = [nval_common' + idx_adj, (pcum_time(nval_common,end) - pcum_time(nval_common,end-1)) / dt];
    % derivative at Second, third, ... time indices
    for idx = 1:lent - 3
        fprintf("%3d, [%d %d]",idx,1,lent - 3)
        nval1 = smooth_minnum(idx)  :smooth_maxnum(idx);
        nval2 = smooth_minnum(idx+1):smooth_maxnum(idx+1);
        nval3 = smooth_minnum(idx+2):smooth_maxnum(idx+2);
        nval4 = smooth_minnum(idx+3):smooth_maxnum(idx+3);
        nval_common = intersect(nval1,nval2);
        nval_common = intersect(nval_common,nval3);
        nval_common = intersect(nval_common,nval4);
        nval_common = nval_common(1:end-1) - idx_adj;
        p_cum_dev{idx+1} = zeros(length(nval_common),2);
        if idx == lent - 3
            p_cum_dev{idx+2} = zeros(length(nval_common),2);
        end
        for idx1=1:length(nval_common)
            pval = pcum_time(nval_common(idx1),idx:idx+3);
            tval = (idx:idx+3).*dt;
            pvalf = csapi(tval, pval);
            dpval = fnder(pvalf,1);
            dpvalf = @(t)fnval(dpval,t);
            p_cum_dev{idx+1}(idx1,1) = nval_common(idx1) + idx_adj;
            p_cum_dev{idx+1}(idx1,2) = dpvalf(tval(2));
            if idx == lent - 3
                p_cum_dev{idx+2}(idx1,1) = nval_common(idx1) + idx_adj;
                p_cum_dev{idx+2}(idx1,2) = dpvalf(tval(3));
            end
        end
        for bb = 1:12
            fprintf("\b")
        end
    end
    %% Additional smoothing
    if ~isempty(smooth_windows)
        % without this process, mean_theory calculated from Calculate_Mean_Var_In_situ.m might be incorrect
        smooth_meannumberdata = smoothdata(smooth_meannumberdata, 'gaussian', smooth_windows(1));
        smooth_varnumberdata = smoothdata(smooth_varnumberdata, 'gaussian', smooth_windows(2));
        smooth_smdata = smooth_varnumberdata + smooth_meannumberdata.^2;
        smooth_mean_pp = csapi(timedata, smooth_meannumberdata);
        smooth_mean_derv = fnder(smooth_mean_pp, 1);
        smooth_dndtmean = fnval(smooth_mean_derv, timedata);
    end
    %% Additional preprocessing for fast computation of Calculate_Cost_In_situ_fast.m
    nval_c = cell(lent, 1);
    pvalue_c = cell(lent, 1);
    pvalue_1_c = cell(lent, 1);
    pcum_val_c = cell(lent, 1);
    ptimedev_c = cell(lent, 1);
    for idx = 1:lent
        p_cum_dev_ = p_cum_dev{idx};
        pcum_ = pcumval{idx};
        nval = p_cum_dev_(:, 1);
        ptimedev = p_cum_dev_(:, 2);
        if options.trim~=0
            trim = options.trim;
            nval = nval(trim:end - trim);
            ptimedev = ptimedev(trim:end - trim);
        end

        idx1 = find(pcum_(:, 1) == nval(1));
        idxList = idx1-1:idx1 + length(nval);
        % pvalue = pcum_(idx1:idx1 + length(nval), 2);
        pvalue = pcum_(idxList(2:end), 2);
        pcum_val = pvalue(1:end - 1);

        if idx1 == 1 % pcum_(1,1) == nval(1)
            pvalue_1 = pcum_(idxList(2:end-1), 2);
            pvalue_1 = [0; pvalue_1];
        else
            pvalue_1 = pcum_(idxList(1:end-1), 2);
        end

        pvalue = -diff(pvalue); %p(n, t)
        pvalue_1 = -diff(pvalue_1); % p(n-1, t)

        nval_c{idx} = nval;
        pvalue_c{idx} = pvalue; %p(n, t)
        pvalue_1_c{idx} = pvalue_1; % p(n-1, t)
        pcum_val_c{idx} = pcum_val; % \psi_n^{p}(t) = \Sigma_{j=n}^{\infty}p(j, t)
        ptimedev_c{idx} = ptimedev; % \partial_t\psi_n^p(t)
    end
    
    expdata.timedata = timedata;
    expdata.radidata = radidata;
    expdata.numberdata = numberdata;
    expdata.maxnum = maxnum;
    expdata.minnum = minnum;
    expdata.meannumberdata = meannumberdata;
    expdata.varnumberdata = varnumberdata;
    expdata.smooth_numberdata = smooth_numberdata;
    expdata.smooth_meannumberdata = smooth_meannumberdata;
    expdata.smooth_varnumberdata = smooth_varnumberdata;
    expdata.smooth_dndtmean = smooth_dndtmean;
    expdata.smooth_smdata = smooth_smdata;
    expdata.p_cum_dev = p_cum_dev;
    expdata.pcumval = pcumval;
    expdata.pntdata = pntdata;
    expdata.nval_c = nval_c;
    expdata.pvalue_c = pvalue_c;
    expdata.pvalue_1_c = pvalue_1_c;
    expdata.pcum_val_c = pcum_val_c;
    expdata.ptimedev_c = ptimedev_c;
    
    save(filename, '-struct', 'expdata')
end