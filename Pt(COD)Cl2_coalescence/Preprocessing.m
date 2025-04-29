clear
coal_data = readmatrix('coalescence_growth_cod.csv');
[lent,colnum] = size(coal_data);
%%
% Radius data for coalescence
timedata = 2*(1:lent)';
nanidx = isnan(coal_data(1,:));
numtrajcoal = sum(nanidx);
nancolidx = 1:length(nanidx);
nancolidx = nancolidx(nanidx);
nancolidx = [nancolidx, length(nanidx)+1];
coal_eventnum = (diff(nancolidx)-1)/3;
tot_coal_eventnum = sum(coal_eventnum-1);
coal_numtrajobs = sum(coal_eventnum);
coal_radidata = zeros(lent,coal_numtrajobs);
coal_xdata = zeros(lent,coal_numtrajobs);
coal_ydata = zeros(lent,coal_numtrajobs);
coal_time = zeros(1,tot_coal_eventnum);
coal_size = zeros(3,tot_coal_eventnum);
coal_data(coal_data == 0) = nan;
idx2 = 1;
for idx=1:numtrajcoal
    ncoal = coal_eventnum(idx);
    coalidx = nancolidx(idx);
    for idx1=1:ncoal
        coal_radidata(:,idx2) = coal_data(:,coalidx+1+(idx1-1)*3);
        coal_xdata(:,idx2) = (coal_data(:,coalidx+2+(idx1-1)*3) + 0.19140977611823007.*(1:lent)')*0.53;
        coal_ydata(:,idx2) = (coal_data(:,coalidx+3+(idx1-1)*3) + 0.40770301236991163.*(1:lent)')*0.53;
%         if ncoal == 2
%             coal_radidata_1(:,idx0) = coal_radidata(:,idx2);
%             idx0 = idx0 + 1;
%         end
        idx2 = idx2 + 1;
    end
end
coal_radidata = sqrt((0.53*0.53*coal_radidata) ./ pi); % radius
coal_radidata(coal_radidata == 0) = nan;
coal_number = round(265.077.*coal_radidata.^3);
coal_rel_std = sqrt(var(coal_number,1,2,'omitnan'))./mean(coal_number,2,'omitnan');
%% Number data
nanoparticle_data = readmatrix(fullfile('..","Pt(COD)Cl2","Monomeric_2.csv'));
monomeric_radidata = sqrt((0.53*0.53*nanoparticle_data(:,1:3:end)) ./ pi);
monomeric_number = round(265.077.*monomeric_radidata.^3);
coal_number = round(265.077.*coal_radidata.^3);
% Total number distribution
numberdata = [monomeric_number, coal_number];
numberdata(numberdata == 0) = nan;
[lent, ntrj] = size(numberdata);
smooth_numberdata = nan(lent, ntrj);
smooth_dndtdata = nan(lent, ntrj);
% FFT filter method
for idx=1:ntrj
    ridx = ~isnan(numberdata(:,idx)); % find the time indices where idx-th trajectory of numberdata is nan
    rnum = sum(ridx); % number of time points that are not nan
    rdata = numberdata(ridx,idx); 
    ftrdata = [rdata(end:-1:1);rdata;rdata(end:-1:1)];
    tdata = timetable(seconds((1:length(ftrdata))'), ftrdata);
    smooth_tdata = lowpass(tdata, 0.01);
    sm_data = smooth_tdata.ftrdata(rnum+1:2*rnum);
    smooth_numberdata(ridx, idx) = sm_data;
    smpp = csapi(timedata(ridx),sm_data);
    smfprime = fnder(smpp,1);
    smooth_dndtdata(ridx,idx) = fnval(smfprime, timedata(ridx));
end
smooth_numberdata = round(smooth_numberdata);
smooth_maxnum = max(smooth_numberdata,[],2,'omitnan');
smooth_minnum = min(smooth_numberdata,[],2,'omitnan');
smooth_meannumberdata = mean(smooth_numberdata,2,'omitnan');
smooth_varnumberdata = var(smooth_numberdata,1,2,'omitnan');
smooth_smdata = smooth_varnumberdata + smooth_meannumberdata.^2;
smooth_dndtmean = mean(smooth_dndtdata,2,'omitnan');
%% Obtain cdf P(n <= z, t) and pdf P(n = z, t) for nbin number of ranges
nbin = 30;
edges = zeros(nbin + 1, lent);
value = zeros(nbin, lent);
smooth_value = zeros(nbin, lent);
pntdata = cell(lent, 1);
for idx = 1:lent
    num_data = smooth_numberdata(idx,:);
    num_data = num_data(~isnan(num_data));
    % [value1, edge] = histcounts(num_data, nbin, 'Normalization', 'cdf'); % P(n <= z)
    % [value2,  ~  ] = histcounts(num_data, nbin, 'Normalization', 'pdf'); % P(n = z)
    % value(:,idx) = 1 - value1 + value2; % P(n >= z)
    [value(:,idx), edges(:,idx)] = histcounts(num_data, nbin, 'Normalization', 'pdf');
    % edges(:,idx) = edge;
    smooth_value(:,idx) = movmean(value(:,idx),5); % moving averaged pdf values
    pntdata{idx} = [edges(1:end-1,idx) + (edges(2,idx) - edges(1,idx)) / 2, value(:,idx)];
end
save('PtCODCl2_coalescence.mat', 'coal_rel_std', 'pntdata')