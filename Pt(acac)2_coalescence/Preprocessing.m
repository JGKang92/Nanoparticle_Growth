clear
coal_data = readmatrix('coalescence_pixel.csv');
pixel_len = 0.53;           % pixel size in nm
%% Radius data for coalescence
lent = size(coal_data,1);
timedata = coal_data(:,1);
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
coal_rel_xpos = zeros(lent,tot_coal_eventnum); 
coal_rel_ypos = zeros(lent,tot_coal_eventnum);
coal_rel_dist = zeros(lent,tot_coal_eventnum);
coal_data(coal_data == 0) = nan;
idx2 = 1;
for idx=1:numtrajcoal
    ncoal = coal_eventnum(idx);
    coalidx = nancolidx(idx);
    for idx1=1:ncoal
        coal_radidata(:,idx2) = coal_data(:,coalidx+1+(idx1-1)*3);
        coal_xdata(:,idx2) = (coal_data(:,coalidx+2+(idx1-1)*3) + 0.19140977611823007.*(1:lent)')*pixel_len;
        coal_ydata(:,idx2) = (coal_data(:,coalidx+3+(idx1-1)*3) + 0.40770301236991163.*(1:lent)')*pixel_len;
        idx2 = idx2 + 1;
    end
end
coal_xdata(coal_xdata == 0) = nan;
coal_ydata(coal_ydata == 0) = nan;
%
idx2 = 1;
idx3 = 1;
for idx=1:numtrajcoal
    ncoal = coal_eventnum(idx);
    coalidx = nancolidx(idx);
    xdata = zeros(lent,ncoal);
    ydata = zeros(lent,ncoal);
    for idx1=1:ncoal
        xdata(:,idx1) = (coal_data(:,coalidx+2+(idx1-1)*3) + 0.19140977611823007.*(1:lent)')*pixel_len;
        ydata(:,idx1) = (coal_data(:,coalidx+3+(idx1-1)*3) + 0.40770301236991163.*(1:lent)')*pixel_len;
    end
    idx2 = idx2 + 1;
    for idx1=1:ncoal-1
        sval = xdata(8:end,idx1) == xdata(8:end,idx1+1);
        sidx_all = find(sval == 1);
        for idx4=1:length(sidx_all)
            sdiff = diff(sidx_all(idx4:end));
            if sum(sdiff ~= 1) == 0
                sidx = sidx_all(idx4);
                break
            end
        end
        init_idx = xdata(8:end,idx1) ~= 0 & xdata(8:end,idx1+1) ~= 0 & ~isnan(xdata(8:end,idx1));
        init_idx = find(init_idx == 1,1) + 7;
        coal_rel_dist(:,idx3) = sqrt((xdata(:,idx1) - xdata(:,idx1+1)).^2 + (ydata(:,idx1) - ydata(:,idx1+1)).^2);
        coal_rel_xpos(:,idx3) = xdata(:,idx1) - xdata(:,idx1+1);
        coal_rel_ypos(:,idx3) = ydata(:,idx1) - ydata(:,idx1+1);
%         coal_time(idx3) = sidx+7 - (init_idx - 8);
        coal_time(idx3) = sidx + 7;
        coal_size(1,idx3) = coal_data(sidx+6, coalidx+1+(idx1-1)*3);
        coal_size(2,idx3) = coal_data(sidx+6, coalidx+1+(idx1)*3);
        coal_size(3,idx3) = coal_data(sidx+7, coalidx+1+(idx1-1)*3);
        coal_before_radi(idx3) = sqrt(coal_rel_xpos(sidx+6,idx3).^2 + coal_rel_ypos(sidx+6,idx3).^2);
        coal_radidata(sidx+7:end,idx2) = 0;
        idx2 = idx2 + 1;
        idx3 = idx3 + 1;
    end
end
coal_radidata = pixel_len*sqrt((coal_radidata) ./ pi); % radius
coal_radidata(coal_radidata == 0) = nan;
coal_size = pixel_len*sqrt(coal_size./pi);
coal_radius = coal_size(1,:) + coal_size(2,:);
%% Coalescence data modification
coal_xdata_modified = coal_xdata;
coal_ydata_modified = coal_ydata;
coal_radidata_modified = coal_radidata;
idx2=1;
idx3=1;
for idx=1:numtrajcoal
    ncoal = coal_eventnum(idx);
    for idx1=1:ncoal-1
        coal_idx = coal_time(idx2);
        coal_xdata_modified(coal_idx:end,idx3+idx1) = nan;
        coal_ydata_modified(coal_idx:end,idx3+idx1) = nan;
        coal_radidata_modified(coal_idx:end,idx3+idx1) = nan;
        idx2 = idx2 + 1;
        idx3 = idx3 + 1;
    end
    idx3 = idx3 + 1;
end
%% Monomeric growth data
monomeric_data = readmatrix('../Pt(acac)2_monomeric/monomeric_pixel.csv');
monomeric_radidata = monomeric_data(:,3:4:end);
monomeric_xdata = monomeric_data(:,4:4:end);
monomeric_ydata = monomeric_data(:,5:4:end);
monomeric_radidata = pixel_len*sqrt(monomeric_radidata./pi);
monomeric_radidata(monomeric_radidata == 0) = nan;
monomeric_xdata(monomeric_xdata == 0) = nan;
monomeric_ydata(monomeric_ydata == 0) = nan;
monomeric_xdata = (monomeric_xdata+0.19140977611823007.*(1:lent)')*pixel_len;
monomeric_ydata = (monomeric_ydata+0.40770301236991163.*(1:lent)')*pixel_len;
mean_monomeric_radius = mean(monomeric_radidata,2,'omitnan');
%% Total data 
total_xdata = [monomeric_xdata, coal_xdata_modified];
total_ydata = [monomeric_ydata, coal_ydata_modified];
total_radidata = [monomeric_radidata, coal_radidata_modified];
mean_total_radius = mean(total_radidata,2,'omitnan');
total_xdata_modified = total_xdata;
total_ydata_modified = total_ydata;
for idx=1:length(total_xdata)
    nidx1= find(~isnan(total_xdata(:,idx)) == 1, 1);
    total_xdata_modified(:,idx) = total_xdata(:,idx) - total_xdata(nidx1,idx);
    total_ydata_modified(:,idx) = total_ydata(:,idx) - total_ydata(nidx1,idx);
end
msd_total = (mean(total_xdata_modified.^2,2,'omitnan') + mean(total_ydata_modified.^2,2,'omitnan'))./(4*(1:180)');
%% x position displacement, y position displacement
tot_pair_num = length(total_xdata).*(length(total_xdata)-1)/2;
xdisp_pair = nan(tot_pair_num,173);
ydisp_pair = nan(tot_pair_num,173);
radius_pair = nan(tot_pair_num,173);
% xdisp_pair and ydisp_pair
for idx2=1:173
    idx3 = 1;
    for idx=1:length(total_xdata)
        nidx1 = find(~isnan(total_xdata(:,idx)) == 1, idx2);
        if length(nidx1) ~= idx2
            idx3 = idx3 + length(total_xdata) - idx;
            continue
        elseif length(nidx1) > 1
            nidx1 = nidx1(end);
        end
        for idx1=idx+1:length(total_xdata)
            nidx2 = find(~isnan(total_xdata(:,idx1)) == 1, idx2);
            if length(nidx2) ~= idx2
                idx3 = idx3 + 1;
                continue
            elseif length(nidx2) > 1
                nidx2 = nidx2(end);
            end
            nidx = max(nidx1,nidx2);
%             xdisp_pair{:,idx2} = [xdisp_pair{:,idx2}, (total_xdata(nidx,idx) - total_xdata(nidx,idx1))];
%             ydisp_pair{:,idx2} = [ydisp_pair{:,idx2}, (total_ydata(nidx,idx) - total_ydata(nidx,idx1))];
%             radius_pair{:,idx2} = [radius_pair{:,idx2}, total_radidata(nidx,idx) + total_radidata(nidx,idx1)];
            xdisp_pair(idx3,idx2) = total_xdata(nidx,idx) - total_xdata(nidx,idx1);
            ydisp_pair(idx3,idx2) = total_ydata(nidx,idx) - total_ydata(nidx,idx1);
%             radius_pair(idx3,idx2) = total_radidata(nidx,idx) + total_radidata(nidx,idx1);
            idx3 = idx3 + 1;
        end
    end
end
% radius_pair part
for idx=1:173
    idx3=1;
    for idx1=1:length(total_xdata)
        for idx2=idx1+1:length(total_xdata)
            radius_pair(idx3,idx) = total_radidata(idx+7,idx1) + total_radidata(idx+7,idx2);
            idx3 = idx3+1;
        end
    end
end
%% Initial r0 distribution
% tidx = floor(linspace(1,173,64));
% cvec = hsv(173);
% figure(1)
% clf
% for idx=1:64
% %    tidxval = tidx(idx);
%    r0data = sqrt(xdisp_pair(:,tidxval).^2 + ydisp_pair(:,tidxval).^2);
%    r0data = r0data(~isnan(r0data));
%    [N, edges]=histcounts(r0data,150);
%    [N_norm, ~] = histcounts(r0data,150,'Normalization','pdf');
%    edges_new = edges(2:end) - (edges(2) - edges(1))/2;
%    subplot(1,2,1)
%    plot(edges_new,movmean(N_norm,10), '-', 'Color', cvec(tidxval,:))
%    hold on
%    subplot(1,2,2)
%    plot(edges_new,movmean(N_norm,10)./(2.*pi*edges_new), '-', 'Color', cvec(tidxval,:))
% %    title(sprintf('%d sec',tidxval))
%     hold on
% end
% hold off
%%
rhotinit = sqrt(xdisp_pair(:,1).^2 + ydisp_pair(:,1).^2);
rhotinit = rhotinit(~isnan(rhotinit));
[N_norm, edges_norm] = histcounts(rhotinit,150,'Normalization','pdf');
edges_new = edges_norm(2:end) - (edges_norm(2) - edges_norm(1))/2;
rhotdata = (movmean(N_norm,10)./(2.*pi*edges_new))./(2.265e-5);
% figure(1)
% clf
% plot(edges_new,rhotdata, '-')
% save('rhot_init.mat','rhotdata','edges_new')
%% pair of within critical distance, 7.6810
%crit_pair_num  = zeros(1,173);
%for idx=1:173
%    crit_pair_num(idx) = sum(sqrt(xdisp_pair(:,idx).^2 + ydisp_pair(:,idx).^2) < 7.6810);
%end
%figure(2)
%clf
%plot(8:lent, crit_pair_num)
%% Survival probability
fst_nan = zeros(1,184);
AtA0_exp = zeros(1,173);
ntraj_obs = zeros(1,173);
n0val = sum(~isnan(total_radidata(8,:)));
for idx=1:184
    fst_nan(idx) = find(~isnan(total_radidata(:,idx)) == 1, 1);
end
for idx=1:173
    ntraj_obs(idx) = sum(fst_nan <= idx + 7);
    AtA0_exp(idx) = sum(~isnan(total_radidata(idx+7,:)));
end
AtA0_exp = AtA0_exp ./ AtA0_exp(1);
% figure(3)
% clf
% subplot(1,2,1)
% plot(8:lent, ntraj_obs)
% subplot(1,2,2)
% plot(8:180, survt)
%figure(3)
%plot(8:lent,survt, '-')
%save('AtA0value.dat','survt','-ascii')
% save('ntraj_obs.dat', 'ntraj_obs', '-ascii')
%% sigma(t) behavior
sigmat = mean(radius_pair,'omitnan');
% figure(5)
% clf
% plot(8:lent, sigmat, 'b-', 1:lent, 2*mean_monomeric_radius, 'r-', 1:lent, 2*mean_total_radius, 'g--')
% sigmat_save = zeros(173,2);
% sigmat_save(:,1) = 8:lent;
% sigmat_save(:,2) = sigmat;
% save('sigmat_20210128.dat','sigmat_save','-ascii')
%% Dr(t)
delta_x_pair = xdisp_pair - xdisp_pair(:,1);
delta_y_pair = ydisp_pair - ydisp_pair(:,1);
delta_x_2_pair = mean(delta_x_pair.^2,'omitnan') - mean(delta_x_pair,'omitnan').^2;
delta_y_2_pair = mean(delta_y_pair.^2,'omitnan') - mean(delta_y_pair,'omitnan').^2;
dt_disp_pair = (delta_x_2_pair + delta_y_2_pair) ./ (4*(1:173));
%figure(6)
%clf
%subplot(1,2,1)
%plot(8:lent, dt_disp_pair)
%subplot(1,2,2)
%loglog(8:lent, dt_disp_pair)
%dt_save = zeros(173,2);
%dt_save(:,1) = 0:172;
%dt_save(:,2) = dt_disp_pair;
%save('dt_20210129.dat','dt_save','-ascii')
%% Coalescence event time 0 analysis
% surface to surface
coal_surf_dist = zeros(lent,tot_coal_eventnum);
idx2 = 1;
idx3 = 1;
for idx=1:numtrajcoal
    ncoal = coal_eventnum(idx);
    for idx1=1:ncoal-1
        radius_sum = coal_radidata(:,idx3) + coal_radidata(:,idx3+1);
        if sum(~isnan(radius_sum)) == 0
            radius_sum = coal_radidata(:,idx3-1) + coal_radidata(:,idx3+1);
        end
        coal_surf_dist(:,idx2) = coal_rel_dist(:,idx2) - radius_sum;
        idx2 = idx2 + 1;
        idx3 = idx3 + 1;
    end
    idx3 = idx3 + 1;
end
coal_surf_dist_modified = nan(lent,tot_coal_eventnum);
for idx=1:tot_coal_eventnum
    nanidx = ~isnan(coal_surf_dist(:,idx));
    fidx = find(nanidx==1);
    coal_surf_dist_modified(lent-(fidx(end)-fidx(1)):end,idx) = coal_surf_dist(nanidx,idx);
end
%figure(7)
%plot(1:lent, coal_surf_dist_modified, '-')
%% Number data
monomeric_number = round(265.077.*monomeric_radidata.^3);
coal_number = round(265.077.*coal_radidata.^3);
monomeric_rel_std = sqrt(var(monomeric_number,1,2,'omitnan'))./mean(monomeric_number,2,'omitnan');
coal_rel_std = sqrt(var(coal_number,1,2,'omitnan'))./mean(coal_number,2,'omitnan');
coal_rel_std = coal_rel_std(8:180);
% Total number distribution
numberdata = [monomeric_number(8:180,:), coal_number(8:180,:)];
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
%% Save variables
save('Ptacac2_coalescence.mat', 'coal_rel_std', 'dt_disp_pair', 'AtA0_exp', 'sigmat', 'coal_time', ...
    'total_radidata', 'pntdata', 'timedata', 'rhotdata', 'edges_new', 'coal_radius', 'ntraj_obs', 'coal_radidata')