function [timedata_1, coal_exp_radidata, coal_theory_radidata, afit, b, ...
    mean_coal, ycalc] =  Smoluchowski_coagulation_kinetics(data)
coal_radidata = data.coal_radidata;
coal_time = data.coal_time;
init_p = 60;
final_p = 150;
mean_coal = mean(coal_radidata,2,'omitnan');
logmean = log(mean_coal);
logtime = log(1:180);
X = [ones(length(logtime(init_p:final_p)),1),logtime(init_p:final_p)'];
b = X\logmean(init_p:final_p);
ycalc = X*b;
mincoal_time = min(coal_time);
maxcoal_time = max(coal_time);
timedata_1 = mincoal_time:maxcoal_time;
afit = zeros(1,length(timedata_1));
nbin = 10;
idx1 = 1;
options = optimoptions('lsqnonlin', 'Display', 'none');
xval = 0:0.01:1.8;
coal_exp_radidata = cell(1,length(timedata_1));
coal_theory_radidata = cell(1,length(timedata_1));
for idx = timedata_1
    radi_data = coal_radidata(idx, :) ./ mean_coal(idx);
    [value, edges] = histcounts(radi_data, nbin, 'Normalization', 'pdf');
    edges = edges(1:end-1) + (edges(2)-edges(1))/2;
    fun = @(a) Smoluchowski_radius_dist(a, edges) - value;
    afit(idx1) = lsqnonlin(fun, 1, 0, inf, options);
    coal_exp_radidata{idx1} = [edges', value'];
    coal_theory_radidata{idx1} = [xval', Smoluchowski_radius_dist(afit(idx1), xval)'];
    idx1 = idx1 + 1;
end
% figure(1)
% subplot(2,2,1)
% plot(logtime, logmean, 'bo', logtime(init_p:final_p), logmean(init_p:final_p), 'go', logtime(init_p:final_p), ycalc, 'r-')
% hold on
% plot(logtime, logmean_mono, 'ko', logtime(init_p:final_p), logmean_mono(init_p:final_p), 'go', logtime(init_p:final_p), ycalc1, 'r-')
% savedata = [timedata, mean_coal, mean_mono];
% save('scaled_mean_Radii_coal_mono.dat', 'savedata', '-ascii')
% savedata = [exp(logtime(init_p:final_p))', exp(ycalc), exp(ycalc1)];
% save('fitting_scaled_mean_Radii_coal_mono.dat','savedata', '-ascii')