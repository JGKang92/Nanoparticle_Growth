restoredefaultpath
clear
clc
close all   
addpath(fullfile("..","functions"))

PtA = load(fullfile("..","Pt_acac_2.mat"));
PtC = load(fullfile("..","Pt_COD_Cl_2.mat"));
Au = load(fullfile("..","AuHCl_4.mat"));
Cd = load(fullfile("..","CdSe.mat"));
Fe = load(fullfile("..","Fe_xO_y.mat"));
PtE = load(fullfile("..","Pt_ex.mat"));
cmap = [0.1484    0.2852    0.5469;
        0.2227    0.2812    0.5273;
        0.2812    0.2812    0.5078;
        0.3320    0.2773    0.4844;
        0.3789    0.2734    0.4648;
        0.4219    0.2695    0.4414;
        0.4609    0.2656    0.4180;
        0.5078    0.2578    0.3945;
        0.5430    0.2500    0.3750;
        0.5781    0.2422    0.3516;
        0.6133    0.2305    0.3281;
        0.6484    0.2188    0.3086;
        0.6797    0.2070    0.2852;
        0.7148    0.1953    0.2617;
        0.7422    0.1758    0.2422;
        0.7773    0.1523    0.2188;
        0.8047    0.1289    0.1992;
        0.8320    0.0977    0.1758;
        0.8594    0.0508    0.1523;
        0.8828    0.0039    0.1328];
npList = [PtA, PtC, Au, PtE, Fe, Cd];
%% Bootstrapping
samplingNumber = 1e3;
for ii = 1:3
    np = npList(ii);
    % bootstrap sampling of mean
    bootstrap_mean = bootstrap_sampling(np.result.numberdata_with_error, ...
        @(x)mean(x), samplingNumber);
    np.result.bootstrap_mean = bootstrap_mean;
    % bootstrap sampling of variance
    bootstrap_var = bootstrap_sampling(np.result.numberdata_with_error, ...
        @(x)var(x,1), samplingNumber);
    np.result.bootstrap_var = bootstrap_var;
    % bootstrap sampling of higher order moments
    for order = 3:10
        bootstrap_order = bootstrap_sampling(np.result.numberdata_with_error, ...
            @(x)mean(x.^order), samplingNumber);
        fieldName = sprintf("bootstrap_%d",order);
        np.result.(fieldName) = bootstrap_order;
    end
    npList(ii) = np;
end
PtA = npList(1);
PtC = npList(2);
Au = npList(3);
%% Main Figures 
Figure1(PtA, PtC, Au, save = [false, true], cmap = cmap)
%% 
Figure2(PtA, PtC, Au, save = true)
%% scale 30%
Figure3(PtA, PtC, Au, save = true)
%%
Figure4(PtA, PtC, Au, save = true, cmap = cmap)
%%
Figure5(PtA, PtC, Au, PtE, Fe, Cd, save = true)
%%
function bootstrap_Func = bootstrap_sampling(numberdata, Func, samplingNumber, options)
    arguments
        numberdata
        Func
        samplingNumber
        options.seed = "default"
    end
    rng(options.seed)
    if samplingNumber > 0
        lent = size(numberdata,1);
        bootstrap_Func = zeros(lent, samplingNumber);
        for tt = 1:lent
            data = numberdata(tt,:);
            data = data(~isnan(data));
            numData = numel(data);
            for mm = 1:samplingNumber
                resample = data(randi(numData, 1, numData));
                bootstrap_Func(tt,mm) = Func(resample);
            end
        end
    end
end