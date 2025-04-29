%% This script finds pareto front using genetic algorithm
restoredefaultpath
clear
clc
addpath(fullfile("..","functions"))
is_test = false;
target_NP = "Au";
situ = "in";
switch target_NP
    case "Au"
        stt_idx = 1;
        reduced_factor = 1;
        rho_s = 59.0088;    % moAnomer densadfadity in atom / nm^3
        sigma_s = (3 / (4 * pi * rho_s))^(1/3);
        lbub = [-100,0; ...lnqfqb
                -100,700; ...lnqeqb
                4,inf; ...alpha+4
                0,inf; ...kappa_a*rho_{1,inf}
                0,inf; ...kappa_a*sigma_s/D_1
                sigma_s, sigma_s; ... sigma_s
                0,inf; ...supersaturation ratio
                1,2 ... shape index 
                ]';
        if is_test == true
            bs = 40;
        else
            bs = 1;
        end
        expdata_original = load("expdata_Au.mat");
end
if is_test == true
    dirName = "test_results";
    fileName = string(datetime("now",Format = "yyMMdd")) + "test";
    popSize = 2^4;
    MaxGen = 2e1;
    FunTol = 1e-5;
else
    dirName = "?cpu";
    fileName = [];
    popSize = 2^8;
    MaxGen = 1e4;
    FunTol = 1e-6;
end
expdata = sampling_expdata(expdata_original, stt_idx, bs);
lent = numel(expdata.timedata);
%% Initial Setting  
galbub = [-6,-1; ...lnqfqb
           -6,6; ...lnqeqb
            lbub(1, 3),10; ...alpha+4
            0,10; ...kappa_a*rho_{1,inf}
            1e-4,1; ...kappa_a*sigma_s/D_1
            sigma_s, sigma_s; ... sigma_s
            0,100; ...supersaturation ratio
            1,2 ... shape index
            ]';
CCCCCCCC = {1, 2, 3, 4, 5, 6, 7, 8};
CCCCCCTT = {1, 2, 3, 4, 5, 6, 7:6+lent, 7+lent:6+2*lent};

other_opt = struct('sigma_s',sigma_s, ...
            'indices',{CCCCCCTT}, ...
            'factor', reduced_factor, ...
            'loss',"mse", ...
            'scale',"minmax", ...
            'weights',linspace(20,1, lent), ...linspace(20,1, lent)
            'collect',"mean", ...
            'loss_var',"mse", ...
            'scale_var',"minmax", ...
            'track',"latest", ...
            'startPoints',1);
ga_opt = optimoptions("ga", ...
    CreationFcn = @(g, f, o) gacuGuess2(g, f, o, expdata, other_opt, situ), ...
    Display = "iter", ...
    MaxGenerations = MaxGen, ...
    FunctionTolerance = 1e-10, ...
    InitialPopulationMatrix = [], ...
    PopulationSize = popSize, ...
    UseParallel = true);
gamulti_opt = optimoptions("gamultiobj", ...
    Display = "iter", ...
    MaxGenerations = MaxGen, ...
    PopulationSize = popSize, ...
    FunctionTolerance = FunTol, ...
    UseParallel = true);
%% (gamultiobj) (CCCCCCTT) Fit Jn & Variance simultaneousely 
other_opt.track = "every_100";
other_opt.indices = CCCCCCTT;
gamulti_opt.CreationFcn = @(gl, ff, o)gacuGuess2(gl, ff, o, expdata, other_opt, situ);
[jnvarfit1, jnvarcost1, jnvarresults1] = optim_and_save_in( ...
    @(x)[Calculate_Cost_In_situ_fast(x, expdata, other_opt); Var_loss(x, expdata, other_opt)],...
    [], ...
    [], ...
    change_TC(galbub(1,:), CCCCCCCC, other_opt.indices, []), ...
    change_TC(galbub(2,:), CCCCCCCC, other_opt.indices, []), ...
    gamulti_opt, other_opt, dirName, fileName);
function cost = Var_loss(x, expdata, options)
    [~, var_theory] = Calculate_Mean_Var_In_situ(x, expdata, options);
    var_exp = options.factor*expdata.smooth_varnumberdata;
    loss = options.loss_var;
    scale = options.scale_var;
    cost = lossReg(var_theory, var_exp, loss, scale);
end
