%% This script finds pareto front using genetic algorithm
restoredefaultpath
clear
addpath(fullfile("..","functions"))
is_test = false;
situ = "in";
expdata_original = load("expdata_PtCOD.mat");

stt_idx = 1;
reduced_factor = 1.1;
sigma_s = 265.077.^(-1/3);
dirName = fullfile("..","Pt(acac)2_monomeric","??"); % Update directory name to which fitdata for Pt(acac)2 is saved
PtA = load(fullfile(dirName, "GamultiobjOptions????????????.mat")); % Update filename to which fitdata for Pt(acac)2 is saved

lbub = [-100,0; ...lnqfqb
        -100,700; ...lnqeqb
        4,inf; ...alpha+4
        0,inf; ...kappa_a*rho_{1,inf}
        0,inf; ...kappa_a*sigma_s/D_1
        sigma_s, sigma_s; ... sigma_s
        0,inf; ...supersaturation ratio
        1,2 ... shape index
        ]';
if ~isempty(PtA)
    [x, ~] = pick_param(PtA, 1e-5);
    lbub(:,1:5) = [x(1:5);x(1:5)];
end

if is_test == true
    bs = 18;
    dirName = "test_results";
    fileName = string(datetime("now",Format = "yyMMdd")) + "test";
    popSize = 2^4;
    MaxGen = 2e1;
    FunTol = 1e-5;
else
    bs = 1;
    dirName = "?cpu";
    fileName = [];
    popSize = 2^8;
    MaxGen = 1e4;
    if ~isempty(PtA)
        popSize = PtA.optimopt.PopulationSize;
        MaxGen = PtA.optimopt.MaxGenerations;
    end
    FunTol = 1e-6;
end
expdata = sampling_expdata(expdata_original, stt_idx, bs);
%% Initial Setting
galbub = [-6, -1; -6  ,6  ; lbub(1,3),10; 0,10; 1e-4,1e-1; lbub(1,6),lbub(2,6); 0,100;  1, 2]'; % non-spherical theory

if ~isempty(PtA)
    galbub(:,1:6) = lbub(:,1:6);
end
lent = numel(expdata.timedata);
weights = ones(1, lent);
numStartPoints = 1;
CCCCCCCC = {1, 2, 3, 4, 5, 6, 7, 8};
CCCCCCTT = {1, 2, 3, 4, 5, 6, 7:6+lent, 7+lent:6+2*lent};

other_opt = struct('sigma_s',sigma_s, ...
            'indices',{CCCCCCTT}, ...
            'factor', reduced_factor, ...
            'loss',"mse", ...
            'scale',"minmax", ...
            'collect',"mean", ...
            'loss_var',"mse", ...
            'scale_var',"minmax", ...
            'weights',weights, ...
            'track',"latest", ...
            'startPoints',numStartPoints);
if ~isempty(PtA)
    other_opt.Pt_acac_2data = fullfile(PtA.dirName, PtA.filename);
else
    other_opt.Pt_acac_2data = [];
end
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
if isempty(PtA)
    gamulti_opt.CreationFcn = @(gl, ff, o)gacuGuess2(gl, ff, o, expdata, other_opt, situ);
else
    gamulti_opt.CreationFcn = [];
end
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