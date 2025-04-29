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
is_save = true;
is_show = true;
%% Nanoparticle size distributions at various times for Group-A and the whole set including both Group-A and Group-B
plot_save_pnGroupAB(PtA, "ptacac", 101, is_save, is_show);
plot_save_pnGroupAB(PtC, "ptcod", 102, is_save, is_show);
ptacac2coaldata = load(fullfile("..","Pt(acac)2_coalescence","Ptacac2_coalescence.mat"));
if is_save
    varNames = ["Particle size", "Particle size distribution"];
    p = ptacac2coaldata.pntdata;
    for ii = 1:numel(p)
        t = PtA.expdata.timedata(ii);
        p_ = p{ii};
        writetable(array2table(p_, VariableNames = varNames), ...
            "ptacacA_pnt(S10).xls", Sheet = sprintf("%.1fsec",t))
        disp(ii)
    end
end
if is_save
    ptcodcoaldata = load(fullfile("..","Pt(COD)Cl2_coalescence","PtCODCl2_coalescence.mat"));
    p = ptcodcoaldata.pntdata;
    for ii = 1:numel(p)
        t = PtC.expdata.timedata(ii);
        p_ = p{ii};
        writetable(array2table(p_, VariableNames = varNames), ...
            "ptcodA_pnt(S10).xls", Sheet = sprintf("%.1fsec",t))
        disp(ii)
    end
end
%% Geometry dependence of nanoparticle size statistics
if is_show
    figtile(11);
end
plot_save_GeomDep(PtA, "ptacac.xls", [nexttile(1), nexttile(3)], is_save, is_show);
plot_save_GeomDep(PtC, "ptcod.xls",  [nexttile(2), nexttile(4)], is_save, is_show);
%% Rate coefficients of monomer association and dissociation
n = round(logspace(log10(20), 4, 1e2))';
x = PtA.fitdata.xfit;
% x(1:6) = ones(1, 6)*1000;
[kna_per_area, knd_per_area] = Calculate_kna_knd_per_area(x, n, PtA.expdata, PtA.fitdata.other_opt);
if is_show
    figtile(12);
    nexttile
    colororder({'r','b'})
    yyaxis left
    semilogx(n, kna_per_area,'-o', LineWidth = 2)
    yyaxis right
    semilogx(n, knd_per_area,'-o', LineWidth = 2)
    xlim tight
end
if is_save
    writetable(array2table([n, kna_per_area, knd_per_area], ...
            VariableNames = ["particle size","kna", "knd"]), ...
        "kn_per_area.csv")
end
%% Effect of ratio between number densities of n-mers and (n+1)-mers on nanoparticle growth rate
if is_save
    for ii = 1:numel(PtA.expdata.timedata)
        a = PtA.result.ratio_comparison{ii};
        tbl = array2table(a, ...
            VariableNames = ["Particle Size","approximation","original"]);
        writetable(tbl, "ratio_comparison(S13).xls", ...
            Sheet = sprintf("%.1fsec",PtA.expdata.timedata(ii)))
        disp(ii)
    end
end
if is_show
    f = figtile(13);
    for ii = 1:5:numel(PtA.result.ratio_comparison)
        nexttile
        p = PtA.result.ratio_comparison{ii};
        plot(p(:, 1), p(:, 2))
        hold on
        plot(p(:, 1), p(:, 3))
    end
end
%%  Time profile of the supersaturation ratio extracted from our analysis of the experimental data
if is_save
    figtile(14);
end
plot_save_Supsat(PtA, "ptacac", nexttile(1), is_save, is_show)
plot_save_Supsat(PtC, "ptcod", nexttile(2), is_save, is_show)
plot_save_Supsat(Au, "Au", nexttile(3), is_save, is_show)
plot_save_Supsat(Fe, "Fe", nexttile(4), is_save, is_show)
plot_save_Supsat(Cd, "Cd", nexttile(5), is_save, is_show)
plot_save_Supsat(PtE, "ptEx", nexttile(6), is_save, is_show)

%% Size-dependent growth rate and size distributions of iron oxide and CdSe nanoparticles at various times
plot_save_pnJn(Fe, "Fe", 151, [2, 3, 5, 10, 20], is_save, is_show)
plot_save_pnJn(Cd, "Cd", 152, [7, 9, 10, 18, 28], is_save, is_show)
plot_save_pnJn(PtE, "PtEx", 153, PtE.expdata.timedata, is_save, is_show)
%% Time profiles of scaled cross correlations between nanoparticle size and growth rate coefficient for in-situ nanoparticle growth data
if is_show
    figtile(19);
    ax1 = nexttile(1); ax1.NextPlot = "add";
    ax2 = nexttile(2); ax2.NextPlot = "add";
    ax3 = nexttile(3); ax3.NextPlot = "add";
else
    ax1 = [];
    ax2 = [];
    ax3 = [];
end
plot_save_njnCorr(PtA, "ptacac", ax1, is_save, is_show)
plot_save_njnCorr(PtC, "ptcod", ax2, is_save, is_show)
plot_save_njnCorr(Au, "Au", ax3, is_save, is_show)
%% Approximate theoretical results for the relative variance of the nanoparticle size, calculated with eq M12
if is_show
    f = figtile(20);
    f.Position(2:4) = [1, 1600, 700];
end
plot_save_RelvarApprox(PtA, "ptacac", is_save, is_show)
plot_save_RelvarApprox(PtC, "ptcod", is_save, is_show)
plot_save_RelvarApprox(Au, "Au", is_save, is_show)
plot_save_RelvarApprox(Fe, "Fe", is_save, is_show)
plot_save_RelvarApprox(Cd, "Cd", is_save, is_show)
%% Theoretical results for relative variance of nanoparticle size with the optimal parameter values doubled or halved.
plot_save_RelvarDoubledHalved(PtA, "ptacac", 281, "in", is_save, is_show)
plot_save_RelvarDoubledHalved(Au, "Au", 282, "in", is_save, is_show)
plot_save_RelvarDoubledHalved(Fe, "Fe", 283, "ex", is_save, is_show)
plot_save_RelvarDoubledHalved(Cd, "Cd", 284, "ex", is_save, is_show)

if is_show
    figtile(28, [4, 5]);
end
Plot_param_sens(PtA.expdata, PtA.result, 0)
Plot_param_sens(Au.expdata, Au.result, 5)
Plot_param_sens(Fe.expdata, Fe.result, 10)
Plot_param_sens(Cd.expdata, Cd.result, 15)
for ii = 1:20
    nexttile(ii)
    axis([0 inf 0 1])
end
%% Phase-dependent free energy change associated with the formation of an n-mer from n monomer
if is_show
    f = figtile(32, [2, 5]);
end
plot_save_PhaseDependentFreeEnergyChange(PtA, "ptacac", [17,27; 28,35; 41,80; 107,154; 155,180] - 7, is_save, is_show)
plot_save_PhaseDependentFreeEnergyChange(PtC, "ptcod", [1,21; 22,150; 151,168], is_save, is_show)
plot_save_PhaseDependentFreeEnergyChange(Au, "Au", [1,86;87,387], is_save, is_show)

%% Plot jn per unit area
if is_show
    figtile(1);
    nexttile(1)
    plot_jn_per_area(PtA, 1, 14)
    title("Pt(acac)_2")
    ylim([-0.9 3])
    
    nexttile(2)
    plot_jn_per_area(PtC, 2, 14)
    title("Pt(COD)Cl_2")
    ylim([-0.9 3])
    
    nexttile(3)
    plot_jn_per_area(Au, 2, 50)
    title("AuHCl_4")
    ylim([-0.9 3])
end
%% Export All Figures as .pdf format
figs_all = findall(groot, Type = "figure");
figsFileName = sprintf("Figures_SI_%s.pdf", datetime("now", Format = "yyMMddHHmm"));
[~, idx] = sort([figs_all.Number]);
figs_all = figs_all(idx);
for ii = 1:numel(figs_all)
    fig = figs_all(ii);
    if ii == 1
        a = false;
    else
        a = true;
    end
    exportgraphics(fig, figsFileName, Append = a)
    fprintf("%2d/%2d exported as pdf file\n",ii,numel(figs_all))
end
%%
function plot_save_pnGroupAB(np, filename, figNum, is_save, is_show)
    varNames = ["Particle size", "Particle size distribution"];
    if is_save
        writematrix([[nan, np.expdata.timedata(:)']; ...
                    [1:1:max(np.expdata.maxnum);np.result.pnt_guess]'], ...
            filename + "_pntguess.csv")
        for ii = 1:numel(np.expdata.timedata)
            t = np.expdata.timedata(ii);
            p = np.expdata.pntdata{ii};
            writetable(array2table(p, VariableNames = varNames), ...
                filename + "_pnt_exp.xls", Sheet = sprintf("%.1fsec",t))
            disp(ii)
        end
    end
    if is_show
        f = figtile(figNum);
        title(f.Children, filename)
        for ii = 1:10:numel(np.expdata.timedata)
            experiment = np.expdata.pntdata{ii};
            nexttile
            plot(experiment(:,1), experiment(:,2),'o')
            hold on
            pntguess = [1:1:max(np.expdata.maxnum);np.result.pnt_guess(ii,:)]';
            plot(pntguess(:, 1), pntguess(:, 2),'-')
        end
    end
end
function plot_save_GeomDep(np, filename, axes, is_save, is_show)
    varNames = ["Time", "Variance(experiment)", "CV(experiment)",...
        "Variance_CubOcta", "Variance_TruncatedOcta", "Variance_Octa", "Variance_1.8", ...
        "CV_CubOcta", "CV_TruncatedOcta", "CV_Octa","CV_1.8"];
    a = [np.expdata.timedata(:), ...
        np.expdata.varnumberdata, ...
        sqrt(np.expdata.varnumberdata)./np.expdata.meannumberdata, ...
        np.result.var_theory_cOh, ...
        np.result.var_theory_tOh, ...
        np.result.var_theory_Oh, ...
        np.result.var_theory_18, ...
        sqrt(np.result.var_theory_cOh)./np.result.mean_theory_cOh, ...
        sqrt(np.result.var_theory_tOh)./np.result.mean_theory_tOh, ...
        sqrt(np.result.var_theory_Oh)./np.result.mean_theory_Oh, ...
        sqrt(np.result.var_theory_18)./np.result.mean_theory_18, ...
        ];
    if is_save
        writetable(...
            array2table(a, ...
                VariableNames = varNames), ...
            filename, ...
            Sheet = "MeanVarShapeIndex", ...
            WriteMode = "overwrite")
    end
    if is_show
        ax1 = axes(1);
        ax1.NextPlot = 'add';
        plot(ax1, a(:, 1), a(:, 4),'b-', LineWidth = 2, DisplayName = "Cuboctahedron")
        plot(ax1, a(:, 1), a(:, 5),'r-', LineWidth = 2, DisplayName = "Truncated ctahedron")
        plot(ax1, a(:, 1), a(:, 6),'y-', LineWidth = 2, DisplayName = "Octahedron")
        plot(ax1, a(:, 1), a(:, 7),'g-', LineWidth = 2, DisplayName = "1.8")
        plot(ax1, a(:, 1), a(:, 2),'ko', LineWidth = 2, DisplayName = "Experiment")
        legend(Location = "northwest", Box = "off")
        ax2 = axes(2);
        ax2.NextPlot = 'add';
        plot(ax2, a(:,1), a(:, 3),'ko', LineWidth = 2)
        hold on
        plot(ax2, a(:, 1), a(:, 8),'b-', LineWidth = 2)
        plot(ax2, a(:, 1), a(:, 9),'r-', LineWidth = 2)
        plot(ax2, a(:, 1), a(:, 10),'y-', LineWidth = 2)
        plot(ax2, a(:, 1), a(:, 11),'g-', LineWidth = 2)
    end
end
function plot_save_Supsat(np, filename, axes, is_save, is_show)
    if isfield(np.expdata, "tval")
        t = np.expdata.tval;
    else
        t = np.expdata.timedata;
    end
    a = [t(:), np.fitdata.xfit(np.fitdata.other_opt.indices{7})'];
    if is_save
        varNames = ["Time", "Supersaturation ratio"];
        writetable(array2table(a, VariableNames = varNames), ...
            filename + ".xls", Sheet = "Supersaturation ratio")
    end
    if is_show
        plot(axes, a(:,1), a(:,2),'o')
    end
end
function plot_save_pnJn(np, filename, figNum, timePoints, is_save, is_show)
    t = np.expdata.timedata;
    lent = numel(t);
    % saving Jntilde.xls
    if is_save
        for ii = 1:lent
            fprintf("%3d/%3d\n",ii,lent)
            theory = np.result.Jntilde{ii};
            experiment = np.result.psi_n_val{ii};
            if ~isequal(theory(:, 1), experiment(:,1))
                error("something wrong in S15")
            end
            writetable(table(theory(:, 1), theory(:, 2), experiment(:, 2), ...
                    VariableNames = ["Particle Size","Jntilde(Theory)","Jntilde(experiment)"]), ...
                filename + "_Jntilde.xls", ...
                Sheet = sprintf("%.1fmin",t(ii)))
            for jj = 1:8
                fprintf("\b")
            end
        end
    end
    
    experiment = np.expdata.pntdata';
    if size(experiment,1) == numel(t)
        experiment = experiment;
    elseif size(experiment,2) == numel(t)
        experiment = experiment';
    else
        error("something wrong in pntdata")
    end
    if is_save % saving pnt.xls
        for ii = 1:lent
            fprintf("%3d/%3d\n",ii,lent)
            writetable(...
                table( (1:length(experiment))', experiment(ii,:)', ...
                    VariableNames = ["Particle Size","pnt(experiment)"]), ...
                filename + "_pnt.xls", ...
                Sheet = sprintf("%.1fmin(experiment)", t(ii)))
            writetable(...
                table( (1:max(np.expdata.maxnum))', np.result.pnt_guess(ii, :)', ...
                    VariableNames = ["Particle Size","pnt(experiment)"]), ...
                filename + "_pnt.xls", ...
                Sheet = sprintf("%.1fmin(theory)",t(ii)))
            for jj = 1:8
                fprintf("\b")
            end
        end
    end
    if is_show % Plotting part
        figtile(figNum);
        idxList = arrayfun(@(x)find(t == x),timePoints);
        for ii = 1:numel(idxList)
            idx = idxList(ii);
            nexttile
            theory = np.result.Jntilde{idx};
            experiment = np.result.psi_n_val{idx};
            plot(theory(:, 1), theory(:, 2),'r-')
            hold on
            skip = 100;
            plot(experiment(1:skip:end, 1), experiment(1:skip:end, 2),'ro')
            title(t(idx))
            ylim([-5 15].*1e-4)
        end
        for ii = 1:numel(idxList)
            idx = idxList(ii);
            nexttile
            skip = 10;
            experiment = np.expdata.pntdata;
            if size(experiment,1) == lent
                exp = experiment;
            elseif size(experiment, 2) == lent
                exp = experiment';
            else
                error("Something wrong...")
            end
            plot(1:skip:length(experiment), exp(idx,1:skip:end),'bs')
            hold on
            plot(1:max(np.expdata.maxnum), np.result.pnt_guess(idx, :),'b-')
            title(t(idx))
        end
    end
end
function plot_save_njnCorr(np, filename, axes, is_save, is_show)
    a1 = [np.expdata.timedata(:), np.result.deltan1_deltajn(:), np.result.dt_var_mean(:)];
    if is_show
        plot(axes, a1(:,1), a1(:,2),'o')
        hold on
        plot(axes, a1(:,1), a1(:,3),'r-', LineWidth = 2)
        yline(axes, 0, Color = 'k', Alpha = 1)
        ylim(axes, [-5 10])
    end
    if is_save
        varNames = ["Time", "symbol", "line"];
        writetable( ...
            array2table(a1, VariableNames = varNames), ...
            filename + ".xls", Sheet = "dndj_n")
    end
end
function plot_save_RelvarApprox(np, filename, is_save, is_show)
    t = np.expdata.timedata(:);
    a1 = [t, np.expdata.relvar(:)];
    relvar = np.result.eq_m12_relvar(:);
    t0idx = numel(t) - numel(relvar) + 1;
    a2 = [t(t0idx:end), relvar];
    if is_save
        varNames1 = ["Time","experiment"];
        varNames2 = ["Time","Approximate"];
        writetable(array2table(a1, VariableNames = varNames1), ...
            filename + ".xls", Sheet = "RelVarExperiment")
        writetable(array2table(a2, VariableNames = varNames2), ...
            filename + ".xls", Sheet = "RelVarApprox")
    end
    if is_show
        nexttile
        plot(a1(:,1), a1(:,2),'o')
        hold on
        plot(a2(:,1), a2(:,2), 'k--', LineWidth = 2)
        title(filename)
    end
    
end
function plot_save_RelvarDoubledHalved(np, filename, figNum, situ, is_save, is_show)
    varNames = ["Time","experiment","original", ...
        "ef_doubled", "ee_doubled", "a_doubled","kar1inf_doubled", "kassd1_doubled", ...
        "ef_halved", "ee_halved", "a_halved","kar1inf_halved", "kassd1_halved"];
    switch situ
        case "in"
            t = np.expdata.timedata;
            v = np.expdata.varnumberdata;
            m = np.expdata.meannumberdata;
            relvar = v./m.^2;
        case "ex"
            t = np.expdata.tval;
            v = np.expdata.varr3data;
            m = np.expdata.meanr3data;
            relvar = v./m.^2;
            if numel(t)~=numel(relvar)
                relvar_new = nan(numel(t),1);
                for ii = 1:numel(t)
                    tidx = find(t(ii)==np.expdata.timedata);
                    if isempty(tidx)
                        continue
                    else
                        relvar_new(ii) = relvar(tidx);
                    end
                end
            end
            relvar = relvar_new;
    end
    if is_save
        a1 = [t(:), relvar, ...
            np.result.var_theory./np.result.mean_theory.^2, ...
            cell2mat(np.result.var_theory_twice)./cell2mat(np.result.mean_theory_twice).^2, ...
            cell2mat(np.result.var_theory_half)./cell2mat(np.result.mean_theory_half).^2];
        writetable(array2table(a1, ...
            VariableNames = varNames),filename + ".xls",Sheet = "Doubled_Halved")
    end
    if is_show
        figtile(figNum);
        Plot_param_sens(np.expdata, np.result, 0);
    end
end
function plot_save_PhaseDependentFreeEnergyChange(np, filename, tidxList, is_save, is_show)
    for ii = 1:size(tidxList, 1)
        jjList = tidxList(ii,1):tidxList(ii,2);
        a = np.result.free_energy(:, jjList);
        if is_show
            nexttile
            hold on
            for jj = 1:size(a, 2)
                plot(1:1e4, a(:,jj))
            end
        end
        if is_save
            fprintf("%3d/%3d",ii,size(tidxList,1))
            t = np.expdata.timedata(:);
            varNames = arrayfun(@(x)sprintf("%.1fsec",x),t(jjList));
            varNames = ['Particle Size'; varNames];
            tbl = array2table([(1:1e4)',a], VariableNames=varNames);
            if size(tbl,2) > 200
                tbl = tbl(:, 1:2:end);
            end
            writetable(tbl, filename + "_free_energy.xls", Sheet = sprintf("Phase%d",ii))
            for bb = 1:7
                fprintf("\b")
            end
        end
    end
end
function plot_jn_per_area(np, stt, skip)
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
    jn = Calculate_jn(np.fitdata.xfit,np.expdata, np.fitdata.other_opt);
    colors = generate_colors_from_cmap(cmap, numel(jn));
    for ii = stt:skip:numel(jn)
        hold on
        plot(jn{ii}(:,1), jn{ii}(:,3), ...
            Color = colors(ii,:), ...
            LineWidth = 2)
        txt = sprintf("%.0f",np.expdata.timedata(ii));
        text(jn{ii}(end,1), jn{ii}(end,3), txt)
    end
    set(gca,XScale = "log")
    
    for ii = stt:skip:numel(jn)
        hold on
        idx = find(jn{ii}(:,1)>np.result.mean_theory(ii),1);
        scatter(jn{ii}(idx,1), jn{ii}(idx,3),200,'p','filled', ...
            MarkerFaceColor = [0.96,0.68,0.50])
        % pause(0.1)
        % txt = sprintf("%.0f",np.expdata.timedata(ii));
        % text(jn{ii}(idx,1)*1.05, jn{ii}(idx,3), txt)
    end
    xlabel("n")
    axis tight
    yline(0)
end