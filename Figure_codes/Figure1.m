function [f1c, f1e] = Figure1(PtA, PtC, Au, options)
    arguments
        PtA
        PtC
        Au
        options.save = [true, true]
        options.cmap = [];
        options.bootstrap = true
    end
    if isempty(options.cmap)
        cmap = colormap;
    else
        cmap = options.cmap;
    end
    %% Figure 1c
    f1c = figtile(1, [1, 3]);
    f1c.Position = [1, 1, 1500, 500];
    nexttile(1)
    plot(PtA.expdata.timedata, PtA.expdata.radidata,'.')
    title("Pt(acac)_2")
    
    nexttile(2)
    plot(PtC.expdata.timedata, PtC.expdata.radidata,'.')
    title("Pt(COD)Cl_2")
    
    nexttile(3)
    plot(Au.expdata.timedata, Au.expdata.radidata,'.')
    title("HAuCl_4")
    %% Figure 1e
    clc
    f1e = figtile(2, [2, 3]); f1e.Position = [1, 1, 840, 600]; % reduce size by 46%
    bootstrap = options.bootstrap;
    Plot_Mean_Var(PtA, [nexttile(1), nexttile(4)], cmap, bootstrap)
    Plot_Mean_Var(PtC, [nexttile(2), nexttile(5)], cmap, bootstrap)
    Plot_Mean_Var(Au, [nexttile(3), nexttile(6)], cmap, bootstrap)
    ax = nexttile(6);
    set(ax, YLim = [0 6.6])
    yticks(ax, 0:2:6)
    for ii = 1:6
        ax = nexttile(ii);
        set(ax, Box = "off", LineWidth = 1.1, ...
            XMinorTick = "on", YMinorTick = "on", FontSize = 16, ...
            TickLength = [0.02, 0.01], ...
            XColor = 'k', YColor = 'k', Color = "none")
        if ii <=3
            yticks(0:2:4)
            ax.YAxis.MinorTickValues = 1:2:5;
        elseif ii <= 5
            yticks(0:1:2)
            ax.YAxis.MinorTickValues = [0.5, 1.5];
            % ax.YAxis.Exponent = 1e3;
        else
            yticks(0:2:6)
            ax.YAxis.MinorTickValues = 1:2:5;
        end
        switch mod(ii, 3)
            case 1
                xticks(0:40:160)
                mtv = 20:40:180;
                xlim([0 175])
                xline([28, 41, 107, 155], 'k--', Alpha = 1, LineWidth = 1.2)
            case 2
                xticks(0:45:90)
                mtv = [45, 90]./2;
                xlim([0 90])
                xline([11, 75.5],'k--', Alpha = 1, LineWidth = 1.2)
            case 0
                xticks(0:100:200)
                mtv = [50, 150];
                xlim([0 200])
                xline(43,'k--', Alpha = 1, LineWidth = 1.2)
        end
        ax.XAxis.MinorTickValues = mtv;
        xl = xlim;
        xline(xl(end), Color = zeros(1, 3), Alpha = 1, LineWidth = 1.1)
        yl = ylim;
        yline(yl(end), Color = zeros(1, 3), Alpha = 1, LineWidth = 1.1)
        set(ax, XColor = "k", YColor = "k")
    end
    
    if options.save(1)
        exportgraphics(f1c, "Figure1c.pdf", ContentType = "vector", BackgroundColor = "none")
    end
    if options.save(2)
        exportgraphics(f1e, "Figure1e.pdf", ContentType = "vector", BackgroundColor = "none")
    end
    % ExportPrettyVectorGraphics(f, [])
    function Plot_Mean_Var(np, axesList, cmap, bootstrap)
        if bootstrap
            X = [t; flip(t)];
            bootstrap_means = np.result.bootstrap_mean;
            bootstrap_vars = np.result.bootstrap_var;

            mean_estimate = mean(bootstrap_means, 2);
            mean_error = std(bootstrap_means, [], 2);
            Y_mean = [mean_estimate - mean_error; flip(mean_estimate + mean_error)]./factor(1);
            var_estimate = mean(bootstrap_vars, 2);
            var_error = std(bootstrap_vars, [], 2);
            Y_var = [var_estimate - var_error; flip(var_estimate + var_error)]./factor(2);
            fill(axesList(1), X, Y_mean,'k', FaceAlpha = 0.3, EdgeColor = "none")
            fill(axesList(2), X, Y_var,'k', FaceAlpha = 0.3, EdgeColor = "none")
            
            % CI = [5, 95];
            % var_CI = prctile(bootstrap_vars, CI, 2);
            % mean_CI = prctile(bootstrap_means, CI, 2);
            % Y_mean_CI = [mean_CI(:,1); flip(mean_CI(:,2))]./factor(1);
            % Y_var_CI = [var_CI(:,1); flip(var_CI(:,2))]./factor(2);
            % fill(axesList(1), X, Y_mean_CI,'k', FaceAlpha = 0.5, EdgeColor = "none")
            % fill(axesList(2), X, Y_var_CI,'k', FaceAlpha = 0.5, EdgeColor = "none")
        end
        factor = [1e3, 1e6];
        t = np.expdata.timedata;
        tt = [t(:), t(:)];
        z = zeros(size(tt));
        axesList(1).NextPlot = "replace";
        plot(axesList(1), t, np.expdata.meannumberdata./factor(1),'ko', ...
            LineWidth = 1.2)
        axesList(1).NextPlot = "add";
        surf(axesList(1), tt, repmat(np.result.mean_theory./factor(1), 1, 2), z, tt, ...
            LineWidth = 3, ...
            EdgeColor = "interp")
        colormap(cmap)
        
        % ylabel(axesList(1), "mean")
        view(axesList(1), 2)
        grid(axesList(1), "off")
        ylim(axesList(1), [0 5.5])
        yticks([0, 2, 4]);
        yticklabels({"0","2","4"})
        axesList(1).XAxis.MinorTickValues = [1, 3, 5];
        set(axesList(1), Box = "off", ...
            XTickLabels = "", ...
            YMinorTick = "on", ...
            LineWidth = 1.1);
        
        axesList(2).NextPlot = "replace";
        plot(axesList(2), t, np.expdata.varnumberdata./factor(2),'ko', ...
            LineWidth = 1.2)
        axesList(2).NextPlot = "add";
        surf(axesList(2), tt, repmat(np.result.var_theory, 1, 2)./factor(2), z, tt, ...
            LineWidth = 3, ...
            EdgeColor = "interp")
        colormap(cmap)
        % ylabel(axesList(2), "var")
        view(axesList(2), 2)
        grid(axesList(2), "off")
        ylim([0 2.2])
        
        
    end
end