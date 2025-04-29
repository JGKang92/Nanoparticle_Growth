function Figure5(PtA, PtC, Au, PtE, Fe, Cd, options)
arguments
    PtA
    PtC
    Au
    PtE
    Fe
    Cd
    options.save = true
end
    ms = 8; lw = 2; lw1 = 1.3;
    y = [255, 222, 33]./255; 
    g = [46, 169, 69]./255;
    k = [0,0,0];
    b = [48, 100, 142]./255;
    r = [221, 19, 35]./255;
    
    %% Figure 5b (Marker)
    axesList = cell(1,5);
    markersToShow = 20;
    f5b = figtile(5002,[1, 5]);
    f5b.Position = [1,350, 1500, 305];
    t = f5b.Children;
    ax1 = axes(t);
    si = PtA.fitdata.xfit(PtA.fitdata.other_opt.indices{8});
    skip = floor(numel(si)./markersToShow);
    plot(ax1, PtA.expdata.timedata(1:skip:end), si(1:skip:end),'o', ...
        Color = r, MarkerSize = ms, LineWidth = lw)
    xlabel("Time (sec)")
    axesList{1} = ax1;

    ax2 = axes(t);
    si = PtC.fitdata.xfit(PtC.fitdata.other_opt.indices{8});
    skip = floor(numel(si)./markersToShow);
    plot(ax2, PtC.expdata.timedata(1:skip:end), si(1:skip:end),'s', ...
        Color = b, MarkerSize = ms, LineWidth = lw)
    ax2.XAxisLocation = "top";
    % ax2.YAxisLocation = "right";
    ax2.YAxis.Color = "none";
    % title("Pt nanoparticle")
    axesList{2} = ax2;

    ax = nexttile(2);
    si = Au.fitdata.xfit(Au.fitdata.other_opt.indices{8});
    skip = floor(numel(si)./markersToShow);
    plot(Au.expdata.timedata(1:skip:end), si(1:skip:end),'o', ...
        Color = r, MarkerSize = ms, LineWidth = lw)
    ylim([0.9, 2.1])
    xlabel("Time (sec)")
    % title("Au nanoparticle")
    axesList{3} = ax;
    
    ax = nexttile(3);
    si = PtE.fitdata.xfit(PtE.fitdata.other_opt.indices{8});
    skip = floor(numel(si)./markersToShow);
    plot(PtE.expdata.tval(1:skip:end), si(1:skip:end),'o', ...
        Color = r, MarkerSize = ms, LineWidth = lw)
    xlabel("Time (min)")
    % title("Pt nanoparticle")
    axesList{4} = ax;

    ax = nexttile(4);
    si = Fe.fitdata.xfit(Fe.fitdata.other_opt.indices{8});
    skip = floor(numel(si)./markersToShow);
    plot(Fe.expdata.tval(1:skip:end), si(1:skip:end),'o', ...
        Color = r, MarkerSize = ms, LineWidth = lw)
    xlabel("Time (min)")
    % title("Iron Oxide nanoparticle")
    axesList{5} = ax;

    ax = nexttile(5);
    si = Cd.fitdata.xfit(Cd.fitdata.other_opt.indices{8});
    skip = floor(numel(si)./markersToShow);
    plot(Cd.expdata.tval(1:skip:end), si(1:skip:end),'o', ...
        Color = r, MarkerSize = ms, LineWidth = lw)
    xlabel("Time (min)")
    % title("CdSe nanoparticle")
    axesList{6} = ax;
    for ii = 1:numel(axesList)
        ax = axesList{ii};
        yline(ax, [1, 4/3, 1.7, 2], '--', Color = "k", Alpha = 1, LineWidth = 2)
        ylim(ax, [0.9, 2.1])
        make_pretty_axes(ax);
        
        switch ii
            case 1
                ax.XColor = r;
            case 2
                ax.XColor = b;
            otherwise
                ax.XColor = "k";
        end
        ax.YRuler.TickLabelGapOffset = -2;
        ax.XRuler.TickLabelGapOffset = -4;
        if ii~= 1
            ax.YTickLabels = {};
        else
            ax.YTick = [1, 2];
        end
    end
    %% Figure 5 for exporting figure scale:40%
    xl = [20, 1e4]; yl1 = [-0.3 3]; yl2 = [0.9, 2.1];
    f5a = figtile(5001, [1, 5]); f5a.Position = [1, 1, 1500, 268];
    
    numToShow = 100;
    % ============================ Chemical potential part ============================
    % ______________________________  In-situ part ______________________________
    %--------------------- Chemical potential of Pt nanoparticle ---------------------
    nexttile(1)
    [exp_mun, ~, ~] = Extract_mun_In_situ_For_Figure5(PtC.fitdata.xfit, PtC.expdata, PtC.fitdata.other_opt);
    x = exp_mun(:,1);
    TF = ismember(x, round(logspace(log10(x(1)),log10(x(end)),numToShow)));
    plot(x(TF), exp_mun(TF, 2),'ks', MarkerSize = ms, LineWidth = lw1)
    hold on
    [exp_mun, muns_mn, muns_log] = Extract_mun_In_situ_For_Figure5(PtA.fitdata.xfit, PtA.expdata, PtA.fitdata.other_opt);
    x = exp_mun(:,1);
    TF = ismember(x,round(logspace(log10(x(1)), log10(x(end)),numToShow)));
    plot(x(TF), exp_mun(TF, 2),'ko', MarkerSize = ms, LineWidth = lw1)
    plot(muns_log(:, 1), muns_log(:, 2),'s-', ...
        MarkerSize = ms, LineWidth = lw, Color = g)
    plot(muns_mn(:, 1), muns_mn(:,2),'k-', LineWidth = lw)
    x = 20:1e4; factor = 9;
    plot(x, factor.*x.^(-1/3),'k--', LineWidth = 2)
    [muns_max, idx] = max(muns_mn(:,2));
    plot(muns_mn(idx,1),muns_max,'p', MarkerSize = ms*2, ...
        MarkerFaceColor = y, MarkerEdgeColor = "none")
    
    % --------------------- Chemical potential of Au nanoparticle ---------------------
    nexttile(2)
    [exp_mun, muns_mn, muns_log] = Extract_mun_In_situ_For_Figure5(Au.fitdata.xfit, Au.expdata, Au.fitdata.other_opt);
    x = exp_mun(:, 1);
    TF = ismember(x,round(logspace(log10(x(1)), log10(x(end)),numToShow)));
    plot(x(TF), exp_mun(TF, 2),'ko', MarkerSize = ms, LineWidth = lw1)
    hold on
    plot(muns_log(:, 1), muns_log(:, 2), 's-', ...
        MarkerSize = ms, LineWidth = lw, Color = g)
    plot(muns_mn(:, 1), muns_mn(:, 2), 'k-', LineWidth = lw)
    [muns_max, idx] = max(muns_mn(:,2));
    x = 20:1e4; factor = 15;
    plot(x, factor.*x.^(-1/3),'k--', LineWidth = 2)
    plot(muns_mn(idx,1),muns_max,'p', MarkerSize = ms*2, ...
        MarkerFaceColor = y, MarkerEdgeColor = "none")
    % ______________________________  In-situ part ______________________________
    
    % ______________________________  Ex-situ part ______________________________
    %--------------------- Chemical potential of Pt nanoparticle ---------------------
    nexttile(3)
    [exp_mun, muns_mn, muns_log] = Extract_mun_Ex_situ_For_Figure5(PtE.fitdata.xfit, PtE.expdata, PtE.fitdata.other_opt);
    x = exp_mun(:, 1);
    TF = ismember(x,round(logspace(log10(x(1)), log10(x(end)),numToShow)));
    plot(x(TF), exp_mun(TF, 2),'ko', MarkerSize = ms, LineWidth = lw1)
    hold on
    plot(muns_log(:, 1), muns_log(:, 2), 's-', ...
        MarkerSize = ms, LineWidth = lw, Color = g)
    plot(muns_mn(:, 1), muns_mn(:, 2), 'k-', LineWidth = lw)
    [muns_max, idx] = max(muns_mn(:,2));
    plot(muns_mn(idx,1),muns_max,'p', MarkerSize = ms*2, ...
        MarkerFaceColor = y, MarkerEdgeColor = "none")
    x = 20:1e4; factor = 10;
    plot(x, factor.*x.^(-1/3),'k--', LineWidth = 2)
    
    %--------------------- Chemical potential of Fe nanoparticle ---------------------
    nexttile(4)
    [exp_mun, muns_mn, muns_log] = Extract_mun_Ex_situ_For_Figure5(Fe.fitdata.xfit, Fe.expdata, Fe.fitdata.other_opt);
    x = exp_mun(:, 1);
    TF = ismember(x,round(logspace(log10(x(1)), log10(x(end)),numToShow)));
    plot(x(TF), exp_mun(TF, 2),'ko', MarkerSize = ms, LineWidth = lw1)
    hold on
    plot(muns_log(:, 1), muns_log(:, 2), 's-', ...
        MarkerSize = ms, LineWidth = lw, Color = g)
    plot(muns_mn(:, 1), muns_mn(:, 2), 'k-', LineWidth = lw)
    x = 20:1e4; factor = 5.5;
    plot(x, factor.*x.^(-1/3),'k--', LineWidth = 2)
    [muns_max, idx] = max(muns_mn(:,2));
    plot(muns_mn(idx,1),muns_max,'p', MarkerSize = ms*2, ...
        MarkerFaceColor = y, MarkerEdgeColor = "none")
    
    %--------------------- Chemical potential of Cd nanoparticle ---------------------
    nexttile(5)
    [exp_mun, muns_mn, muns_log] = Extract_mun_Ex_situ_For_Figure5(Cd.fitdata.xfit, Cd.expdata, Cd.fitdata.other_opt);
    x = exp_mun(:, 1);
    TF = ismember(x,round(logspace(log10(x(1)), log10(x(end)),numToShow)));
    plot(x(TF), exp_mun(TF, 2),'ko', MarkerSize = ms, LineWidth = lw1)
    hold on
    plot(muns_log(:, 1), muns_log(:, 2), 's-', ...
        MarkerSize = ms, LineWidth = lw, Color = g)
    plot(muns_mn(:, 1), muns_mn(:, 2), 'k-', LineWidth = lw)
    x = 20:1e4; factor = 3.5;
    plot(x, factor.*x.^(-1/3),'k--', LineWidth = 2)
    [muns_max, idx] = max(muns_mn(:,2));
    plot(muns_mn(idx,1),muns_max,'p', MarkerSize = ms*2, ...
        MarkerFaceColor = y, MarkerEdgeColor = "none")
    % ______________________________  Ex-situ part ______________________________
    % ============================ Chemical potential part ============================
    
    % ============================ Export graphics part ============================
    for ii = 1:numel(f5a.Children.Children)
        ax = f5a.Children.Children(ii);
        yline(ax, 0, Color = k, Alpha = 1)
        set(ax, ylim = yl1, xlim = xl, XScale = "log", ...
            XTick = 10.^[2, 3, 4], YTick = 0:1:3, XColor = k, ...
            LineWidth = 1, Box = "off")
        xline(ax, xl(end), LineWidth = 1, Color = k, Alpha=1)
        yline(ax, yl1(end), LineWidth = 1, Color = k, Alpha=1)
        ax.YRuler.TickLabelGapOffset = -2;
        ax.XRuler.TickLabelGapOffset = -4;
        set(ax, XColor = 'k', YColor = 'k', Color = "none", ...
            LineWidth = 1, FontSize = 16)
        if ii ~= 5
            set(ax, YTickLabels = {})
        end
    end
    %%
    if options.save
        exportgraphics(f5a, "Figure5a.pdf", ContentType="vector", BackgroundColor="none")
        exportgraphics(f5b, "Figure5b.pdf", ContentType="vector", BackgroundColor="none")
    end
    % ============================ Export graphics part ============================
    
    function Plot_muns_(ax, np, factor, numToShow)
        re = np.result;
        exp_mun_ = re.exp_mun;
        muns_mn_ = re.muns_mn;
        muns_log_ = re.muns_log;
        
        [muns_max_, muns_max_pos] = max(muns_mn_(:,2));
        
        m = exp_mun_(1,1);
        n_ = size(exp_mun_,1);
        semilogx(ax, exp_mun_(:,1), exp_mun_(:,2), 'ko', ...
            LineWidth = 2, ...
            MarkerSize = 15, MarkerIndices = unique(round(logspace(log10(m),log10(n_),numToShow))))
        ax.NextPlot = 'add';
        semilogx(ax, muns_log_(:, 1), muns_log_(:, 2), 'gs-', ...
            LineWidth = 2, ...
            DisplayName = 'ln(k_{n+1}^d/k_n^a\rho_{1,\infty})', ...
            MarkerSize = 15)
        semilogx(ax, muns_mn_(:, 1), muns_mn_(:, 2), 'k-', ...
            DisplayName = 'Equation (1)', ...
            LineWidth = 2)
        semilogx(ax, muns_mn_(muns_max_pos, 1), muns_max_, 'yp', ...
            MarkerFaceColor = 'y', ...
            LineWidth = 2, ...
            DisplayName = '\mu_{n^+}^s', ...
            MarkerSize = 15)
        x_ = 20:1e4;
        semilogx(x_, factor*x_.^(-1/3), 'k--', ...
            LineWidth = 2, ...
            DisplayName = '~n^{-1/3}')
        legend(Color = 'none', FontSize = 12)
        axis([20, 1e4, -0.5, 3])
    end
end