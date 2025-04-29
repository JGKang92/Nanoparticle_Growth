function Figure4(PtA, PtC, Au, options)
    arguments
        PtA
        PtC
        Au
        options.PtATime = [17, 35, 80, 140, 160];
        options.PtCTime = [5, 45, 80];
        options.AuTime = [10, 80];
        options.save = true;
        options.cmap = [];
    end
    PtATime = options.PtATime;
    PtCTime = options.PtCTime;
    AuTime = options.AuTime;
    if isempty(options.cmap)
        cmap = colormap;
    else
        cmap = options.cmap;
    end
    %% scale 28%
    clc
    height = 540; 
    nColumn = 5;
    f4a = figtile(41,[2, 5]);
    f4a.Position = [1, 1, 1500*nColumn/5, height];
    % title(f.Children, "Pt(acac)_2")
    Plot_mun_pnjn_(PtA.expdata, PtA.result, ...
        arrayfun(@(x)find(PtA.expdata.timedata == x), PtATime), ...
        nColumn, [0 , nColumn], cmap);
    
    nColumn = 3;
    f4b = figtile(42, [2, nColumn]);
    f4b.Position = [1, 1, 1500*nColumn/5, height];
    % title(f.Children, "Pt(COD)Cl_2")
    Plot_mun_pnjn_(PtC.expdata, PtC.result, ...
        arrayfun(@(x)find(PtC.expdata.timedata == x), PtCTime), ...
        nColumn, [0, nColumn], cmap);
    
    nColumn = 2;
    f4c = figtile(43, [2, nColumn]);
    f4c.Position = [1500*3/5, 1, 1500*2/5, height];
    % title(f.Children, "HAuCl_4")
    Plot_mun_pnjn_(Au.expdata      , Au.result      , ...
        arrayfun(@(x)find(Au.expdata.timedata==x),AuTime), ...
        nColumn, [0, nColumn], cmap);
    if options.save
        exportgraphics(f4a, "Figure4a.pdf", ContentType = "vector", BackgroundColor = "none")
        exportgraphics(f4b, "Figure4b.pdf", ContentType = "vector", BackgroundColor = "none")
        exportgraphics(f4c, "Figure4c.pdf", ContentType = "vector", BackgroundColor = "none")
    end
    function Plot_mun_pnjn_(expdata, result, fidx, col, ab, cmap)
        t = expdata.timedata; lent = numel(t);
        colors = generate_colors_from_cmap(cmap, lent);
        % In-situ only
        % plot mun
        for idx = 1:col
            fi = fidx(idx);
            % subplot(4,5,ab(1)+idx)
            ax = nexttile(ab(1) + idx);
            plot(20:1e4, result.mun_s(20:end,fi), ...
                LineWidth = 4, Color = colors(fi,:))
            yline(result.monomer_chemical_potential(fi), 'k-', ...
                Alpha = 1, LineWidth = 2)
            xline(expdata.meannumberdata(fi), 'k-', ...
                Alpha = 1, LineWidth = 1)
            stdvalue = sqrt(expdata.varnumberdata(fi));
            xline(expdata.meannumberdata(fi) - stdvalue, 'k--', Alpha = 1, LineWidth = 1)
            xline(expdata.meannumberdata(fi) + stdvalue, 'k--', Alpha = 1, LineWidth = 1)
            yline(0,'k-', Alpha = 1)
            axis([20, 1e4, -0.5, 2.5])
            xticks(10.^(2:4))
            yticks([0, 1, 2])
            xl = xlim; xline(xl(2), "k-", Alpha = 1, LineWidth = 1.2)
            yl = ylim; yline(yl(2), "k-", Alpha = 1, LineWidth = 1.2)
            set(gca, FontSize = 20, XScale = "log", ...
                Color = "none", XColor = "k", YColor = "k", ...
                Box = "off", LineWidth = 1.2)
            ax.XRuler.TickLabelGapOffset = -6;
            ax.YRuler.TickLabelGapOffset = -2;
            % ylabel(sprintf("Monomer chemical\npotential(k_{B} T)"))
        end
        % plot pnjn
        for idx = 1:col
            fi = fidx(idx);
            % subplot(4,5,ab(2)+idx)
            ax = nexttile(ab(2) + idx);
            a = expdata.p_cum_dev{fi};
            rhoct_f = pchip(expdata.timedata, result.rhoct);
            rhoct_1 = fnder(rhoct_f, 1);
            rhoct_1f = @(t)fnval(rhoct_1, t);
            aa = expdata.pcumval{fi};
            [~, ~, ib] = intersect(a(:,1), aa(:,1));
            aa_ = aa(ib, 2);
            numToShow = 20;
            skip = floor(size(a,1)./numToShow);
            plot(a(:,1), a(:,2) + aa_.*rhoct_1f(expdata.timedata(fi)), 'ko', ...
                MarkerIndices = 1:skip:size(a,1), MarkerSize = 10, ...11
                LineWidth = 2)
            % plot(a(1:2*(50+fi):end,1), a(1:2*(50+fi):end,2), 'ko')
            hold on
            b = result.pw_growth_rate{fi};
            plot(b(:,1), b(:,2), ...
                LineWidth = 4, ...
                Color = colors(fi,:))
            xline(expdata.meannumberdata(fi), 'k-', Alpha = 1, LineWidth = 1)
            yline(0,'k-', Alpha = 1, LineWidth = 1)
            hold off
            axis([20, 1e4, -2.5e-2, 7.5e-2])
            xticks(10.^(2:4))
            yticks([0,3,6].*1e-2)
            xl = xlim; xline(xl(2), "k-", Alpha = 1, LineWidth = 1.2)
            yl = ylim; yline(yl(2), "k-", Alpha = 1, LineWidth = 1.2)
            set(gca, FontSize = 20, XScale = "log", ...
                Color = "none", XColor = "k", YColor = "k", ...
                Box = "off", LineWidth = 1.2)
            ax.YAxis.Exponent = -2;
            ax.XRuler.TickLabelGapOffset = -6;
            ax.YRuler.TickLabelGapOffset = -2;
            % xlabel("Particle Size(n)")
            % ylabel(sprintf("Probability-weighted\ngrowth rate coefficient(s^-1)"))
        end
    end
end