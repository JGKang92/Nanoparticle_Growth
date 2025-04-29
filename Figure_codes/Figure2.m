function Figure2(PtA, PtC, Au, options)
arguments
    PtA
    PtC
    Au
    options.save = true
end
    %% Figure 2c-e reduce size by 32%
    clc
    % xranges = 1e3.*[repelem([0, 1],2,1); repelem([0,2],10,1); repelem([0,2.5],3,1); repelem([0,4],3,1); repelem([0, 8],2,1); repelem([0,9],5,1)]; 
    xranges = 1e3.*[repelem([0,2],10,1); repelem([0,2.5],5,1); 0,3;0,4;0,5;repelem([0,8],2,1); repelem([0,9],5,1)];
    yranges = 1e-4.*repelem([-1, 2], 25,1);
    f = Plot_psi_Jn_(PtA.expdata, PtA.result, ...
        [10:4:26, 30:2:38, 50:10:90, 110:10:150, 160:5:180] - PtA.expdata.timedata(1) + 1, ...
        xranges, ...
        yranges, ...
        [5, 5], ...
        reshape(1:25, 5, 5)', ...
        23, 2, 12);
    f.Position(2:4) = [1, 1375, 915];
    if options.save
        exportgraphics(f, "Figure2c.pdf", ContentType = "vector", BackgroundColor = "none")
    end
    %
    f = Plot_psi_Jn_(PtC.expdata, PtC.result, ...
        [2:2:10, 25:10:65, 76:2:82, 83] * 2, ...
        1e3.*[0,2.5;0,2.5;0,3;0,3;0,5;0,5;0,6;0,6;0,6;0,6;0,7;0,7;0,7;0,8;0,8], ...
        1e-4.*repelem([-1, 2], 15, 1), ...
        [5, 3], ...
        reshape(1:15,3,5)', ...
        24,2, 12);
    f.Position(2:4) = [1, 825, 915];
    if options.save
        exportgraphics(f, "Figure2d.pdf", ContentType = "vector", BackgroundColor = "none")
    end

    f = Plot_psi_Jn_(Au.expdata, Au.result, ...
        [10, 20:20:120, 160:80:320], ...
        1e3.*[0, 2; 0, 2.5; 0, 3; 0, 3; 0, 3; 0, 4; 0, 4;0, 6; 0, 6; 0, 8], ...
        1e-4.*repelem([-1, 2], 12, 1), ...
        [5, 2], ...
        reshape(1:10,2,5)', ...
        25, 2, 15);
    f.Position(2:4) = [1, 550, 915];
    if options.save
        exportgraphics(f, "Figure2e.pdf", ContentType = "vector", BackgroundColor = "none")
    end

    function f = Plot_psi_Jn_(expdata, result, ...
        fidx, xranges, yranges, rowcol, indices, figNum, skip, numToShow)
        f = figtile(figNum, rowcol);
        ms = 10;
        lw = 2;
        factor = [1e-4, 1e-3];
        p = '#e2007c';
        r = '#dd1323';
        b = '#364981';
        for idx = 1:prod(rowcol)
            ax = nexttile(indices(idx));
            fii = fidx(idx);
            psi_n_val_ = result.psi_n_val{fii};
            Jntilde = result.Jntilde{fii};
            Jntilde_pnt_guess = result.Jntilde_pnt_guess{fii};
            pnt = expdata.pntdata{fii};
            colororder({r,b})
            
            yyaxis left
            % numToShow = 15;
            skip_ = floor(numel(psi_n_val_(:,1))./numToShow);
            plot(psi_n_val_(1:skip_:end,1), ...
                 psi_n_val_(1:skip_:end,2)./factor(1), 'o', ...
                 Color = r, ...
                 MarkerSize = ms, ...
                 LineWidth = lw, ...
                 DisplayName = "Experiment")
            hold on
            plot(Jntilde(:,1), Jntilde(:,2)./factor(1), '--', ...
                Color = p, ...
                LineWidth = lw, ...
                DisplayName = "Theory 2")
            plot(Jntilde_pnt_guess(:,1), Jntilde_pnt_guess(:,2)./factor(1), 'r-', ...
                 LineWidth = lw, ...
                 DisplayName = "Theory 1")
            hold off
            if isempty(xranges)
                xlim([psi_n_val_(1,1), psi_n_val_(end,1)])
            else
                xlim(xranges(idx,:))
            end
            if isempty(yranges)
                ylim padded
            else
                ylim(yranges(idx,:)./factor(1))
            end
            
            yyaxis right
            plot(pnt(1:skip:end,1), pnt(1:skip:end,2)./factor(2), 'bs', ...
                'MarkerSize',ms, ...
                'LineWidth',lw, ...
                'DisplayName',"Experiment")
            hold on
            plot(1:max(expdata.maxnum), result.pnt_guess(fii,1:end)./factor(2), 'b-', ...
                'LineWidth',lw, ...
                'DisplayName',"Theory")
            hold off
            ylim([-0.1, 3])
            % ylim padded
            axisrange = axis;
            
            x = axisrange(1) + 0.95*diff(axisrange(1:2));
            y = axisrange(3) + 0.9*diff(axisrange(3:4));
            t = expdata.timedata(fii);
            if floor(t) ~= t
                warning("Time is not integer")
            end
            xt = xticks; xticks(linspace(xt(1),xt(end),3))
            text(x, y, sprintf("%.0fs",t), ...
                FontSize = 18, FontName = "Helvetica", ...
                HorizontalAlignment="right")
            set(ax, FontSize = 18, FontName = "Helvetica", ...
                Color = "none", XColor = 'k', ...
                Box = "off", LineWidth = 1.1)
            xl = xlim;
            xline(xl, Color = "k", LineWidth = 1.1, Alpha = 1)
            yl = ylim;
            yline(yl(end), Color = "k", LineWidth = 1.1, Alpha = 1)
    
            ax.XRuler.TickLabelGapOffset = -2;
            ax.YRuler.TickLabelGapOffset = -2;
        end
        xlabel(f.Children, "Particle size (monomer number, n)")
    end
end