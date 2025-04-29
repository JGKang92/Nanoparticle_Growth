function Figure3(PtA, PtC, Au, options) % scale to 30%
arguments
    PtA
    PtC
    Au
    options.save = true;
    options.bootstrap = true;
end
    orders = 3:10;
    
    lw = 1.3; fs = 22;
    yt = 10.^(10:10:40);
    f = figtile(3, [1, 3]); f.Position = [1, 1, 1400, 420];
    
    nexttile(1)
    plot_moment(PtA, [0 180], [1e8, 1e40], 0:40:160, orders)
    
    nexttile(2)
    plot_moment(PtC, [0 84], [1e8 1e40], 0:20:80, orders)
    
    nexttile(3)
    plot_moment(Au, [0 200], [1e8 1e40], 0:50:200, orders)

    for ii = 1:3
        nexttile(ii)
        yticks(yt)
        set(gca, YScale = "log", FontSize = fs, Box = "off", Color = "none", LineWidth = lw, ...
            XColor = "k", YColor = "k", TickLength = [0.02, 0.01])
        xl = xlim;
        xline(xl(2), Color = "k", LineWidth = lw, Alpha = 1)
        yl = ylim;
        yline(yl(2), Color = "k", LineWidth = lw, Alpha = 1)
    end
    if options.save
        exportgraphics(f,"Figure3.pdf", ContentType = "vector", BackgroundColor = "none")
    end

    function plot_moment(np, xl, yl, xt, orders)
        numToShow = 30; 
        colors = [0, 0, 0; 1,1,1;
            244, 171, 131; %#f4ab83
            175, 170, 113; %##afaa71
            91, 186, 153; %##5bba99
            48, 184, 189; %##30b8bd
            95, 165, 215; %##5fa5d7
            146, 145, 195; %##9291c3
            171, 129, 180; %##ab81b4
            234, 100, 98; %##ea6462
            ]./255;
        t = np.expdata.timedata;
        mmt_theory = np.result.moments;
        mmt_experiment = np.expdata.moment_raw_value;
        hold on
        skip = ceil(numel(t)./numToShow);
        for order = orders
            plot(t, mmt_theory(:, order), ...
                LineWidth = 2, ...
                Color = colors(order,:))
            plot(t(1:skip:end), mmt_experiment(1:skip:end, order),'o', ...
                Color = colors(order,:), ...
                LineWidth = 2, ...
                MarkerSize = 10)
        end
        bootstrap = options.bootstrap;
        if bootstrap == true
            for order = orders
                fieldName = sprintf("bootstrap_%d",order);
                bootstrap_moment = np.result.(fieldName);
                moment_estimate = mean(bootstrap_moment, 2);
                moment_error = std(bootstrap_moment, [], 2);
                X = [t; flip(t)];
                Y = [moment_estimate - moment_error; flip(moment_estimate + moment_error)];
                if ~isempty(find(Y<0,1)) % Avoid negative value for log y-scale
                    Y(Y<0) = nan; 
                    Y = fillmissing(Y, 'linear');
                end

                fill(X, Y,colors(order,:), ...
                    FaceAlpha = 0.5, EdgeColor = "none")
                f_ = figtile(43218);
                nexttile
                fill([t; flip(t)], [moment_estimate - moment_error; flip(moment_estimate + moment_error)],colors(order,:), ...
                    FaceAlpha = 0.5, EdgeColor = "none")
                set(gca, YScale = "log")
                close(f_);
            end
        end
        xlim(xl)
        ylim(yl)
        xticks(xt)
    end
end