function Plot_mun_pnjn(expdata, result, fidx, col, ab, draw)
    % In-situ only
    if draw
        % plot mun
        for idx = 1:col
            fi = fidx(idx);
            % subplot(4,5,ab(1)+idx)
            nexttile(ab(1) + idx)
            semilogx(20:1e4, result.mun_s(20:end,fi), '-', 'LineWidth', 2)
            yline(result.monomer_chemical_potential(fi), 'k-', 'LineWidth', 2)
            xline(expdata.meannumberdata(fi), 'k-')
            stdvalue = sqrt(expdata.varnumberdata(fi));
            xline(expdata.meannumberdata(fi) - stdvalue, 'k--')
            xline(expdata.meannumberdata(fi) + stdvalue, 'k--')
            yline(0)
            axis([20, 1e4, -0.5, 2.5])
            xticks(10.^(2:4))
            yticks([0, 1, 2])
            set(gca, FontSize = 16)
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
            semilogx(a(:,1), a(:,2) + aa_.*rhoct_1f(expdata.timedata(fi)), 'ko', ...
                MarkerIndices = 1:skip:size(a,1), ...
                MarkerSize = 10, LineWidth = 2)
            % semilogx(a(1:2*(50+fi):end,1), a(1:2*(50+fi):end,2), 'ko')
            hold on
            b = result.pw_growth_rate{fi};
            semilogx(b(:,1), b(:,2), '-', ...
                LineWidth = 4)
            xline(expdata.meannumberdata(fi), 'k-')
            yline(0)
            hold off
            axis([20, 1e4, -2.5e-2, 7.5e-2])
            xticks(10.^(2:4))
            ax.YAxis.Exponent = -2;
            set(gca, FontSize = 16)
            % xlabel("Particle Size(n)")
            % ylabel(sprintf("Probability-weighted\ngrowth rate coefficient(s^-1)"))
        end
    end
end