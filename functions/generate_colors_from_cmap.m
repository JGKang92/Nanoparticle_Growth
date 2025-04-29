function colors = generate_colors_from_cmap(cmap, numColors)
    n = 1:1:size(cmap, 1); n = (n-n(1))./(n(end)-n(1));
    fr = griddedInterpolant(n, cmap(:,1));
    fg = griddedInterpolant(n, cmap(:,2));
    fb = griddedInterpolant(n, cmap(:,3));
    colors = nan(numColors,3);
    for ii = 1:numColors
        c = (ii-1)./(numColors-1);
        colors(ii,:) = [fr(c), fg(c), fb(c)];
    end
        
end