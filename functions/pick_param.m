function [xfit, idx] = pick_param(fitdata, threshold)
    if isempty(fitdata)
        xfit = [];
        idx = nan;
    else
        s = fitdata.scores;
        i1 = find(s(:,2) < threshold);
        if isempty(i1)
            [~, i2] = min(s(:,2));
            idx = i2;
        else
            [~, i2] = min(s(i1,1));
            idx = i1(i2);
        end
        xfit = fitdata.population(idx,:);
    end
    
end