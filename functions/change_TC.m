function paramout = change_TC(paramin, ind_in, ind_out, method)
    % T: time-dependent fitting parameter
    % C: time-independent fitting parameter
    % X: Not fitting parameter
    n_paramout = sum(cellfun(@numel, ind_out));
    n_paramset = size(paramin, 1); % number of set of parameters
    paramout = nan(n_paramset,n_paramout);
    for ii = 1:numel(ind_in)
        iin = ind_in{ii};
        iout = ind_out{ii};
        n_in = numel(iin);
        n_out = numel(iout);
        for jj = 1:n_paramset
            paramout(jj, iout) = change1(paramin(jj, iin), n_in, n_out, method); 
        end
    end
end
function xout = change1(xin, n_in, n_out, method)
    if n_in == n_out
        xout = xin;
    elseif n_in < n_out
        if n_in == 1
            xout = repelem(xin, n_in);
        else
            x = linspace(0,1,n_in);
            y = xin;
            F = griddedInterpolant(x, y, "spline");
            xq = linspace(0,1, n_out);
            xout = F(xq);
        end
    elseif n_out == 1
        if strcmp(method, "mean")
            xout = mean(xin);
        elseif strcmp(method, "median")
            xout = median(xin);
        end
    end
end