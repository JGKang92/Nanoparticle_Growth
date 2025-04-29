function [y, ind_new, lb_new, ub_new] = fill_parameters(x, ind_in, t_old, t_new, lb, ub, IPmethod, EPmethod)
    numParams = sum(cellfun(@numel,ind_in));
    if size(x,2) ~= numParams
        error("x should be a matrix whose size is (%d, ?)\n", numParams);
    end
    if isempty(IPmethod)
        IPmethod = "linear";
    end
    if isempty(EPmethod)
        EPmethod = "linear";
    end
    
    if nargout > 2
        lb_new = cell(1, numel(ind_in));
        ub_new = cell(1, numel(ind_in));
        ind_new = cell(1, numel(ind_in));

        for ii = 1:numel(ind_in)
            iin = ind_in{ii};
            if isscalar(iin)
                lb_new{ii} = lb(iin);
                ub_new{ii} = ub(iin);
                if ii == 1
                    ind_new{ii} = 1;
                else
                    ind_new{ii} = 1 + ind_new{ii-1}(end);
                end
            elseif numel(iin) == numel(t_old)
                Fl = griddedInterpolant(t_old, lb(iin), "linear","linear");
                lb_new{ii} = Fl(t_new);
                Fu = griddedInterpolant(t_old, ub(iin), "linear","linear");
                ub_new{ii} = Fu(t_new);
                if ii == 1
                    ind_new{ii} = 1:1:numel(t_new);
                else
                    iiout = 1:1:numel(t_new);
                    iiout = iiout + ind_new{ii-1}(end);
                    ind_new{ii} = iiout;
                end
            else
                error("%d-th parameter should either be scalar or vector with %d elements\n",ii, numel(t_old));
            end
        end
        lb_new = cell2mat(lb_new);
        ub_new = cell2mat(ub_new);
    end
    
    n = size(x, 1);
    y = cell(n, numel(ind_in));
    for jj = 1:n
        for ii = 1:numel(ind_in)
            iin = ind_in{ii};
            lb_ = unique(lb(iin));
            ub_ = unique(ub(iin));
            if isscalar(iin)
                y{jj,ii} = x(jj,iin);
            elseif numel(iin) == numel(t_old)
                F = griddedInterpolant(t_old, x(jj,iin), IPmethod, EPmethod);
                y_new = F(t_new);
                y_new(y_new>ub_) = ub_;
                y_new(y_new<lb_) = lb_;
                y{jj,ii} = y_new;
            else
                error("%d-th parameter should either be scalar or vector with %d elements\n",ii, numel(t_old));
            end
        end
    end
    y = cell2mat(y);
    
end