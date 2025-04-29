function y_new = makevec(t, param)
    % t: expdata.timedata
    % param: parameter values for fitting
    if isscalar(param) % time-independent parameter
        y_new = repmat(param, size(t));
    elseif numel(t) == numel(param) % time-dependent parameter
        y_new = param;
    else
        error("wrong input\n");
    end
end
% Made by soar8nalra@gmail.com