function [numberdata, keepidx] = remove_trajectories(numberdata, options)
    arguments
        numberdata
        options.nannum = 3;
        options.verbose = true;
    end
    nannum = options.nannum;
    [lent, ntrj] = size(numberdata);
    verbose = options.verbose;
    parfor ii = 1:ntrj
        data = numberdata(:,ii);
        nanidx = isnan(data);
        partition = []; % array of consecutive number of nans
        jj = 0;
        stop = true;
        for tt = 1:lent
            if nanidx(tt)
                if stop
                    jj = jj + 1;
                    partition(jj) = 1;
                    stop = false;
                else
                    partition(jj) = partition(jj) + 1;
                end
            else
                stop = true;
            end
        end
        if length(find(partition > nannum)) < 2
            keepidx(ii) = true;
        elseif verbose
            fprintf("%d-th trajectory will be removed\n",ii)
        end
    end
    numberdata = numberdata(:, keepidx);
end