function cost = lossReg(prediction, experiment, loss_method, scale_method)
    
    switch scale_method
        case "none"
            scale_factor = 1;
        case "var"
            scale_factor = var(experiment,'omitnan');
        case "std"
            scale_factor = std(experiment);
        case "minmax"
            scale_factor = max(experiment) - min(experiment);
        case "mad" % median absolute deviation
            scale_factor = median(abs(experiment - median(experiment,'omitnan')),'omitnan');
        otherwise
            error("method of scale is not given correctly\n");
    end
    if scale_factor == 0
        scale_factor = 1;
    end
    prediction = prediction./scale_factor;
    experiment = experiment./scale_factor;
    switch loss_method
        case "mse" % mean squared error
            cost = mean((prediction - experiment).^2);
        case "logmse" %log(mean squared error)
            cost = log(mean((prediction - experiment).^2));
        case "sse" % sum of squared error
            cost = sum((prediction - experiment).^2);
        case "logsse" % log(sum of squared error)
            cost = log(sum((prediction - experiment).^2));
        case "meanlogcosh" % mean of log cosh error
            cost = mean(log(cosh(prediction - experiment)));
        case "sumlogcosh" % sum of log cosh error
            cost = sum(log(cosh(prediction - experiment)));
        case "rmse" % relative mean squared error
            % experiment should always have positive value
            % do not scale for this loss function
            cost = mean((prediction - experiment).^2./experiment.^2);
        otherwise
            error("method of loss function is not given correctly\n");
    end
end
% Made by soar8nalra@gmail.com