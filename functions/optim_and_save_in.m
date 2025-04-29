function [xfit, cost, vars_to_save] = optim_and_save_in(...
    obj, nonlcon, x0, lb, ub, optimopt, other_opt, dirName, fileName)
    arguments
        obj
        nonlcon
        x0
        lb
        ub
        optimopt
        other_opt
        dirName
        fileName
    end
    track = other_opt.track;
    
    hwinfo = cpuinfo();
    stt_script = datetime('now', 'Format', 'yyMMddHHmmss');
    info = struct('cpu', hwinfo.CPUName, ... %cpu name
                'OS', feature('GetOS'), ... % OS 
                'version', version); % matlab version
    if isempty(dirName)
        dirName = [];
    elseif strcmp(dirName,"?cpu")
        dirName = info.cpu;
    elseif strcmp(dirName,"?OS")
        dirName = info.OS;
    elseif strcmp(dirName,"?version")
        dirName = info.version;
    end
    if ~isfolder(dirName)
        mkdir(dirName)    
    end
    solver = string(split(class(optimopt),"."));
    solver = solver(end);

    logFileName = sprintf("%s%s.log",fileName,stt_script);
    matFileName = sprintf("%s%s%s.mat",fileName,solver,stt_script);

    vars_to_save.obj = obj;
    vars_to_save.nonlcon = nonlcon;
    vars_to_save.x0 = x0;
    vars_to_save.lb = lb;
    vars_to_save.ub = ub;
    vars_to_save.optimopt = optimopt;
    vars_to_save.other_opt = other_opt;
    vars_to_save.dirName = dirName;
    vars_to_save.filename = matFileName;
    vars_to_save.hwinfo = hwinfo;
    vars_to_save.info = info;
    vars_to_save.stt_script = stt_script;
    vars_to_save.finished = false;
    if isfile(fullfile(dirName, matFileName))
        matFileName = "other_" + matFileName;
    end
    ff = fullfile(dirName, matFileName);
    save(ff,'-struct', 'vars_to_save');
    
    if isfield(other_opt,"seed")
        rng(other_opt.seed)
    else
        rng('default')
    end
    
    switch solver
        case "Fmincon"
            optimopt.OutputFcn = @(x, o, s) logging_fmincon_in(x, o, s, logFileName,dirName,track);
            problem = createOptimProblem("fmincon", ...
                objective = obj, ...
                nonlcon = nonlcon, ...
                x0 = x0, lb = lb, ub = ub, ...
                options = optimopt);
            startPoints = other_opt.startPoints;
            if isa(startPoints, "double") 
                if startPoints == 1
                    [xfit, cost, exitflag, output] = fmincon(problem);
                else
                    [xfit, cost, exitflag, output, solutions] = run(ms, problem, startPoints);
                    vars_to_save.solutions = solutions;
                end
            else
                [xfit, cost, exitflag, output, solutions] = run(ms, problem, startPoints);
                vars_to_save.solutions = solutions;
            end
        case "GaOptions"
            optimopt.OutputFcn = @(o, s, f)logging_ga_in(o, s, f, logFileName,dirName,track);
            if ~isempty(x0) % set nvars
                if ismatrix(x0)
                    nvars = size(x0, 2);
                elseif isvector(x0)
                    nvars = numel(x0);
                else
                    error("x0 should be matrix or vector\n");
                end
                optimopt.InitialPopulationMatrix = x0;
            elseif ~isempty(lb)
                nvars = numel(lb);
            elseif ~isempty(ub)
                nvars = numel(ub);
            else
                error("Can't find out the number of variables to optimize")
            end
            [xfit, cost, exitflag, output, population, scores] = ga(obj, nvars, [],[],[],[],lb, ub, nonlcon, optimopt);
            vars_to_save.population = population;
            vars_to_save.scores = scores;
        case "GamultiobjOptions"
            if ~isempty(x0) % set nvars
                if ismatrix(x0)
                    nvars = size(x0, 2);
                elseif isvector(x0)
                    nvars = numel(x0);
                else
                    error("x0 should be matrix or vector\n");
                end
                optimopt.InitialPopulationMatrix = x0;
            elseif ~isempty(lb)
                nvars = numel(lb);
            elseif ~isempty(ub)
                nvars = numel(ub);
            else
                error("x0 is not given correctly")
            end
            optimopt.OutputFcn = @(o, s, f) logging_gamulti_in(o, s, f, logFileName, dirName, track);
            [xfit, cost, exitflag, output, population, scores] = gamultiobj(obj, nvars,[],[],[],[],lb,ub,nonlcon,optimopt);
            vars_to_save.population = population;
            vars_to_save.scores = scores;
        otherwise
            error("This solver(%s) is not implemented in this function\n", solver)
    end
    end_script = datetime('now', Format = 'yyMMddHHmmss');
    
    fileID = fopen(fullfile(dirName, logFileName),"a");
    fprintf(fileID,"\nOptimization ended with exitflag = %d\n",exitflag);
    fprintf(fileID, "%s\n", output.message);
    fclose(fileID);
    vars_to_save.xfit = xfit;
    vars_to_save.cost = cost;
    vars_to_save.exitflag = exitflag;
    vars_to_save.output = output;
    vars_to_save.end_script = end_script;
    vars_to_save.finished = true; 
    save(ff,'-struct', 'vars_to_save');
    fprintf("saved file %s\n", matFileName)
end
% Made by soar8nalra@gmail.com