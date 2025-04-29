function Population = gacuGuess2(GenomeLength,FitnessFcn,options, expdata, other_opt, situ)
    totalPopulation = sum(options.PopulationSize);
    initPopProvided = size(options.InitialPopulation,1);
    individualsToCreate = totalPopulation - initPopProvided;
    ind = other_opt.indices;
    % Initialize Population to be created
    Population = nan(totalPopulation,GenomeLength);
    % Use initial population provided already
    if initPopProvided > 0
        Population(1:initPopProvided,:) = options.InitialPopulation;
    end
    % Create remaining population
    % problemtype is either 'unconstrained', 'boundconstraints', or
    % 'linearconstraints'. Nonlinear constrained algorithm 'ALGA' does not
    % create or use initial population of its own. It calls sub-problem
    % solvers (galincon/gaunc)
    range = options.PopInitRange;
    lowerBound = range(1,:);
    upperBound = range(2,:);
    for ii = 1:numel(ind)
        ind_ = ind{ii};
        lb = lowerBound(ind_);
        ub = upperBound(ind_);
        if isequal(lb, ub)
            Pop_ii = repmat(lb, individualsToCreate, 1);
        elseif ii == 5
            Pop_ii = repmat(log10(lb), individualsToCreate, 1) + ...
                repmat((log10(ub)-log10(lb)), individualsToCreate, 1).*rand(individualsToCreate, numel(ind_));
            Pop_ii = 10.^Pop_ii;
        elseif ii == 7 % supersaturation ratio
            continue;
        else
            Pop_ii = repmat(lb, individualsToCreate, 1) + ...
                repmat((ub-lb), individualsToCreate, 1).*rand(individualsToCreate, numel(ind_));
        end
        Population(initPopProvided+1:end,ind_) = Pop_ii;
    end
    switch situ
        case "in"
            supsat = guess_supsat_In_situ(Population, expdata, other_opt);
        case "ex2"
            supsat = guess_supsat_Ex_situ2(Population, expdata, other_opt);
        otherwise
    end
    Population(:,ind{7}) = supsat; % !!warning: This overwrites InitialPopulationMatrix
    disp("Initial Population is created\n")
end