function Population = gacuGuess3(GenomeLength,FitnessFcn,options, expdata, other_opt, situ)
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
    pop = gacuGuess2_right(GenomeLength,FitnessFcn,options, expdata, other_opt, situ);
    i7 = ind{7};
    while true % recursively generate initial population whose supersaturation ratio is not out of the range [lb, ub]
        % find individuals whose supersaturation ratio is out of boundary.
        idxList = find(sum(pop(:,i7) > upperBound(i7) | pop(:, i7) < lowerBound(i7),2)~=0);
        if isempty(idxList)
            break
        end
        pop_ = gacuGuess2_right(GenomeLength,FitnessFcn,options, expdata, other_opt, situ);
        idxList2 = find(sum(pop_(:,i7) > upperBound(i7) | pop_(:, i7) < lowerBound(i7),2)==0);
        pop(idxList,:) = pop_(idxList2(1:numel(idxList)),:);
    end
    Population(initPopProvided+1:end,:) = pop(1:individualsToCreate,:);
end