function dn = get_IndicesName(fd)
    nIndices = cellfun(@numel, fd.other_opt.indices);
    indString = "";
    for ii = 1:numel(nIndices)
        if nIndices(ii) == 1 %time-independent
            if fd.lb(fd.other_opt.indices{ii}) == fd.ub(fd.other_opt.indices{ii})
                indString = indString + "X"; % not adjustable parameter
            else
                indString =  indString + "C";% adjustable parameter
            end 
        else % time-dependent
            if isequal(fd.lb(fd.other_opt.indices{ii}), fd.ub(fd.other_opt.indices{ii}))
                indString = indString + "\xi";% not adjustable parameter
            else
                indString = indString + "T";% adjustable parameter
            end
        end
    end
    dn = strjoin(indString,"");
    x0 = fd.optimopt.InitialPopulationMatrix;
    if ~isempty(x0)
        dn = sprintf("%s(%dx0)",dn,size(x0,1));
    end
end