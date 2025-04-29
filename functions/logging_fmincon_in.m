function stop = logging_fmincon_in(x, optimValues, state, filename, dirName, track)
    arguments
        x
        optimValues
        state
        filename
        dirName;
        track = false;
    end
    prefix = split(filename,".");
    prefix = prefix(1);
    parFileName = prefix + ".txt";
    fileID = fopen(fullfile(dirName,filename),"a");
    switch state
        case "init"
            fprintf("log file name: %s\n", filename)
            fprintf("history of parameter values: %s\n",parFileName)
            fprintf(fileID, "==================================================\n");
            fprintf(fileID, "fmincon begins at %s\n", ...
                datetime("now", 'Format', "yyyy-MM-dd HH:mm:ss"));
            fprintf(fileID, "Iter  F-count  FunctionValue  datetime\n");
        case "iter"
            fprintf(fileID, "%04d  ", optimValues.iteration);
            fprintf(fileID, "%07d  ", optimValues.funccount);
            fprintf(fileID, "%+13.6d  ", optimValues.fval);
            fprintf(fileID, "%s \n", datetime("now",'Format', "HH:mm:ss"));
            if strcmp(track, "all")
                writematrix([optimValues.iteration, x(:)'], fullfile(dirName,parFileName), ...
                    'WriteMode',"append")
            elseif strcmp(track, "latest")
                writematrix([optimValues.iteration, x(:)'], fullfile(dirName,parFileName), ...
                    'WriteMode',"overwrite")
            end
        case "done"
            fprintf(fileID, "fmincon finished at %s\n", ...
                datetime("now", 'Format', "yyyy-MM-dd HH:mm:ss"));
            fprintf(fileID, "==================================================\n");
    end
    fclose(fileID);
    stop = false;
end
% Made by soar8nalra@gmail.com