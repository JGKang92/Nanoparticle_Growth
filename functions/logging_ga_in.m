function [state, options, optchanged] = logging_ga_in(options, state, flag, logFileName, dirName, track)
    optchanged = false;
    fileID = fopen(fullfile(dirName,logFileName),"a");
    prefix = split(logFileName,".");
    prefix = prefix(1);
    parFileName = prefix + ".txt";
    trackFullFile = fullfile(dirName, parFileName);
    
    pop = state.Population;
    score = state.Score;
    popSize = size(pop, 1);
    
    switch flag
        case "init"
            fprintf("log file name: %s\n", logFileName)
            fprintf("history of parameter values: %s\n",parFileName)
            fprintf(fileID, "==================================================\n");
            fprintf(fileID, "Genetic Algorithm begins at %s\n", ...
                datetime("now", 'Format', "yyyy-MM-dd HH:mm:ss"));
            fprintf(fileID, "Gen   F-count  BestFuncValue  MeanFuncValue  datetime  LatestTimeOfTheEnd\n");
            if isfile(trackFullFile)
                warning("There already exists file: %s\n", parFileName);
                warning("File will be overwritten\n");
                delete(trackFullFile)
            end
        case {"iter","interrupt"}
            g = state.Generation;
            if g == 0
                fprintf(fileID, "Subproblem of Nonlinearly constrained problem--------------------------------------------------\n");
            end
            fprintf(fileID, "%04d  ",g);
            fprintf(fileID, "%07d  ",state.FunEval);
            fprintf(fileID, "%+13.6e  ",min(score));
            fprintf(fileID, "%+13.6e  ",mean(score,"omitnan"));
            time_now = datetime("now",Format = "HH:mm:ss");
            fprintf(fileID, "%s  ",string(time_now));
            t_toc = toc(state.StartTime)./(1 + g); % time consumed per generation
            time_end = time_now + seconds(t_toc*(options.MaxGenerations - g - 1));
            fprintf(fileID, "%s\n", string(datetime(time_end, Format = "yy/MM/dd HH:mm:ss")));
            writematrix([repelem(g, popSize)', score, pop], ...
                fullfile(dirName,"now_" + parFileName), ...
                WriteMode = "overwrite")
            if strcmp(track, "latest")
                writematrix([repelem(g, popSize)', score, pop], ...
                    trackFullFile, ...
                    WriteMode = "overwrite")
            % elseif strcmp(track,"all") % This makes large file
            %     writematrix([repelem(g, popSize)', pop], trackFullFile, ...
            %         WriteMode = "append")
            % elseif strcmp(track,"bestever") % This makes large file
            %     [~, idx] = min(score);
            %     writematrix([g, pop(idx,:)], trackFullFile, ...
            %         WriteMode = "append")
            elseif strcmp(track,"bestnow")
                [~, idx] = min(score);
                writematrix([g, score(idx,:), pop(idx,:)],trackFullFile, ...
                    WriteMode = "overwrite")
            elseif contains(track, "every")
                n = split(track,"_");
                n = str2double(n(end));
                if isnan(n)
                    error("N in the every_N should be a number\n")
                end
                if mod(g, n) == 0
                    writematrix([repelem(g, popSize)', score, pop], trackFullFile, ...
                        WriteMode = "append")
                end
            end
            % Update the fraction of mutation and crossover after 25 generations.
            % if state.Generation == 25
            %     options.CrossoverFraction = 0.8;
            %     optchanged = true;
            % end
        case "done"
            fprintf(fileID, "fitting finished at %s\n", ...
                datetime("now", 'Format', "yyyy-MM-dd HH:mm:ss"));
            fprintf(fileID, "==================================================\n");
    end
    fclose(fileID);
end
% Made by soar8nalra@gmail.com