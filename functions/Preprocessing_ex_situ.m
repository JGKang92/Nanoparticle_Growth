function expdata = Preprocessing_ex_situ(npdata, filename)
    timedata = npdata(1,:);
    nanoparticle_data = npdata(2:end,:);
    r3data = nanoparticle_data.^3;
    
    % Statistical properties
    meanr3data = mean(r3data, 1, 'omitnan');
    varr3data = var(r3data, 1, 1, 'omitnan');
    tval = timedata(1):timedata(end);
    
    r3data = r3data';
    lent_original = size(r3data, 1);
    dr = min(diff(r3data(:)));
    % Save variables
    time_script = datetime("now");
    
    expdata.r3data = r3data;
    expdata.meanr3data = meanr3data';
    expdata.varr3data = varr3data';
    expdata.timedata = timedata';
    expdata.tval = tval';
    
    expdata.time_script = time_script;
    % filename = "experiment_radius3_iron_oxide_new.mat"
    if ~isempty(filename)
        save(filename, '-struct', 'expdata')
        fprintf("saved file:%s\n", filename)
    end
end