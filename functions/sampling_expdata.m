function expdata_new = sampling_expdata(expdata, stt_idx, bs)
    expdata_new = expdata;
    lent = numel(expdata.timedata);
    fnList = fieldnames(expdata);
    if isempty(stt_idx)
        stt_idx = 1;
    end
    for ii = 1:numel(fnList)
        fn_ = fnList{ii};
        old_data = expdata.(fn_);
        if isscalar(old_data)
            expdata_new.(fn_) = old_data;
        elseif size(old_data, 1) == lent
            expdata_new.(fn_) = old_data(stt_idx:bs:end,:);
        elseif size(old_data, 2) == lent
            expdata_new.(fn_) = old_data(:,stt_idx:bs:end);
        end
    end
end
% Made by soar8nalra@gmail.com