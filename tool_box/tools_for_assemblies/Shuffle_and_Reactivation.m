function results = Shuffle_and_Reactivation(data, X, Y, is, patterns, cond, th , Spks, type, config, normalization)
% Its downsample data to X elements Y times to calculate Reactivation
% data: the original distribution (array or vector)
% X: the number of elements to pick in each iteration
% Y: the number of times to perform the operation
% Morici Juan Facunddo 06/02/2025


% Initialize an array to store results
tmp = [];
for i = 1:Y
    % Shuffle the data
    s = randperm(size(data,1));
    shuffled_data = data(s,:);
    
    % Pick X elements
    selected_data = shuffled_data(1:X);
    selected_data = sort([shuffled_data(:,1) shuffled_data(:,2)]);
    
    is.baseline = InIntervals(Spks(:,1),Restrict(selected_data,is.timestamps.sleep.baseline));
    is.reward = InIntervals(Spks(:,1),Restrict(selected_data,is.timestamps.sleep.reward));
    is.aversive = InIntervals(Spks(:,1),Restrict(selected_data,is.timestamps.sleep.aversive));
    
    [R] = reactivation_strength(patterns , cond , Spks , is , th , type , config , normalization , []);
    tmp(:,:,i) = R;
    clear R r
end


R = [];
for i = 1 : size(tmp,1)
    r = [];
    for ii = 1 : size(tmp,2)
        r = [r , nanmean(squeeze(tmp(i,ii,:)))];
    end
    R = [R ; r]; clear r
end

results = R;

end