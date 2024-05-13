function [Pre , Post] = Reactivation_Mean_peaks(spks , TS)
% This function calculates the mean reactivation for each assemblie during
% Pre and Post sleep


Pre = [];
Post = [];
for i = 1 : size(spks,1)
    peaks = spks{i,1};
    % Pre sleep
    pre = Restrict(peaks,TS.pre);
    l1 = nanmean(pre(:,2));
    Pre = [Pre ; l1];
    
    % Post sleep
    post = Restrict(peaks,TS.post);
    l2 = nanmean(post(:,2));
    Post = [Post ; l2 (l2-l1)/nanmean([pre(:,2);post(:,2)])];
    clear pre post l1 l2
end


end