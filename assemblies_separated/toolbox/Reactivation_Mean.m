function [Pre , Post] = Reactivation_Mean(spks , patterns , TS);
% This function calculates the mean reactivation for each assemblie during
% Pre and Post sleep


bins = spks(:,1);
spks = spks(:,2:end);

a = assembly_activity(patterns , spks');

TS.pre = InIntervals(bins,TS.pre);
TS.post = InIntervals(bins,TS.post);

Pre = [];
Post = [];
for i = 1:size(a,1)
    Pre = [Pre ; nanmean(a(i,TS.pre))];
    Post = [Post ; nanmean(a(i,TS.post))];
end

end