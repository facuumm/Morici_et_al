function [Run1 , Run2] = Activation_Run_Mean(spks , patterns , TS);
% This function calculates the mean reactivation for each assemblie during
% Pre and Post sleep


bins = spks(:,1);
spks = spks(:,2:end);

a = assembly_activity(patterns , spks');

TS.run1 = InIntervals(bins,TS.run1);
TS.run2 = InIntervals(bins,TS.run2);

Run1 = [];
Run2 = [];
for i = 1:size(a,1)
    Run1 = [Run1 ; nanmean(a(i,TS.run1))];
    Run2 = [Run2 ; nanmean(a(i,TS.run2))];
end

end