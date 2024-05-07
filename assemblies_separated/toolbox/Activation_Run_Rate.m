function [Run1 , Run2] = Activation_Run_Rate(pks,TS)
% Run1 is the run where the assemblie was detected, the Run2 is the other
% one


Run1 = [];
Run2 = [];
for i = 1 : size(pks,1)
    peaks = pks{i,1}(:,1);
    % Run1
    l = Restrict(peaks(:,1),TS.run1);
    l = size(l,1)./sum(TS.run1(:,2)-TS.run1(:,1));
    Run1 = [Run1 ; l]; clear l

    % Run2
    l = Restrict(peaks(:,1),TS.run2);
    l = size(l,1)./sum(TS.run2(:,2)-TS.run2(:,1));
    Run2 = [Run2 ; l]; clear l
end

end