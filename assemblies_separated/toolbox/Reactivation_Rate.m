function [Pre , Post] = Reactivation_Rate(pks,TS)

Pre = [];
Post = [];
for i = 1 : size(pks,1)
    peaks = pks{i,1}(:,1);
    % Pre sleep
    l = Restrict(peaks(:,1),TS.pre);
    l = size(l,1)./sum(TS.pre(:,2)-TS.pre(:,1));
    Pre = [Pre ; l]; clear l
    
    % Post sleep
    l = Restrict(peaks(:,1),TS.post);
    l = size(l,1)./sum(TS.post(:,2)-TS.post(:,1));
    Post = [Post ; l]; clear l
end

end