function [average , time] = triggered_average_assemblies_response(events,patterns,cond,spks,limits)

bins = spks(:,1);
dt = bins(2)-bins(1);
limits_pos = [round(limits(1)/dt) round(limits(2)/dt)];
average = [];
time = [limits(1) : dt : limits(2)+dt];

a = assembly_activity(patterns(:,cond) , spks(:,2:end)');
a = zscore(a,1,2);


for i = 1 : size(a,1)
    tmp = [];
    for ii = 1 : size(events,1)
        [~ , index] = min(abs(events(ii)-bins));
        tmp = [tmp , a(i,index+limits_pos(1) : index+limits_pos(2))'];
    end
    average = [average , nanmean(tmp,2)];
end



end