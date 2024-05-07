function pks = assemblies_peaks(spks,patterns,th)
% This funciton calculate assemblies activity and detect peaks.

bins = spks(:,1);
dt = bins(2)-bins(1); % delta time
spks = spks(:,2:end);

a = assembly_activity(patterns , spks');

a = zscore(a,1,2);
R = [];
pks = cell(size(a,1),1);
for i = 1:size(a,1)
    [p,l] = findpeaks(a(i,:),bins,'MinPeakHeight',th);
    pks{i,1} = [l,p']; clear l p
end

end