function [clusters , numberD , numberV , spks , spks_dHPC , spks_vHPC , cellulartype] = load_SU_FM(path,criteria_type,criteria_fr,aversiveTS_run,rewardTS_run)
% This function upload and count the SU from my dataset.
% This is not a fucntion that could be easily used in other data sets.
% I created it just to reduce space in my codes.
% Facundo Morici 08/2024

spks = double([readNPY([path,'\spike_clusters.npy']) readNPY([path,'\spike_times.npy'])]);
K = tdfread([path,'\cluster_group.tsv']); % import clusters ID and groups (MUA,Noise,Good)
Kinfo = tdfread([path,'\cluster_info.tsv']); % import information of clusters
K = [K.cluster_id(K.group(:,1) == 'g') , Kinfo.ch(K.group(:,1) == 'g'),]; % defining good clusters
% Load neuronal classification
load([path,'\Cell_type_classification'])
K = [K , Cell_type_classification(:,6:8)];
group_dHPC = K(K(:,2) > 63,:);
group_vHPC = K(K(:,2) <= 63,:);

%Loop to select dorsal or ventral LFP and SU
% z=1 --> dorsal
% z=2 --> ventral
for z = 1:2
    if z == 1
        spks_dHPC = spks(ismember(spks(:,1),group_dHPC(:,1)),:); %keep spks from good clusters
    else
        spks_vHPC = spks(ismember(spks(:,1),group_vHPC(:,1)),:); %keep spks from good clusters
    end
end
clear z
spks_vHPC(:,2) = double(spks_vHPC(:,2))./20000;
spks_dHPC(:,2) = double(spks_dHPC(:,2))./20000;
spks(:,2) = double(spks(:,2))./20000;

% Selection of celltype to analyze
if criteria_type == 0 %pyr
    cellulartype = [K(:,1) , K(:,4) , K(:,2)];
elseif criteria_type == 1 % int
    cellulartype = [K(:,1) , not(K(:,4)) , K(:,2)];
elseif criteria_type == 2 % all
    cellulartype = [K(:,1) , ones(size(K,1),1) , K(:,2)];
end

%% Counting the Number f SU
numberD = 0;
clusters.all = [];
clusters.dHPC = [];
for ii=1:size(group_dHPC,1)
    cluster = group_dHPC(ii,1);
    Cellulartype = logical(cellulartype(cellulartype(:,1) == cluster,2));
    if Cellulartype
        a = length(Restrict(spks(spks(:,1)==cluster,2),aversiveTS_run)) / ((aversiveTS_run(2)-aversiveTS_run(1)));
        r = length(Restrict(spks(spks(:,1)==cluster,2),rewardTS_run)) / ((rewardTS_run(2)-rewardTS_run(1)));
        if or(a > criteria_fr ,  r > criteria_fr)
            numberD = numberD+1;
            clusters.all = [clusters.all ; cluster];
            clusters.dHPC = [clusters.dHPC ; cluster];
        end
    end
    clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5 r a
end

numberV = 0;
clusters.vHPC = [];
for ii=1:size(group_vHPC,1)
    cluster = group_vHPC(ii,1);
    Cellulartype = logical(cellulartype(cellulartype(:,1) == cluster,2));
    if Cellulartype
        a = length(Restrict(spks(spks(:,1)==cluster,2),aversiveTS_run)) / ((aversiveTS_run(2)-aversiveTS_run(1)));
        r = length(Restrict(spks(spks(:,1)==cluster,2),rewardTS_run)) / ((rewardTS_run(2)-rewardTS_run(1)));
        if or(a > criteria_fr ,  r > criteria_fr)
            numberV = numberV+1;
            clusters.all = [clusters.all ; cluster];
            clusters.vHPC = [clusters.vHPC ; cluster];
        end
    end
    clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5 r a
end

cellulartype = [K(:,1) , K(:,4) , K(:,2)];

clear freq limits
clear camara shock rightvalve leftvalve
clear ejeX ejeY dX dY dX_int dY_int

end