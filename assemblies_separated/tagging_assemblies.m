clear
clc
% close all

%% Parameters
path = {'E:\Rat126\Ephys\in_Pyr';'E:\Rat103\usable';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path



% for SU
criteria_fr = 0; %criteria to include or not a SU into the analysis
criteria_n = [3 3]; % minimal number of neurons from each structure [vHPC dHPC]
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)

% For assemblies activation events
binSize = [0.015]; %for assemblie detection and activity strength
th = 2;

% Variables for CCG construction
sm = 2; %smooth
dur = 2; %total duration
bin = 0.015; % bin size in sec

tag.dHPC = [];
tag.vHPC = [];

%% Main loop, to iterate across sessions
for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    for t = 1 : length(subFolders)-2
        
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)
        
        %Loading TS of the sessions
        disp('Uploading session time stamps')
        load('session_organization.mat')
        
        %% Awake
        disp('Uploading behavioral data')
        load('behavioral_data.mat')
        
        %% Spikes
        % Load Units
        disp('Uploading Spiking activity')
        cd 'Spikesorting'
        [clusters , numberD , numberV , spks , spks_dHPC , spks_vHPC , cellulartype] = load_SU_FM(cd,criteria_type,criteria_fr,aversiveTS_run./1000,rewardTS_run./1000);
        
        
        %% Assemblies detection
        if or(numberD > 3 , numberV > 3)
            
            if numberD > 3
                disp('Constructing SipkeTrains from dHPC')
                limits = [0 segments.Var1(end)/1000];
                events = [];
                [Spikes , bins , Clusters] = spike_train_construction(spks_dHPC, clusters.dHPC, cellulartype, binSize, limits, events, false, false);
                SpikeTrainD = [bins' Spikes];
                clear limits events Spikes bins
            end
            
            if numberV > 3
                disp('Constructing SipkeTrains from vHPC')
                limits = [0 segments.Var1(end)/1000];
                events = [];
                [Spikes , bins , Clusters] = spike_train_construction(spks_vHPC, clusters.vHPC, cellulartype, binSize, limits, events, false, false);
                SpikeTrainV = [bins' Spikes];
                clear limits events Spikes bins
            end
            
            if exist('separated_assemblies_WakeBaseline.mat')
                disp('Loading Assemblies')
                load('separated_assemblies_WakeBaseline.mat')
                
                % Aversive
                if isfield(patterns,'dHPC')
                    disp('Tagging dorsal Assemblies')
                    [P] = assembly_peaks_detection(patterns.dHPC ,SpikeTrainD ,th);
                    for i = 1 : size(P,2)
                        % Rate for averisve
                        Times1 = Restrict(P{i},movement.aversive);
                        Rate1 = length(Times1)/sum(movement.aversive(:,2) - movement.aversive(:,1));
                        % Rate for reward
                        Times2 = Restrict(P{i},movement.reward);
                        Rate2 = length(Times2)/sum(movement.reward(:,2) - movement.reward(:,1));
                        % Surrogates construction
                        [surrogate1 , percentile1] = surrogate_assembly_activity(P{i},movement.aversive);
                        [surrogate2 , percentile2] = surrogate_assembly_activity(P{i},movement.reward);
                        % Store
                        tag.dHPC = [tag.dHPC ; tt t i Rate1>percentile1 Rate2>percentile2 Rate1 nanmean(surrogate1) Rate2 nanmean(surrogate2)];
                        % Clear variables from the loop
                        clear Times1 Times2 Rate1 Rate2 surrogate1 surrogate2 percentile1 percentile2
                    end
                end
                
                if isfield(patterns,'vHPC')
                    disp('Tagging ventral Assemblies')
                    [P] = assembly_peaks_detection(patterns.vHPC ,SpikeTrainV ,th);
                    for i = 1 : size(P,2)
                        % Rate for averisve
                        Times1 = Restrict(P{i},movement.aversive);
                        Rate1 = length(Times1)/sum(movement.aversive(:,2) - movement.aversive(:,1));
                        % Rate for reward
                        Times2 = Restrict(P{i},movement.reward);
                        Rate2 = length(Times2)/sum(movement.reward(:,2) - movement.reward(:,1));
                        % Surrogates construction
                        [surrogate1 , percentile1] = surrogate_assembly_activity(P{i},movement.aversive);
                        [surrogate2 , percentile2] = surrogate_assembly_activity(P{i},movement.reward);
                        % Store
                        tag.vHPC = [tag.vHPC ; tt t i Rate1>percentile1 Rate2>percentile2 Rate1 nanmean(surrogate1) Rate2 nanmean(surrogate2)];
                        % Clear variables from the loop
                        clear Times1 Times2 Rate1 Rate2 surrogate1 surrogate2 percentile1 percentile2
                    end
                end
            end
        end
        disp('   ')
    end
end

figure,
% dHPC
index1 = logical(tag.dHPC(:,4));       % all aversive
index2 = logical(tag.dHPC(:,5));       % all reward
index3 = and(index1 , index2);         % both
index4 = and(index1 , not(index2));    % only aversive
index5 = and(index2 , not(index1));    % only reward
index6 = and(not(index1),not(index2)); % none

% vHPC
index1 = logical(tag.vHPC(:,4));       % all aversive
index2 = logical(tag.vHPC(:,5));       % all reward
index3 = and(index1 , index2);         % both
index4 = and(index1 , not(index2));    % only aversive
index5 = and(index2 , not(index1));    % only reward
index6 = and(not(index1),not(index2)); % none



x = [ones(length(tag.dHPC),1) ; ones(length(tag.dHPC),1)*2];
y = []; 

