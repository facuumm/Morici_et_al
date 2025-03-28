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
sm = 1; %smooth
dur = 1; %total duration
bin = 0.01; % bin size in sec

Peri.aversive.dHPC = [];    Peri.aversive.vHPC = [];
Peri.reward.dHPC = [];      Peri.reward.vHPC = [];


% Behavior
minimal_speed = 7; % minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods

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
                [Spikes , bins , Clusters] = spike_train_construction(spks_dHPC, clusters.dHPC, cellulartype, binSize, limits, events, false, true);
                SpikeTrainD = [bins' Spikes];
                clear limits events Spikes bins
            end
            
            if numberV > 3
                disp('Constructing SipkeTrains from vHPC')
                limits = [0 segments.Var1(end)/1000];
                events = [];
                [Spikes , bins , Clusters] = spike_train_construction(spks_vHPC, clusters.vHPC, cellulartype, binSize, limits, events, false, true);
                SpikeTrainV = [bins' Spikes];
                clear limits events Spikes bins
            end
            
            if exist('separated_assemblies.mat')
                disp('Loading Assemblies')
                load('separated_assemblies.mat')
                
                if isfield(patterns,'aversive')
                    if isfield(patterns.aversive,'dHPC')
                        if numberV > 0
                            disp('CCG of vSU locked to aversive dAssembly activations')
                            [P] = assembly_peaks_detection(patterns.aversive.dHPC ,SpikeTrainD ,th);
                            for i = 1 : size(P,2)
                                Times1 = Restrict(P{i},movement.aversive);
                                
                                for ii = 1 : numberV
                                    Times2 = spks(spks(:,1)==clusters.vHPC(ii),2);
                                    Times2 = Restrict(Times2,movement.aversive);
                                    if and(length(Times2) > 30 , length(Times1)>30)
                                        [p , lags] = PHIST(Times1,Times2,dur,bin,sm);
                                        Peri.aversive.vHPC = [Peri.aversive.vHPC ; p'];
                                    end
                                    clear Times2
                                end
                                clear Times1
                            end
                        end
                    end
                    
                    if isfield(patterns.aversive,'vHPC')
                        if numberD > 0
                            disp('CCG of dSU locked to aversive vAssembly activations')
                            [P] = assembly_peaks_detection(patterns.aversive.vHPC ,SpikeTrainV ,th);
                            for i = 1 : size(P,2)
                                Times1 = Restrict(P{i},movement.aversive);
                                
                                for ii = 1 : numberD
                                    Times2 = spks(spks(:,1)==clusters.dHPC(ii),2);
                                    Times2 = Restrict(Times2,movement.aversive);
                                    if and(length(Times2) > 30 , length(Times1)>30)
                                        [p , lags] = PHIST(Times1,Times2,dur,bin,sm);
                                        Peri.aversive.dHPC = [Peri.aversive.dHPC ; p'];
                                    end
                                    clear Times2
                                end
                                clear Times1
                            end
                        end
                    end                    
                end
                
                
                if isfield(patterns,'reward')
                    if isfield(patterns.reward,'dHPC')
                        if numberV > 0
                            disp('CCG of vSU locked to reward dAssembly activations')
                            [P] = assembly_peaks_detection(patterns.reward.dHPC ,SpikeTrainD ,th);
                            for i = 1 : size(P,2)
                                Times1 = Restrict(P{i},movement.reward);
                                
                                for ii = 1 : numberV
                                    Times2 = spks(spks(:,1)==clusters.vHPC(ii),2);
                                    Times2 = Restrict(Times2,movement.reward);
                                    if and(length(Times2) > 30 , length(Times1)>30)
                                        [p , lags] = PHIST(Times1,Times2,dur,bin,sm);
                                        Peri.reward.vHPC = [Peri.reward.vHPC ; p'];
                                    end
                                    clear Times2
                                end
                                clear Times1
                            end
                        end
                    end
                    
                    if isfield(patterns.reward,'vHPC')
                        if numberD > 0
                            disp('CCG of dSU locked to reward vAssembly activations')
                            [P] = assembly_peaks_detection(patterns.reward.vHPC ,SpikeTrainV ,th);
                            for i = 1 : size(P,2)
                                Times1 = Restrict(P{i},movement.reward);
                                
                                for ii = 1 : numberD
                                    Times2 = spks(spks(:,1)==clusters.dHPC(ii),2);
                                    Times2 = Restrict(Times2,movement.reward);
                                    if and(length(Times2) > 30 , length(Times1)>30)
                                        [p , lags] = PHIST(Times1,Times2,dur,bin,sm);
                                        Peri.reward.dHPC = [Peri.reward.dHPC ; p'];
                                    end
                                    clear Times2
                                end
                                clear Times1
                            end
                        end
                    end                    
                end
                
                
                
            end
        end
        disp('   ')
    end
end

figure,
m = Peri.aversive.dHPC;
[i ii] = max(m');
[i ii] = sort(ii);
imagesc(lags,[1:1:size(Peri.aversive.dHPC,1)] , m(ii,:)),hold on
xline(0)

figure, plot(lags,nanmean(m))

figure,
m = Peri.aversive.vHPC;
[i ii] = max(m');
[i ii] = sort(ii);
imagesc(lags,[1:1:size(Peri.aversive.vHPC,1)] , m(ii,:)),hold on
xline(0)

figure, plot(lags,nanmean(m))