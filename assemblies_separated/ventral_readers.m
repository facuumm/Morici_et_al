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
bin = 0.01; % bin size in sec

Peri.aversive.dHPC = [];    Peri.aversive.vHPC = [];
Peri.reward.dHPC = [];      Peri.reward.vHPC = [];
Peri.all.dHPC = [];         Peri.all.vHPC = [];

Conditional.aversive.dHPC = [];    Conditional.aversive.vHPC = [];
Conditional.reward.dHPC = [];      Conditional.reward.vHPC = [];
Conditional.all.dHPC = [];         Conditional.all.vHPC = [];

modulation.aversive.dHPC = [];     modulation.aversive.vHPC = [];
modulation.reward.dHPC = [];       modulation.reward.vHPC = [];

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
            
            if exist('separated_assemblies_alllineartrack.mat')
                disp('Loading Assemblies')
                load('separated_assemblies_alllineartrack.mat')
                
                % Aversive
                if isfield(patterns,'dHPC')
                    if numberV > 0
                        disp('CCG of vSU locked to aversive dAssembly activations')
                        [P] = assembly_peaks_detection(patterns.dHPC ,SpikeTrainD ,th);
                        for i = 1 : size(P,2)
                            Times1 = Restrict(P{i},movement.aversive);
                            [pInc pDec surp] = AssemblyModulation(movement.aversive,P{i},[0 segments.Var1(end)/1000]);
                            for ii = 1 : numberV
                                Times2 = spks(spks(:,1)==clusters.vHPC(ii),2);
                                Times2 = Restrict(Times2,movement.aversive);
                                if and(length(Times2) > 30 , length(Times1)>30)
                                    [p , lags] = PHIST(Times1,Times2,dur,bin,sm);
                                    Peri.aversive.vHPC = [Peri.aversive.vHPC ; p'];
                                    
                                    % Surrogate
                                    surrogate = different_from_surrogate(Times1,Times2,0.1,dur,bin,sm,[]);
                                    Conditional.aversive.vHPC = [Conditional.aversive.vHPC ; surrogate]; clear surrogate
                                
                                    % Save the modulation of each assemblie for posterior interpretations
                                    modulation.aversive.vHPC = [modulation.aversive.vHPC ; pInc pDec surp];

                                end
                                clear Times2
                            end
                            clear Times1 pInc pDec surp
                        end
                    end
                end
                
                if isfield(patterns,'vHPC')
                    if numberD > 0
                        disp('CCG of dSU locked to aversive vAssembly activations')
                        [P] = assembly_peaks_detection(patterns.vHPC ,SpikeTrainV ,th);
                        for i = 1 : size(P,2)
                            Times1 = Restrict(P{i},movement.aversive);
                            [pInc pDec surp] = AssemblyModulation(movement.aversive,P{i},[0 segments.Var1(end)/1000]);
                           
                            for ii = 1 : numberD
                                Times2 = spks(spks(:,1)==clusters.dHPC(ii),2);
                                Times2 = Restrict(Times2,movement.aversive);
                                if and(length(Times2) > 30 , length(Times1)>30)
                                    [p , lags] = PHIST(Times1,Times2,dur,bin,sm);
                                    Peri.aversive.dHPC = [Peri.aversive.dHPC ; p'];
                                    
                                    % Surrogate
                                    surrogate = different_from_surrogate(Times1,Times2,0.1,dur,bin,sm,[]);
                                    Conditional.aversive.dHPC = [Conditional.aversive.dHPC ; surrogate]; clear surrogate
                                
                                    % Save the modulation of each assemblie for posterior interpretations
                                    modulation.aversive.dHPC = [modulation.aversive.dHPC ; pInc pDec surp];                                
                                end
                                clear Times2
                            end
                            clear Times1 pInc pDec surp
                        end
                    end
                end
                
                % Reward
                if isfield(patterns,'dHPC')
                    if numberV > 0
                        disp('CCG of vSU locked to reward dAssembly activations')
                        [P] = assembly_peaks_detection(patterns.dHPC ,SpikeTrainD ,th);
                        for i = 1 : size(P,2)
                            Times1 = Restrict(P{i},movement.reward);
                            [pInc pDec surp] = AssemblyModulation(movement.reward,P{i},[0 segments.Var1(end)/1000]);
                            for ii = 1 : numberV
                                Times2 = spks(spks(:,1)==clusters.vHPC(ii),2);
                                Times2 = Restrict(Times2,movement.reward);
                                if and(length(Times2) > 30 , length(Times1)>30)
                                    [p , lags] = PHIST(Times1,Times2,dur,bin,sm);
                                    Peri.reward.vHPC = [Peri.reward.vHPC ; p'];
                                    
                                    % Surrogate
                                    surrogate = different_from_surrogate(Times1,Times2,0.1,dur,bin,sm,[]);
                                    Conditional.reward.vHPC = [Conditional.reward.vHPC ; surrogate]; clear surrogate
                                
                                    % Save the modulation of each assemblie for posterior interpretations
                                    modulation.reward.vHPC = [modulation.reward.vHPC ; pInc pDec surp];                                   
                                end
                                clear Times2
                            end
                            clear Times1 pInc pDec surp
                        end
                    end
                end
                
                if isfield(patterns,'vHPC')
                    if numberD > 0
                        disp('CCG of dSU locked to reward vAssembly activations')
                        [P] = assembly_peaks_detection(patterns.vHPC ,SpikeTrainV ,th);
                        for i = 1 : size(P,2)
                            Times1 = Restrict(P{i},movement.reward);
                            [pInc pDec surp] = AssemblyModulation(movement.reward,P{i},[0 segments.Var1(end)/1000]);
                            for ii = 1 : numberD
                                Times2 = spks(spks(:,1)==clusters.dHPC(ii),2);
                                Times2 = Restrict(Times2,movement.reward);
                                if and(length(Times2) > 30 , length(Times1)>30)
                                    [p , lags] = PHIST(Times1,Times2,dur,bin,sm);
                                    Peri.reward.dHPC = [Peri.reward.dHPC ; p'];
                                    
                                    % Surrogate
                                    surrogate = different_from_surrogate(Times1,Times2,0.1,dur,bin,sm,[]);
                                    Conditional.reward.dHPC = [Conditional.reward.dHPC ; surrogate]; clear surrogate
                                    
                                    % Save the modulation of each assemblie for posterior interpretations
                                    modulation.reward.dHPC = [modulation.reward.dHPC ; pInc pDec surp];                                          
                                end
                                clear Times2
                            end
                            clear Times1 pInc pDec surp
                        end
                    end
                end
                
                % All periods
                if isfield(patterns,'dHPC')
                    if numberV > 0
                        disp('CCG of vSU locked to aversive dAssembly activations')
                        [P] = assembly_peaks_detection(patterns.dHPC ,SpikeTrainD ,th);
                        for i = 1 : size(P,2)
                            Times1 = Restrict(P{i},sort([movement.aversive ; movement.reward]));
                            
                            for ii = 1 : numberV
                                Times2 = spks(spks(:,1)==clusters.vHPC(ii),2);
                                Times2 = Restrict(Times2,sort([movement.aversive ; movement.reward]));
                                if and(length(Times2) > 30 , length(Times1)>30)
                                    [p , lags] = PHIST(Times1,Times2,dur,bin,sm);
                                    Peri.all.vHPC = [Peri.all.vHPC ; p'];
                                    
                                    % Surrogate
                                    surrogate = different_from_surrogate(Times1,Times2,0.1,dur,bin,sm,[]);
                                    Conditional.all.vHPC = [Conditional.all.vHPC ; surrogate]; clear surrogate
                                end
                                clear Times2
                            end
                            clear Times1
                        end
                    end
                end
                
                if isfield(patterns,'vHPC')
                    if numberD > 0
                        disp('CCG of dSU locked to aversive vAssembly activations')
                        [P] = assembly_peaks_detection(patterns.vHPC ,SpikeTrainV ,th);
                        for i = 1 : size(P,2)
                            Times1 = Restrict(P{i},sort([movement.aversive ; movement.reward]));
                            
                            for ii = 1 : numberD
                                Times2 = spks(spks(:,1)==clusters.dHPC(ii),2);
                                Times2 = Restrict(Times2,sort([movement.aversive ; movement.reward]));
                                if and(length(Times2) > 30 , length(Times1)>30)
                                    [p , lags] = PHIST(Times1,Times2,dur,bin,sm);
                                    Peri.all.dHPC = [Peri.all.dHPC ; p'];
                                    
                                    % Surrogate
                                    surrogate = different_from_surrogate(Times1,Times2,0.1,dur,bin,sm,[]);
                                    Conditional.all.dHPC = [Conditional.all.dHPC ; surrogate]; clear surrogate
                                end
                                clear Times2
                            end
                            clear Times1
                        end
                    end
                end                
            end
        end
        disp('   ')
    end
end

%% Plots
% Aversive dHPC SU locked to vHPC assemblies
figure,
m = zscore(Peri.aversive.dHPC,1,2);
[~ , i] = min(abs(lags-(-0.1)));
[~ , ii] = min(abs(lags-0.1));
ii = nanmean(m(:,i:ii)');
[~ , ii] = sort(ii,'descend');
imagesc(lags,[1:1:size(Peri.aversive.dHPC,1)] , m(ii,:)),hold on
xline(0),clim([-3 3]), colormap 'jet'

% Aversive vHPC SU locked to sHPC assemblies
figure,
m = zscore(Peri.aversive.vHPC,1,2);
[~ , i] = min(abs(lags-(-0.05)));
[~ , ii] = min(abs(lags-0.05));
ii = nanmean(m(:,i:ii)');
[~ , ii] = sort(ii,'descend');
imagesc(lags,[1:1:size(Peri.aversive.vHPC,1)] , m(ii,:)),hold on
xline(0)
xline(0),clim([-3 3]), colormap 'jet'

% Reward dHPC SU locked to vHPC assemblies
figure,
m = zscore(Peri.reward.dHPC,1,2);
[~ , i] = min(abs(lags-(-0.1)));
[~ , ii] = min(abs(lags-0.1));
ii = nanmean(m(:,i:ii)');
[~ , ii] = sort(ii,'descend');
imagesc(lags,[1:1:size(Peri.reward.dHPC,1)] , m(ii,:)),hold on
xline(0),clim([-3 3]), colormap 'jet'

% Reward vHPC SU locked to dHPC assemblies
figure,
m = zscore(Peri.reward.vHPC,1,2);
[~ , i] = min(abs(lags-(-0.1)));
[~ , ii] = min(abs(lags-0.1));
ii = nanmean(m(:,i:ii)');
[~ , ii] = sort(ii,'descend');
imagesc(lags,[1:1:size(Peri.reward.vHPC,1)] , m(ii,:)),hold on
xline(0),clim([-3 3]), colormap 'jet'

%% Selection of responsive pairs
% Aversive dHPC SU locked to vHPC assemblies
m = zscore(Peri.all.dHPC,1,2);

i   = logical(Conditional.all.dHPC(:,1));
ii  = logical(Conditional.all.dHPC(:,2));
iii = logical(and(not(Conditional.all.dHPC(:,1)) , not(Conditional.all.dHPC(:,2))));

figure,
plot(lags,nanmean(m(i,:)),'b'),hold on
ciplot(nanmean(m(i,:))-nansem(m(i,:)) , nanmean(m(i,:))+nansem(m(i,:)), lags , 'b'), alpha 0.5

plot(lags,nanmean(m(ii,:)),'r'),hold on
ciplot(nanmean(m(ii,:))-nansem(m(ii,:)) , nanmean(m(ii,:))+nansem(m(ii,:)), lags , 'r'), alpha 0.5

plot(lags,nanmean(m(iii,:)),'k'),hold on
ciplot(nanmean(m(iii,:))-nansem(m(iii,:)) , nanmean(m(iii,:))+nansem(m(iii,:)), lags , 'k'), alpha 0.5
xlim([-0.5 0.5])

% Aversive vHPC SU locked to vHPC assemblies
m = zscore(Peri.all.vHPC,1,2);

i   = logical(Conditional.all.vHPC(:,1));
ii  = logical(Conditional.all.vHPC(:,2));
iii = logical(and(not(Conditional.all.vHPC(:,1)) , not(Conditional.all.vHPC(:,2))));

figure,
plot(lags,nanmean(m(i,:)),'b'),hold on
ciplot(nanmean(m(i,:))-nansem(m(i,:)) , nanmean(m(i,:))+nansem(m(i,:)), lags , 'b'), alpha 0.5

plot(lags,nanmean(m(ii,:)),'r'),hold on
ciplot(nanmean(m(ii,:))-nansem(m(ii,:)) , nanmean(m(ii,:))+nansem(m(ii,:)), lags , 'r'), alpha 0.5

plot(lags,nanmean(m(iii,:)),'k'),hold on
ciplot(nanmean(m(iii,:))-nansem(m(iii,:)) , nanmean(m(iii,:))+nansem(m(iii,:)), lags , 'k'), alpha 0.5
xlim([-0.5 0.5])
