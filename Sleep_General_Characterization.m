clear
clc
close all

%% Parameters
path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat103\usable';'E:\Rat132\recordings\in_pyr'};%List of folders from the path

%Sleep
NREM_SpectrumD = []; NREM_SpectrumV = [];
REM_SpectrumD = []; REM_SpectrumV = [];

NREM_SpectrumDB = []; NREM_SpectrumVB = [];
REM_SpectrumDB = []; REM_SpectrumVB = [];

NREM_SpectrumDR = []; NREM_SpectrumVR = [];
REM_SpectrumDR = []; REM_SpectrumVR = [];

NREM_SpectrumDA = []; NREM_SpectrumVA = [];
REM_SpectrumDA = []; REM_SpectrumVA = [];

durations_REM_B = []; durations_REM_R = []; durations_REM_A = [];
durations_NREM_B = []; durations_NREM_R = []; durations_NREM_A = [];


%% Main lo op, to iterate across sessions
for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    for t = 1 : length(subFolders)-2
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)
        
        %load LFP
        if isfile('lfp.mat')
            load('lfp.mat')
        elseif isfile('lfp1.mat')
            load('lfp1.mat')
        end
        
        % Detrend Signal
        [dHPC,p] = Detrend(dHPC);
        clear p
        [vHPC1,p] = Detrend(vHPC1);
        clear p
        
        %Loading TS of the sessions
        x = dir([cd,'\*.cat.evt']);
        segments = readtable([cd,'\',x.name],'FileType','text');
        clear x
        % TimeStamps of begening and end of the sleep and awake trials
        % Reward and Aversive trials
        aversiveTS = [];
        aversiveTS_run = [];
        rewardTS = [];
        rewardTS_run = [];
        for y = 1 : height(segments)
            % Baseline sleep session TS detection
            if y == 1
                baselineTS(1,1) = segments.Var1(y);
            elseif y ==2
                baselineTS(1,2) = segments.Var1(y);
            end
            % Aversive sleep session TS detection
            if strcmp(segments.Var2{y},'aversive')
                if strcmp(segments.Var3{y},'End')
                    aversiveTS(1,1) = segments.Var1(y+1);
                    aversiveTS(1,2) = segments.Var1(y+2);
                    aversiveTS_run(1,1) = segments.Var1(y-1);
                    aversiveTS_run(1,2) = segments.Var1(y);
                    A = y;
                end
                % Rewarded sleep session TS detection
            elseif strcmp(segments.Var2{y},'reward')
                if strcmp(segments.Var3{y},'End')
                    rewardTS(1,1) = segments.Var1(y+1);
                    rewardTS(1,2) = segments.Var1(y+2);
                    rewardTS_run(1,1) = segments.Var1(y-1);
                    rewardTS_run(1,2) = segments.Var1(y);
                    R = y;
                end
            end
        end
        clear y

        %% Sleep
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        
        REM = ToIntervals(states==5);    NREM = ToIntervals(states==3);    WAKE = ToIntervals(states==1);
        clear x states
        
        %keep only WAKE in HomeCage
%         WAKE = Restrict(WAKE, [aversiveTS ; rewardTS ; baselineTS] ./1000);
        
        % NREM events restriction according conditions
        NREM_B = NREM(NREM(:,2)<baselineTS(1,2)/1000,:);
        NREM_A = NREM(NREM(:,2)>aversiveTS(1,1)/1000 & NREM(:,2)<aversiveTS(1,2)/1000,:);
        NREM_R = NREM(NREM(:,2)>rewardTS(1,1)/1000 & NREM(:,2)<rewardTS(1,2)/1000,:);
        
        
        % REM events restriction according conditions
        REM_B = REM(REM(:,2)<baselineTS(1,2)/1000,:);
        REM_A = REM(REM(:,2)>aversiveTS(1,1)/1000 & REM(:,2)<aversiveTS(1,2)/1000,:);
        REM_R = REM(REM(:,2)>rewardTS(1,1)/1000 & REM(:,2)<rewardTS(1,2)/1000,:);
        
        %Calculation of durations acorss conditions
        durations_REM_B = [durations_REM_B ; REM_B(:,2) - REM_B(:,1)];
        durations_REM_R = [durations_REM_R ; REM_R(:,2) - REM_R(:,1)];
        durations_REM_A = [durations_REM_A ; REM_A(:,2) - REM_A(:,1)];
        durations_NREM_B = [durations_NREM_B ; NREM_B(:,2)- NREM_B(:,1)];
        durations_NREM_R = [durations_NREM_R ;  NREM_R(:,2)- NREM_R(:,1)];
        durations_NREM_A = [durations_NREM_A ;  NREM_A(:,2)- NREM_A(:,1)];
        
        %% Spectrograms per sleep phase
        %REM
        %All
        [spectrogram,f,s] = MTSpectrum(Restrict(dHPC,REM),'frequency',1250,'range',[0 40]);
        REM_SpectrumD = [REM_SpectrumD , spectrogram'];
        clear spectrogram f s
        [spectrogram,f,s] = MTSpectrum(Restrict(vHPC1,REM),'frequency',1250,'range',[0 40]);
        REM_SpectrumV = [REM_SpectrumV , spectrogram'];
        clear spectrogram f s
        %Baseline
        [spectrogram,f,s] = MTSpectrum(Restrict(dHPC,REM_B),'frequency',1250,'range',[0 40]);
        REM_SpectrumDB = [REM_SpectrumDB , spectrogram'];
        clear spectrogram f s
        [spectrogram,f,s] = MTSpectrum(Restrict(vHPC1,REM_B),'frequency',1250,'range',[0 40]);
        REM_SpectrumVB = [REM_SpectrumVB , spectrogram'];
        clear spectrogram f s
        %Reward
        [spectrogram,f,s] = MTSpectrum(Restrict(dHPC,REM_R),'frequency',1250,'range',[0 40]);
        REM_SpectrumDR = [REM_SpectrumDR , spectrogram'];
        clear spectrogram f s
        [spectrogram,f,s] = MTSpectrum(Restrict(vHPC1,REM_R),'frequency',1250,'range',[0 40]);
        REM_SpectrumVR = [REM_SpectrumVR , spectrogram'];
        clear spectrogram f s       
        %Aversive
        [spectrogram,f,s] = MTSpectrum(Restrict(dHPC,REM_A),'frequency',1250,'range',[0 40]);
        REM_SpectrumDA = [REM_SpectrumDA , spectrogram'];
        clear spectrogram f s
        [spectrogram,f,s] = MTSpectrum(Restrict(vHPC1,REM_A),'frequency',1250,'range',[0 40]);
        REM_SpectrumVA = [REM_SpectrumVA , spectrogram'];
        clear spectrogram f s              
        
        %NREM
        [spectrogram,f,s] = MTSpectrum(Restrict(dHPC,NREM),'frequency',1250,'range',[0 40]);
        NREM_SpectrumD = [NREM_SpectrumD , spectrogram'];
        clear spectrogram f s
        [spectrogram,f,s] = MTSpectrum(Restrict(vHPC1,NREM),'frequency',1250,'range',[0 40]);
        NREM_SpectrumV = [NREM_SpectrumV , spectrogram'];
        clear spectrogram s f
        %Baseline
        [spectrogram,f,s] = MTSpectrum(Restrict(dHPC,NREM_B),'frequency',1250,'range',[0 40]);
        NREM_SpectrumDB = [NREM_SpectrumDB , spectrogram'];
        clear spectrogram f s
        [spectrogram,f,s] = MTSpectrum(Restrict(vHPC1,NREM_B),'frequency',1250,'range',[0 40]);
        NREM_SpectrumVB = [NREM_SpectrumVB , spectrogram'];
        clear spectrogram f s
        %Reward
        [spectrogram,f,s] = MTSpectrum(Restrict(dHPC,NREM_R),'frequency',1250,'range',[0 40]);
        NREM_SpectrumDR = [NREM_SpectrumDR , spectrogram'];
        clear spectrogram f s
        [spectrogram,f,s] = MTSpectrum(Restrict(vHPC1,NREM_R),'frequency',1250,'range',[0 40]);
        NREM_SpectrumVR = [NREM_SpectrumVR , spectrogram'];
        clear spectrogram f s       
        %Aversive
        [spectrogram,f,s] = MTSpectrum(Restrict(dHPC,NREM_A),'frequency',1250,'range',[0 40]);
        NREM_SpectrumDA = [NREM_SpectrumDA , spectrogram'];
        clear spectrogram f s
        [spectrogram,f,s] = MTSpectrum(Restrict(vHPC1,NREM_A),'frequency',1250,'range',[0 40]);
        NREM_SpectrumVA = [NREM_SpectrumVA , spectrogram'];
        clear spectrogram s 
        
        t
        clear ripplesD ripplesV
        clear aversiveTS aversiveTS_run rewardTS rewardTS_run baselineTS
        clear coordinated coordinatedA coordinatedB coordinatedR
        clear coordinatedV coordinatedA_V coordinatedB_V coordinatedR_V
        clear coordinatedA_V_non_refined coordinatedB_V_non_refined coordinatedR_V_non_refined
        clear uncoordinated uncoordinatedA uncoordinatedA_V uncoordinatedB uncoordinatedB_V
        clear uncoordinatedR uncoordinatedR_V uncoordinatedV
        clear REM REM_A REM_B REM_R NREM NREM_A NREM_B NREM_R WAKE
    end
    tt
end

%
[N,EDGES] = histcounts(durations_NREM_B,30,'BinLimits',[20 350],'Normalization','probability');
plot(EDGES(2:end),Smooth(N,1),'k','LineWidth',1),hold on
[N,EDGES] = histcounts(durations_NREM_R,30,'BinLimits',[0 350],'Normalization','probability');
plot(EDGES(2:end),Smooth(N,1),'b','LineWidth',1),hold on
[N,EDGES] = histcounts(durations_NREM_A,30,'BinLimits',[0 350],'Normalization','probability');
plot(EDGES(2:end),Smooth(N,1),'r','LineWidth',1)

% ----- REM FIGURE -----
figure,
x = REM_SpectrumD;  
m = mean(x,2);
s = std(x,0,2)/sqrt(length(x));
subplot(241),plot(f,m,'k','LineWidth',1),hold on
ciplot(m-s , m+s , f , 'k'), alpha 0.05,ylim([0 4]),xlim([0 20])
clear x m s

x = REM_SpectrumDB;
m = mean(x,2);
s = std(x,0,2)/sqrt(length(x));
subplot(242),plot(f,m,'k','LineWidth',1),hold on
ciplot(m-s , m+s , f , 'k'), alpha 0.05,ylim([0 4]),xlim([0 20])
clear x m s

x = REM_SpectrumDR;
m = mean(x,2);
s = std(x,0,2)/sqrt(length(x));
subplot(243),plot(f,m,'b','LineWidth',1),hold on
ciplot(m-s , m+s , f , 'b'), alpha 0.05,ylim([0 4]),xlim([0 20])
clear x m s

x = REM_SpectrumDA;
m = mean(x,2);
s = std(x,0,2)/sqrt(length(x));
subplot(244),plot(f,m,'r','LineWidth',1),hold on
ciplot(m-s , m+s , f , 'r'), alpha 0.05,ylim([0 4]),xlim([0 20])
clear x m s

x = REM_SpectrumV;
m = mean(x,2);
s = std(x,0,2)/sqrt(length(x));
subplot(245),plot(f,m,'k','LineWidth',1),hold on
ciplot(m-s , m+s , f , 'k'), alpha 0.05,ylim([0 4]),xlim([0 20])
clear x m s

x = REM_SpectrumVB;
m = mean(x,2);
s = std(x,0,2)/sqrt(length(x));
subplot(246),plot(f,m,'k','LineWidth',1),hold on
ciplot(m-s , m+s , f , 'k'), alpha 0.05,ylim([0 4]),xlim([0 20])
clear x m s

x = REM_SpectrumVR;
m = mean(x,2);
s = std(x,0,2)/sqrt(length(x));
subplot(247),plot(f,m,'b','LineWidth',1),hold on
ciplot(m-s , m+s , f , 'b'), alpha 0.05,ylim([0 4]),xlim([0 20])
clear x m s

x = REM_SpectrumVA;
m = mean(x,2);
s = std(x,0,2)/sqrt(length(x));
subplot(248),plot(f,m,'r','LineWidth',1),hold on
ciplot(m-s , m+s , f , 'r'), alpha 0.05,ylim([0 4]),xlim([0 20])
clear x m s

% ----- NREM FIGURE -----
figure,
x = NREM_SpectrumD;
m = mean(x,2);
s = std(x,0,2)/sqrt(length(x));
subplot(241),plot(f,m,'k','LineWidth',1),hold on
ciplot(m-s , m+s , f , 'k'), alpha 0.05,ylim([0 4]),xlim([0 20])
clear x m s

x = NREM_SpectrumDB;
m = mean(x,2);
s = std(x,0,2)/sqrt(length(x));
subplot(242),plot(f,m,'k','LineWidth',1),hold on
ciplot(m-s , m+s , f , 'k'), alpha 0.05,ylim([0 4]),xlim([0 20])
clear x m s

x = NREM_SpectrumDR;
m = mean(x,2);
s = std(x,0,2)/sqrt(length(x));
subplot(243),plot(f,m,'b','LineWidth',1),hold on
ciplot(m-s , m+s , f , 'b'), alpha 0.05,ylim([0 4]),xlim([0 20])
clear x m s

x = NREM_SpectrumDA;
m = mean(x,2);
s = std(x,0,2)/sqrt(length(x));
subplot(244),plot(f,m,'r','LineWidth',1),hold on
ciplot(m-s , m+s , f , 'r'), alpha 0.05,ylim([0 4]),xlim([0 20])
clear x m s

x = NREM_SpectrumV;
m = mean(x,2);
s = std(x,0,2)/sqrt(length(x));
subplot(245),plot(f,m,'k','LineWidth',1),hold on
ciplot(m-s , m+s , f , 'k'), alpha 0.05,ylim([0 4]),xlim([0 20])
clear x m s

x = NREM_SpectrumVB;
m = mean(x,2);
s = std(x,0,2)/sqrt(length(x));
subplot(246),plot(f,m,'k','LineWidth',1),hold on
ciplot(m-s , m+s , f , 'k'), alpha 0.05,ylim([0 4]),xlim([0 20])
clear x m s

x = NREM_SpectrumVR;
m = mean(x,2);
s = std(x,0,2)/sqrt(length(x));
subplot(247),plot(f,m,'b','LineWidth',1),hold on
ciplot(m-s , m+s , f , 'b'), alpha 0.05,ylim([0 4]),xlim([0 20])
clear x m s

x = NREM_SpectrumVA;
m = mean(x,2);
s = std(x,0,2)/sqrt(length(x));
subplot(248),plot(f,m,'r','LineWidth',1),hold on
ciplot(m-s , m+s , f , 'r'), alpha 0.05,ylim([0 4]),xlim([0 20])
clear x m s
