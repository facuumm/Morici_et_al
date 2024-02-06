clear
clc
close all

%% Parameters
path = {'E:\Rat103\usable';'E:\Rat126\Ephys\in-Pyr';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path

% for speed selection
minimal_speed = 7; % minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods

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
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)
        
        %load LFP
        if isfile('lfp.mat')
            load('lfp.mat')
        elseif isfile('lfp1.mat')
            load('lfp1.mat')
        end
        
        %Loading TS of the sessions
        disp('Uploading session time stamps')
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
        disp('Uploading sleep scoring')
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
        
        %% Awake
        disp('Uploading digital imputs')
        % Load digitalin.mat
        load('digitalin.mat')
        
        % periods of movment during eacj condition
        camara = ((camara(:,2)-camara(:,1))/2)+camara(:,1);
        if rewardTS_run(1) < aversiveTS_run(1)
            load('laps1.mat','posx','posy');
            [camaraR,~] = find((camara(:,1)-rewardTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
            pos = [camara(camaraR : camaraR+length(posx)-1),posx,posy];
            %interpolation of dropped frames
            ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
            ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
            pos(:,2) =dX_int; pos(:,3) =dY_int;%saving corrected pos
            behavior.pos.reward = [pos];
            behavior.speed.reward = LinearVelocity(pos,0);
            behavior.speed.reward(:,2) = smoothdata(behavior.speed.reward(:,2),'gaussian',1,'SamplePoints',behavior.speed.reward(:,1));
            behavior.quiet.reward = QuietPeriods( behavior.speed.reward , minimal_speed , minimal_speed_time);
            clear pos camaraR posx posy
            load('laps2.mat','posx','posy');
            [camaraA,~] = find((camara(:,1)-aversiveTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
            pos = [camara(camaraA : camaraA+length(posx)-1),posx,posy];
            %interpolation of dropped frames
            ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
            ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
            pos(:,2) =dX_int; pos(:,3) =dY_int;%saving corrected pos
            behavior.pos.aversive = [pos];
            behavior.speed.aversive = LinearVelocity(pos,0);
            behavior.speed.aversive(:,2) = smoothdata(behavior.speed.aversive(:,2),'gaussian',1,'SamplePoints',behavior.speed.aversive(:,1));
            behavior.quiet.aversive = QuietPeriods(behavior.speed.aversive , minimal_speed , minimal_speed_time);
            clear pos camaraR
        else
            load('laps2.mat','posx','posy');
            [camaraR,~] = find((camara(:,1)-rewardTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
            %         [camaraR2,~] = find((camara(:,1)-rewardTS_run(2)/1000)<0,1,'last'); %TimeStamp of the ending of aversive
            pos = [camara(camaraR : camaraR+length(posx)-1),posx,posy];
            %interpolation of dropped frames
            ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
            ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
            pos(:,2) =dX_int; pos(:,3) =dY_int;%saving corrected pos
            behavior.pos.reward = [pos];
            behavior.speed.reward = LinearVelocity(pos,0);
            behavior.speed.reward(:,2) = smoothdata(behavior.speed.reward(:,2),'gaussian',1,'SamplePoints',behavior.speed.reward(:,1));
            behavior.quiet.reward = QuietPeriods( behavior.speed.reward , minimal_speed , minimal_speed_time);
            clear pos camaraR posx posy
            load('laps1.mat','posx','posy');
            [camaraA,~] = find((camara(:,1)-aversiveTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
            pos = [camara(camaraA : camaraA+length(posx)-1),posx,posy];
            %interpolation of dropped frames
            ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
            ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
            pos(:,2) =dX_int; pos(:,3) =dY_int; %saving corrected pos
            behavior.pos.aversive = [pos];
            behavior.speed.aversive = LinearVelocity(pos,0);
            behavior.speed.aversive(:,2) = smoothdata(behavior.speed.aversive(:,2),'gaussian',1,'SamplePoints',behavior.speed.aversive(:,1));
            behavior.quiet.aversive = QuietPeriods(behavior.speed.aversive , minimal_speed , minimal_speed_time);
            clear pos camaraA posx posy
        end
        
%         %% Selection of channel for subsequent analyses
%         disp('Channel Selection for Time Series Analysis')
%         [spectrum1,f,s] = MTSpectrum(Restrict(vHPC1,REM));
%         [spectrum2,f,s] = MTSpectrum(Restrict(vHPC2,REM));
%         [spectrum3,f,s] = MTSpectrum(Restrict(vHPC3,REM));
%         [spectrum4,f,s] = MTSpectrum(Restrict(vHPC4,REM));
%         [spectrum5,f,s] = MTSpectrum(Restrict(vHPC5,REM));
%         
%         [m mm] = min(abs(f-6));
%         [m mmm] = min(abs(f-10));
%         
%         [m s] = max([mean(spectrum1(mm : mmm)) mean(spectrum2(mm : mmm)) mean(spectrum3(mm : mmm)) mean(spectrum4(mm : mmm)) mean(spectrum5(mm : mmm))]);
%         
%         if s == 1
%             vHPC = vHPC1;
%         elseif s == 2
%             vHPC = vHPC2;
%         elseif s == 3
%             vHPC = vHPC3;
%         elseif s == 4
%             vHPC = vHPC4;
%         elseif s == 5
%             vHPC = vHPC5;
%         end
  
        %% Spectrograms per sleep phase
        disp('Power-Spectrum calculation')
        if exist('dHPC')
            
            % Detrend Signal
            [dHPC,p] = Detrend(dHPC);
            clear p
            
            %All
            dHPC1 = Restrict(dHPC,REM);
            dHPC1(:,2) = dHPC1(:,2)./max(dHPC1(:,2));
            
            [spectrogram,f,s] = MTSpectrum(dHPC1,'frequency',1250,'range',[0 40]);
            REM_SpectrumD = [REM_SpectrumD , spectrogram'];
            clear spectrogram s
            %Baseline
            [spectrogram,f,s] = MTSpectrum(Restrict(dHPC1,REM_B),'frequency',1250,'range',[0 40]);
            REM_SpectrumDB = [REM_SpectrumDB , spectrogram'];
            clear spectrogram s
            %Reward
            [spectrogram,f,s] = MTSpectrum(Restrict(dHPC1,REM_R),'frequency',1250,'range',[0 40]);
            REM_SpectrumDR = [REM_SpectrumDR , spectrogram'];
            clear spectrogram s
            %Aversive
            [spectrogram,f,s] = MTSpectrum(Restrict(dHPC1,REM_A),'frequency',1250,'range',[0 40]);
            REM_SpectrumDA = [REM_SpectrumDA , spectrogram'];
            clear spectrogram s
            
            %NREM
            [spectrogram,f,s] = MTSpectrum(Restrict(dHPC,NREM),'frequency',1250,'range',[0 40]);
            NREM_SpectrumD = [NREM_SpectrumD , spectrogram'];
            clear spectrogram s
            %Baseline
            [spectrogram,f,s] = MTSpectrum(Restrict(dHPC,NREM_B),'frequency',1250,'range',[0 40]);
            NREM_SpectrumDB = [NREM_SpectrumDB , spectrogram'];
            clear spectrogram s
            %Reward
            [spectrogram,f,s] = MTSpectrum(Restrict(dHPC,NREM_R),'frequency',1250,'range',[0 40]);
            NREM_SpectrumDR = [NREM_SpectrumDR , spectrogram'];
            clear spectrogram s
            %Aversive
            [spectrogram,f,s] = MTSpectrum(Restrict(dHPC,NREM_A),'frequency',1250,'range',[0 40]);
            NREM_SpectrumDA = [NREM_SpectrumDA , spectrogram'];
            clear spectrogram s
        end
        
        
        if exist('vHPC')
            
            % Detrend Signal
            [vHPC1,p] = Detrend(vHPC1);
            clear p
            
            vHPC = Restrict(vHPC1,REM);
            vHPC(:,2) = vHPC(:,2)./max(vHPC(:,2));
            
            %REM
            %All
            [spectrogram,f,s] = MTSpectrum(vHPC,'frequency',1250,'range',[0 40]);
            REM_SpectrumV = [REM_SpectrumV , spectrogram'];
            clear spectrogram s
            %Baseline
            [spectrogram,f,s] = MTSpectrum(Restrict(vHPC,REM_B),'frequency',1250,'range',[0 40]);
            REM_SpectrumVB = [REM_SpectrumVB , spectrogram'];
            clear spectrogram s
            %Reward
            [spectrogram,f,s] = MTSpectrum(Restrict(vHPC,REM_R),'frequency',1250,'range',[0 40]);
            REM_SpectrumVR = [REM_SpectrumVR , spectrogram'];
            clear spectrogram s
            %Aversive
            [spectrogram,f,s] = MTSpectrum(Restrict(vHPC,REM_A),'frequency',1250,'range',[0 40]);
            REM_SpectrumVA = [REM_SpectrumVA , spectrogram'];
            clear spectrogram s
            
            %NREM
            [spectrogram,f,s] = MTSpectrum(Restrict(vHPC1,NREM),'frequency',1250,'range',[0 40]);
            NREM_SpectrumV = [NREM_SpectrumV , spectrogram'];
            clear spectrogram s
            %Baseline
            [spectrogram,f,s] = MTSpectrum(Restrict(vHPC1,NREM_B),'frequency',1250,'range',[0 40]);
            NREM_SpectrumVB = [NREM_SpectrumVB , spectrogram'];
            clear spectrogram s
            %Reward
            [spectrogram,f,s] = MTSpectrum(Restrict(vHPC1,NREM_R),'frequency',1250,'range',[0 40]);
            NREM_SpectrumVR = [NREM_SpectrumVR , spectrogram'];
            clear spectrogram s
            %Aversive
            [spectrogram,f,s] = MTSpectrum(Restrict(vHPC1,NREM_A),'frequency',1250,'range',[0 40]);
            NREM_SpectrumVA = [NREM_SpectrumVA , spectrogram'];
            clear spectrogram s
        end
        
        clear ripplesD ripplesV
        clear aversiveTS aversiveTS_run rewardTS rewardTS_run baselineTS
        clear coordinated coordinatedA coordinatedB coordinatedR
        clear coordinatedV coordinatedA_V coordinatedB_V coordinatedR_V
        clear coordinatedA_V_non_refined coordinatedB_V_non_refined coordinatedR_V_non_refined
        clear uncoordinated uncoordinatedA uncoordinatedA_V uncoordinatedB uncoordinatedB_V
        clear uncoordinatedR uncoordinatedR_V uncoordinatedV
        clear REM REM_A REM_B REM_R NREM NREM_A NREM_B NREM_R WAKE
        clear dHPC vHPC vHPC1 vHPC2 vHPC3 vHPC4 vHPC5 dHPC1
    end
    disp(['-- Finishing analysis from rat #',num2str(tt) , ' --'])
    disp('  ')
end

%
[N,EDGES] = histcounts(durations_NREM_B,30,'BinLimits',[20 350],'Normalization','probability');
plot(EDGES(2:end),Smooth(N,1),'k','LineWidth',1),hold on
[N,EDGES] = histcounts(durations_NREM_R,30,'BinLimits',[0 350],'Normalization','probability');
plot(EDGES(2:end),Smooth(N,1),'b','LineWidth',1),hold on
[N,EDGES] = histcounts(durations_NREM_A,30,'BinLimits',[0 350],'Normalization','probability');
plot(EDGES(2:end),Smooth(N,1),'r','LineWidth',1)


[m mm] = min(abs(f-6));
[m mmm] = min(abs(f-10));

tmpB = [];
tmpR = [];
tmpA = [];

for i = 1 : size(REM_SpectrumDB,2)
    tmpB = [tmpB ; max(REM_SpectrumDB(mm:mmm , i)) max(REM_SpectrumVB(mm:mmm , i))];
    tmpR = [tmpR ; max(REM_SpectrumDR(mm:mmm , i)) max(REM_SpectrumVR(mm:mmm , i))];
    tmpA = [tmpA ; max(REM_SpectrumDA(mm:mmm , i)) max(REM_SpectrumVA(mm:mmm , i))];
end


data = [tmpB(:,1) , ones(length(tmpB(:,1)),1) ; tmpR(:,1) , ones(length(tmpR(:,1)),1)*2 ; tmpA(:,1) , ones(length(tmpA(:,1)),1)*3] ; 

kruskalwallis(data(:,1) , data(:,2))
figure,
scatter(data(:,2),data(:,1),"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 4]),hold on
scatter([1 2 3],[nanmean(tmpB(:,1)) nanmean(tmpR(:,1)) nanmean(tmpA(:,1))],"filled"),xlim([0 4]),hold on

% scatter(data(:,2) , data(:,1),'filled'),xlim([0 4]),hold on
% scatter([1 2 3] , [median(tmpB(:,1)) median(tmpR(:,1)) median(tmpA(:,1))],'filled','MarkerFaceColor','k'),xlim([0 4]),hold on

data = [tmpB(:,2) , ones(length(tmpB(:,1)),1) ; tmpR(:,2) , ones(length(tmpR(:,1)),1)*2 ; tmpA(:,2) , ones(length(tmpA(:,1)),1)*3] ; 

kruskalwallis(data(:,1) , data(:,2))
figure,
scatter(data(:,2),data(:,1),"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 4]),hold on
scatter([1 2 3],[nanmean(tmpB(:,2)) nanmean(tmpR(:,2)) nanmean(tmpA(:,2))],"filled"),xlim([0 4]),hold on


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
