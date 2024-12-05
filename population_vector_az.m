  %% Script to calculate population vector correlation of 2D ratemaps 
% Sessions need to have at least 5 neurons to be included in the pv corr
% analysis
% Silva 2023

%%
clear
clc
close all
%% Parameters
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'};%List of folders from the path
% path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat103\usable';'E:\Rat132\recordings\in_pyr'};%List of folders from the path
n_session = 5; % number of neurons per session to consider a session as valid   

%% Accumulated population vector correlation 
% Output 
pv_corr.dhpc = []; 
pv_corr.vhpc = []; 

mean_pv_dhpc = [];
mean_pv_vhpc = [];

for tt = 1:length(path)
    
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    
    temp_dhpc= [];
    temp_vhpc= [];
    
    for t = 1 : length(subFolders)-2
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)
        cd 'Spikesorting'
        disp(['-- Loding dHPC --'])
        
        try load('dHPC_pc_skaggs_circular.mat');
            if size(dHPC,2)>=n_session 
                %Create rate map stacks 
                for n=1:size(dHPC,2)
                    pv_ave_dhpc(:,:,n) = dHPC{n}.frMap_ave; 
                    pv_rew_dhpc(:,:,n) = dHPC{n}.frMap_rew;
                end
%                 if doplot==1
%                     figure(1)
%                     subplot(2,2,1)
%                     imagesc(squeeze(pv_ave_dhpc)')
%                     ylabel('Neurons')
%                     subplot(2,2,2)
%                     imagesc(squeeze(pv_rew_dhpc)')
%                     xlabel('Space')
%                     
%                 end
                %%
                % Vector bin correlation - dorsal
                sesion_corr = NaN(60,1);
                for c=1:size(pv_ave_dhpc,2) %Iterate through all the bins
                    %Create a bin vector (with the data of all neurons in that bin)
                    pv_ave = squeeze(pv_ave_dhpc(1,c,:));
            
                    pv_rew = squeeze(pv_rew_dhpc(1,c,:));
            
                    %Correlation between ave vector bin and rew vector bin
                    cor = corrcoef(pv_ave, pv_rew,'Rows','pairwise');
                    sesion_corr(c,1) = cor(1,2);    
                end 
                temp_dhpc = [temp_dhpc;sesion_corr];
                mean_pv_dhpc = [mean_pv_dhpc;nanmean(sesion_corr)];

                
%                 clear pv_ave_dhpc pv_rew_dhpc pv_ave pv_rew
                
                
            end
        catch 
            disp(['No dHPCfile in ',session]);
        end     
 
        disp(['-- Loding vHPC --'])
        
        try load('vHPC_pc_skaggs_circular.mat');
            if size(vHPC,2)>=5
                %Create rate map stacks 
                for n=1:size(vHPC,2)
                    pv_ave_vhpc(:,:,n) = vHPC{n}.frMap_ave; 
                    pv_rew_vhpc(:,:,n) = vHPC{n}.frMap_rew;
                end 
                
                % Vector bin correlation - ventral
                sesion_corr = NaN(60,1);
                for c=1:size(pv_ave_vhpc,2) %Iterate through all the bins
                    %Create a bin vector (with the data of all neurons in that bin)
                    pv_ave = squeeze(pv_ave_vhpc(1,c,:));
            
                    pv_rew = squeeze(pv_rew_vhpc(1,c,:));
            
                    %Correlation between ave vector bin and rew vector bin
                    cor = corrcoef(pv_ave, pv_rew,'Rows','pairwise');
                    sesion_corr(c,1) = cor(1,2);    
                end 
                temp_vhpc = [temp_vhpc;sesion_corr];  
                mean_pv_vhpc = [mean_pv_vhpc;nanmean(sesion_corr)];
                clear pv_ave_vhpc pv_rew_vhpc pv_ave pv_rew
            end
        catch 
            disp(['No vHPC file in ',session]);
         
        end 
              
    end
    
    pv_corr.dhpc = [pv_corr.dhpc;temp_dhpc];
    pv_corr.vhpc = [pv_corr.vhpc;temp_vhpc];
    disp(['-------- Finished rat#' , num2str(tt) , ' --------'])
    disp(' ')
end %iterate rats 

% Plots accumulated: 
figure(1);clf;hold on;subplot(1,2,1);  
sgtitle('Population vector: green-dHPC blue-vHPC')
h = cdfplot(pv_corr.dhpc);hold on;
h.Color = [0 1 0.5];
h.LineWidth = 2;
h=cdfplot(pv_corr.vhpc);hold on;
h.Color = [0.1 0.2 1];
h.LineWidth = 2;
ylabel('Cumulative freq');
xlabel('Spatial bin corr. coef.(r)')
title('Mean')

%Plots mean
c = [.3,.3,.3];
figure(1);clf;hold on;sgtitle('Mean population vector')
x=[ones(size(mean_pv_dhpc,1),1);2*ones(size(mean_pv_vhpc,1),1)];y=[mean_pv_dhpc;mean_pv_vhpc];
scatter(x,y,30,c,"filled",'jitter','on', 'jitterAmount',0.3,'MarkerFaceAlpha',.5); xlim([0 3]); 
ylabel('Mean pv bin corr');xticks([1 2]);xticklabels({'dhpc', 'vhpc'});hold on;
s1=scatter(1,nanmedian(mean_pv_dhpc), "filled");s1.MarkerFaceColor = [0 0 0];
s2=scatter(2,nanmedian(mean_pv_vhpc), "filled");s2.MarkerFaceColor = [1 0.1 0.2];

sgtitle('Mean population vector');boxplot(y,x);
ylabel('Mean pv bin corr');xticks([1 2]);xticklabels({'dhpc', 'vhpc'});
p = ranksum(mean_pv_dhpc,mean_pv_vhpc)

%% Population vector per session 
subplot(1,2,2); 
for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    
    temp_dhpc= [];
    temp_vhpc= [];
    
    for t = 1 : length(subFolders)-2
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)
        cd 'Spikesorting'
        
       try load('dHPC_pc_skaggs_circular.mat');
        
            if size(dHPC,2)>=n_session 
                %Create rate map stacks 
                for n=1:size(dHPC,2)
                    pv_ave_dhpc(:,:,n) = dHPC{n}.frMap_ave; 
                    pv_rew_dhpc(:,:,n) = dHPC{n}.frMap_rew;
                end
             
                % Vector bin correlation - dorsal
                sesion_corr = NaN(60,1);
                for c=1:size(pv_ave_dhpc,2) %Iterate through all the bins
                    %Create a bin vector (with the data of all neurons in that bin)
                    pv_ave = squeeze(pv_ave_dhpc(1,c,:));
            
                    pv_rew = squeeze(pv_rew_dhpc(1,c,:));
            
                    %Correlation between ave vector bin and rew vector bin
                    cor = corrcoef(pv_ave, pv_rew);
                    sesion_corr(c,1) = cor(1,2);  
                end 
                h = cdfplot(sesion_corr);hold on;
                h.Color = [0 1 0.5];
                h.LineWidth = 1;
                clear pv_ave_dhpc pv_rew_dhpc pv_ave pv_rew
                
            end
       catch
           disp('No dhpc file')
        end 
        
      try load('vHPC_pc_skaggs_circular.mat');
            if size(vHPC,2)>=n_session 
                %Create rate map stacks 
                for n=1:size(vHPC,2)
                    pv_ave_vhpc(:,:,n) = vHPC{n}.frMap_ave; 
                    pv_rew_vhpc(:,:,n) = vHPC{n}.frMap_rew;
                end 
                
                % Vector bin correlation - dorsal
                sesion_corr = NaN(60,1);
                for c=1:size(pv_ave_vhpc,2) %Iterate through all the bins
                    %Create a bin vector (with the data of all neurons in that bin)
                    pv_ave = squeeze(pv_ave_vhpc(1,c,:));
            
                    pv_rew = squeeze(pv_rew_vhpc(1,c,:));
            
                    %Correlation between ave vector bin and rew vector bin
                    cor = corrcoef(pv_ave, pv_rew);
                    sesion_corr(c,1) = cor(1,2);    
                    
                end 
                h = cdfplot(sesion_corr);hold on;
                h.Color = [0.1 0.2 1];
                h.LineWidth = 1;
                clear pv_ave_vhpc pv_rew_vhpc pv_ave pv_rew
            end
      catch
          disp('No vhpc file')
      end 
              
    end
    
    disp(['-------- Finished rat#' , num2str(tt) , ' --------'])
    disp(' ')
end %iterate rats 
ylabel('Cumulative freq');
xlabel('Spatial bin corr. coef.(r)')
title('Per session')

%% Accumulated population vector correlation 
% keeping same amount of neurons in each session
count=0;
for tt = 1:length(path)
    
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    
    temp_dhpc= [];
    temp_vhpc= [];
    
    for t = 1 : length(subFolders)-2
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)
        cd 'Spikesorting'
        disp(['-- Loding dHPC --'])
        
        try 
            load('dHPC_pc_skaggs_circular.mat');
            load('vHPC_pc_skaggs_circular.mat');
            if size(dHPC,2)>=size(vHPC,2)
               total_n=size(vHPC,2); 
            else 
               total_n=size(dHPC,2); 
            end 
            
            if total_n>=n_session
                
                %%%%%%%% dHPC %%%%%%%%%%
                %Create rate map stacks 
                for n=1:total_n
                    pv_ave_dhpc(:,:,n) = dHPC{n}.frMap_ave; 
                    pv_rew_dhpc(:,:,n) = dHPC{n}.frMap_rew;
                end
                % Vector bin correlation - dorsal
                sesion_corr = NaN(60,1);
                for c=1:size(pv_ave_dhpc,2) %Iterate through all the bins
                    %Create a bin vector (with the data of all neurons in that bin)
                    pv_ave = squeeze(pv_ave_dhpc(1,c,:));
            
                    pv_rew = squeeze(pv_rew_dhpc(1,c,:));
            
                    %Correlation between ave vector bin and rew vector bin
                    cor = corrcoef(pv_ave, pv_rew,'Rows','pairwise');
                    sesion_corr(c,1) = cor(1,2);    
                end
                
                temp_dhpc = [temp_dhpc;sesion_corr];
                
                %%%%%%%% vHPC %%%%%%%%%%
                %Create rate map stacks 
                for n=1:total_n
                    pv_ave_vhpc(:,:,n) = vHPC{n}.frMap_ave; 
                    pv_rew_vhpc(:,:,n) = vHPC{n}.frMap_rew;
                end 
                
                % Vector bin correlation - ventral
                sesion_corr = NaN(60,1);
                for c=1:size(pv_ave_vhpc,2) %Iterate through all the bins
                    %Create a bin vector (with the data of all neurons in that bin)
                    pv_ave = squeeze(pv_ave_vhpc(1,c,:));
            
                    pv_rew = squeeze(pv_rew_vhpc(1,c,:));
            
                    %Correlation between ave vector bin and rew vector bin
                    cor = corrcoef(pv_ave, pv_rew,'Rows','pairwise');
                    sesion_corr(c,1) = cor(1,2);    
                end 
                temp_vhpc = [temp_vhpc;sesion_corr];  
            
            end 
        catch 
            count=count+1; disp(['Not enough neurons. Discarded sessions ', num2str(c)])
        end 

              
    end
    
    pv_corr.dhpc = [pv_corr.dhpc;temp_dhpc];
    pv_corr.vhpc = [pv_corr.vhpc;temp_vhpc];
    disp(['-------- Finished rat#' , num2str(tt) , ' --------'])
    disp(' ')
end %iterate rats 

% Plots: 
figure(1);clf;

h = cdfplot(pv_corr.dhpc);hold on;
h.Color = [0 1 0.5];
h.LineWidth = 2;
h=cdfplot(pv_corr.vhpc);hold on;
h.Color = [0.1 0.2 1];
h.LineWidth = 2;
ylabel('Cumulative freq');
xlabel('Spatial bin corr. coef.(r)')
title('Mean population vector')
legend('dHPC','vHPC')


%% PPOPULATION VECTOR BETWEEN-WITHIN
clear
clc
close all
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
    '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path

% for SU
criteria_fr = 0.01; %criteria to include or not a SU into the analysis
criteria_n = 6; % minimal number of neurons from each structure
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
n_SU_V = [];
n_SU_D = [];
Xedges = 60; %number of bins for RateMap construction - 2.5cm bin (because of the extremes cutting)
sigma = 2;%round(15/(180/Xedges)); %defined for gauss kernel of 15cm
binSize = 0.001; % bin size for replay events detection
bin_size = 1; % to bin pos ans spks in between/within  
% Behavior
minimal_speed = 2.5;% minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods
n_session = 5;


temp_dhpc_bet = [];
temp_dhpc_within_ave = [];
temp_dhpc_within_rew = [];

temp_vhpc_bet = [];
temp_vhpc_within_ave = [];
temp_vhpc_within_rew = [];

pv_corr.dhpc_bet = [];
pv_corr.dhpc_within_ave = [];
pv_corr.dhpc_within_rew = [];

pv_corr.vhpc_bet = [];
pv_corr.vhpc_within_ave = [];
pv_corr.vhpc_within_rew = [];

%% MAIN LOOP ORIGINAL, to iterate across sessions  - within/between by bins 
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
        
        %Session output
        dHPC = {};% one cell per pc
        vHPC = {}; 
        
        %% Awake
        disp('Uploading digital imputs')
        % Load digitalin.mat
        load('digitalin.mat')
        %Shocks selection
        Shocks_filt = Restrict(shock,aversiveTS_run ./1000);
        % Keep only the first shock of each TTL (first from 20)
        count = 1;
        deff = [];
        for i = 1:length(Shocks_filt)
            if count == 1
                deff = [deff ; Shocks_filt(i,1)];
            end
            if count ==20
                count = 0;
            end
            count = count + 1;
        end
        Shocks_filt = deff;
        clear count deff shock i
        %Rewards selection
        Rewards_filt = Restrict([leftvalve ; rightvalve],rewardTS_run ./1000);
        disp('Uploading DLC outputs')
        camara = ((camara(:,2)-camara(:,1))/2)+camara(:,1);
        % periods of movment during each condition
        if rewardTS_run(1) < aversiveTS_run(1)
            load('laps1.mat','posx','posy');
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
        %Video sampling rate 
        dt = (mean(diff(behavior.pos.aversive(:,1)))); 
        
        % Generation of no-movements periods
        % Reward
        start = behavior.speed.reward(1,1);   stop = behavior.speed.reward(end,1);
        movement.reward = InvertIntervals(behavior.quiet.reward , start , stop); %keep only those higher than criteria
        movement.reward(movement.reward(:,2) - movement.reward(:,1) <1,:)=[];
        clear tmp start stop
        % Aversive
        start = behavior.speed.aversive(1,1);   stop = behavior.speed.aversive(end,1);
        movement.aversive = InvertIntervals(behavior.quiet.aversive , start , stop);%keep only those higher than criteria
        movement.aversive(movement.aversive(:,2) - movement.aversive(:,1) <1,:)=[];
        clear tmp start stop
        
        %Restrict pos to inside maze and only movement periods
        %---Aversive----
        pos_ave = behavior.pos.aversive(:,1:2); 
        pos_ave(:,2) = pos_ave(:,2)-min(pos_ave(:,2)); pos_ave(:,2) = pos_ave(:,2)/max(pos_ave(:,2)); %normalization of position 
        in_maze = ToIntervals(pos_ave(:,1),and(pos_ave(:,2)>0.09 , pos_ave(:,2)<1-0.09));% eliminating the extrems of the maze
                    
       %Check lenght of in_maze segments to define real laps
        in_lapA = []; % keeps only full laps (maybe too restrictive)
        for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_ave(pos_ave(:,1)==ini,2); 
                        fin_pos = pos_ave(pos_ave(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>0.6
                             in_lapA = [in_lapA;in_maze(ix,:)]; 
                        end
        end
        
        pos_ave = Restrict(pos_ave,in_lapA);pos_ave = Restrict(pos_ave, movement.aversive); % Restrict pos to movement periods and laps
       
         %---Rewarded----
        pos_rew = behavior.pos.reward(:,1:2);
        pos_rew(:,2) = pos_rew(:,2)-min(pos_rew(:,2)); pos_rew(:,2) = pos_rew(:,2)/max(pos_rew(:,2)); %normalization of position 
        in_maze = ToIntervals(pos_rew(:,1),and(pos_rew(:,2)>0.09 , pos_rew(:,2)<1-0.09));% eliminating the extrems of the maze (10cm-reward box)
               
        %Check lenght of in_maze segments to define real laps
        in_lapR = []; 
        for ix=1:size(in_maze,1)
                        ini = in_maze(ix,1);
                        fin = in_maze(ix,2);
                        ini_pos = pos_rew(pos_rew(:,1)==ini,2); 
                        fin_pos = pos_rew(pos_rew(:,1)==fin,2); 
                        if abs(fin_pos-ini_pos)>0.6
                             in_lapR = [in_lapR;in_maze(ix,:)]; 
                        end
        end 
        pos_rew = Restrict(pos_rew,in_lapR); pos_rew = Restrict(pos_rew , movement.reward);
        
        
        %% Spikes
        %Load Units
        disp('Uploading Spiking activity')
        cd 'Spikesorting'
        spks = double([readNPY('spike_clusters.npy') readNPY('spike_times.npy')]);
        K = tdfread('cluster_group.tsv'); % import clusters ID and groups (MUA,Noise,Good)
        Kinfo = tdfread('cluster_info.tsv'); % import information of clusters
        K = [K.cluster_id(K.group(:,1) == 'g') , Kinfo.ch(K.group(:,1) == 'g'),]; % defining good clusters
        % Load neuronal classification
        load('Cell_type_classification')
        K = [K , Cell_type_classification(:,6:7)];
        group_dHPC = K(K(:,2) > 63,:);
        group_vHPC = K(K(:,2) <= 63,:);
        %Loop to select dorsal or ventral LFP and SU
        % z=1 --> dorsal
        % z=2 --> ventral
        for z = 1:2
            if z == 1
                n_SU_D = [n_SU_D ; length(group_dHPC)];
                spks_dHPC = spks(ismember(spks(:,1),group_dHPC(:,1)),:); %keep spks from good clusters
            else
                n_SU_V = [n_SU_V ; length(group_vHPC)];
                spks_vHPC = spks(ismember(spks(:,1),group_vHPC(:,1)),:); %keep spks from good clusters
            end
        end
        clear z
        spks_vHPC(:,2) = double(spks_vHPC(:,2))./20000;
        spks_dHPC(:,2) = double(spks_dHPC(:,2))./20000;
        % Selection of celltype to analyze
        if criteria_type == 0 %pyr
            cellulartype = [K(:,1) , K(:,4)];
        elseif criteria_type == 1 % int
            cellulartype = [K(:,1) , not(K(:,4))];
        elseif criteria_type == 2 % all
            cellulartype = [K(:,1) , ones(length(K),1)];
        end
        % Definition of the emotional transition
        if aversiveTS_run(1)<rewardTS_run(1)
            cond = 1; % 1 if the order was Aversive -> Reward
        else
            cond = 2;% 1 if the order was Reward -> Aversive
        end
        
        %% Firing maps calculation
        disp('dHPC Firing rate map calculation')
        
        %%%%%%%% DHPC %%%%%%%%%%%
        try
            load('dHPC_pc_skaggs_circular.mat');
            dhpc_pc = [];
            for d=1:size(dHPC,2)
                pc = dHPC{d}; 
                dhpc_pc = [dhpc_pc;pc.id]; 
            end
            group_dHPC = group_dHPC(ismember(group_dHPC(:,1),dhpc_pc),:); % keep only pc 
            
            if size(group_dHPC,1)>=n_session 
                %Create rate map stacks between  
                for n=1:size(group_dHPC,1)
                    pv_ave_dhpc(:,:,n) = dHPC{n}.frMap_ave; 
                    pv_rew_dhpc(:,:,n) = dHPC{n}.frMap_rew;
                end
                
                %Create rate map stacks within  
                for n=1:size(group_dHPC,1)
                    cluster = group_dHPC(n,1);
                    spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster

                    % --- Aversive ---
                    spks_ave = spks; 
                    %Restrict spk and position to laps and movement periods
                    spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps
               
                    % --- Reward ---
                    spks_rew = spks; 
                    %Restrict spk  to laps and movement periods
                    spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
               
                    [~,~,ave_1,ave_2] = Within_pc(pos_ave,spks_ave,1,sigma,Xedges);
                    %Within reward 
                    [~,~,rew_1,rew_2] = Within_pc(pos_rew,spks_rew,1,sigma,Xedges);
                    
                    pv_ave1_dhpc(:,:,n) = ave_1; 
                    pv_ave2_dhpc(:,:,n) = ave_2;
                    
                    pv_rew1_dhpc(:,:,n) =rew_1; 
                    pv_rew2_dhpc(:,:,n) =rew_2;
                    
                end
                
                %Vector bin correlation
                
                %Between
                sesion_corr = NaN(60,1);
                for c=1:size(pv_ave_dhpc,2) %Iterate through all the bins
                    %Create a bin vector (with the data of all neurons in that bin)
                    pv_ave = squeeze(pv_ave_dhpc(1,c,:));
            
                    pv_rew = squeeze(pv_rew_dhpc(1,c,:));
            
                    %Correlation between ave vector bin and rew vector bin
                    cor = corrcoef(pv_ave, pv_rew,'Rows','pairwise');
                    sesion_corr(c,1) = cor(1,2);    
                end 
                temp_dhpc_bet = [temp_dhpc_bet;sesion_corr];
                
                %Within ave
                sesion_corr = NaN(60,1);
                for c=1:size(pv_ave_dhpc,2) %Iterate through all the bins
                    %Create a bin vector (with the data of all neurons in that bin)
                    pv_ave1 = squeeze(pv_ave1_dhpc(1,c,:));
            
                    pv_ave2 = squeeze(pv_ave2_dhpc(1,c,:));
            
                    %Correlation between ave vector bin and rew vector bin
                    cor = corrcoef(pv_ave1, pv_ave2,'Rows','pairwise');
                    sesion_corr(c,1) = cor(1,2);    
                end 
                temp_dhpc_within_ave = [temp_dhpc_within_ave;sesion_corr];
                
                %Within rew
                sesion_corr = NaN(60,1);
                for c=1:size(pv_rew_dhpc,2) %Iterate through all the bins
                    %Create a bin vector (with the data of all neurons in that bin)
                    pv_rew1 = squeeze(pv_rew1_dhpc(1,c,:));
            
                    pv_rew2 = squeeze(pv_rew2_dhpc(1,c,:));
            
                    %Correlation 
                    cor = corrcoef(pv_rew1, pv_rew2,'Rows','pairwise');
                    sesion_corr(c,1) = cor(1,2);    
                end 
                temp_dhpc_within_rew = [temp_dhpc_within_rew;sesion_corr];
            end 
  
        catch 
           disp(['No dHPC_pc.mat file in ',session]);
        end 
        pv_corr.dhpc_bet = [pv_corr.dhpc_bet;temp_dhpc_bet];
        pv_corr.dhpc_within_ave = [pv_corr.dhpc_within_ave;temp_dhpc_within_ave];
        pv_corr.dhpc_within_rew = [pv_corr.dhpc_within_rew;temp_dhpc_within_rew];
       
        %%%%%%%% VHPC %%%%%%%%%%%   
        disp('vHPC Firing rate map calculation')
         
        try
            load('vHPC_pc_skaggs_circular.mat');
            vhpc_pc = [];
            for d=1:size(vHPC,2)
                pc = vHPC{d}; 
                vhpc_pc = [vhpc_pc;pc.id]; 
            end
            group_vHPC = group_vHPC(ismember(group_vHPC(:,1),vhpc_pc),:); % keep only pc 
            
            if size(group_vHPC,1)>=n_session 
                %Create rate map stacks between  
                for n=1:size(group_vHPC,1)
                    pv_ave_vhpc(:,:,n) = vHPC{n}.frMap_ave; 
                    pv_rew_vhpc(:,:,n) = vHPC{n}.frMap_rew;
                end
                
                %Create rate map stacks within  
                for n=1:size(group_vHPC,1)
                    cluster = group_vHPC(n,1);
                    spks = spks_vHPC(spks_vHPC(:,1)==cluster,2); % select tspk from cluster

                    % --- Aversive ---
                    spks_ave = spks; 
                    %Restrict spk and position to laps and movement periods
                    spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps
               
                    % --- Reward ---
                    spks_rew = spks; 
                    %Restrict spk  to laps and movement periods
                    spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
               
                    [~,~,ave_1,ave_2] = Within_pc(pos_ave,spks_ave,1,sigma,Xedges);
                    %Within reward 
                    [~,~,rew_1,rew_2] = Within_pc(pos_rew,spks_rew,1,sigma,Xedges);
                    
                    pv_ave1_vhpc(:,:,n) = ave_1; 
                    pv_ave2_vhpc(:,:,n) = ave_2;
                    
                    pv_rew1_vhpc(:,:,n) =rew_1; 
                    pv_rew2_vhpc(:,:,n) =rew_2;
                    
                end
                
                %Vector bin correlation
                
                %Between
                sesion_corr = NaN(60,1);
                for c=1:size(pv_ave_vhpc,2) %Iterate through all the bins
                    %Create a bin vector (with the data of all neurons in that bin)
                    pv_ave = squeeze(pv_ave_vhpc(1,c,:));
            
                    pv_rew = squeeze(pv_rew_vhpc(1,c,:));
            
                    %Correlation between ave vector bin and rew vector bin
                    cor = corrcoef(pv_ave, pv_rew,'Rows','pairwise');
                    sesion_corr(c,1) = cor(1,2);    
                end 
                temp_vhpc_bet = [temp_vhpc_bet;sesion_corr];
                
                %Within ave
                sesion_corr = NaN(60,1);
                for c=1:size(pv_ave_vhpc,2) %Iterate through all the bins
                    %Create a bin vector (with the data of all neurons in that bin)
                    pv_ave1 = squeeze(pv_ave1_vhpc(1,c,:));
            
                    pv_ave2 = squeeze(pv_ave2_vhpc(1,c,:));
            
                    %Correlation between ave vector bin and rew vector bin
                    cor = corrcoef(pv_ave1, pv_ave2,'Rows','pairwise');
                    sesion_corr(c,1) = cor(1,2);    
                end 
                temp_vhpc_within_ave = [temp_vhpc_within_ave;sesion_corr];
                
                %Within rew
                sesion_corr = NaN(60,1);
                for c=1:size(pv_rew_vhpc,2) %Iterate through all the bins
                    %Create a bin vector (with the data of all neurons in that bin)
                    pv_rew1 = squeeze(pv_rew1_vhpc(1,c,:));
            
                    pv_rew2 = squeeze(pv_rew2_vhpc(1,c,:));
            
                    %Correlation 
                    cor = corrcoef(pv_rew1, pv_rew2,'Rows','pairwise');
                    sesion_corr(c,1) = cor(1,2);    
                end 
                temp_vhpc_within_rew = [temp_vhpc_within_rew;sesion_corr];
            end 
  
        catch 
           disp(['No vHPC_pc.mat file in ',session]);
        end 
        pv_corr.vhpc_bet = [pv_corr.vhpc_bet;temp_vhpc_bet];
        pv_corr.vhpc_within_ave = [pv_corr.vhpc_within_ave;temp_vhpc_within_ave];
        pv_corr.vhpc_within_rew = [pv_corr.vhpc_within_rew;temp_vhpc_within_rew];
        

   
  

        clear specificity_dHPC_A specificity_dHPC_R specificity_vHPC_A specificity_vHPC_R
        clear field_dHPC_A field_dHPC_R field_vHPC_A field_vHPC_R
        clear peak_dHPC_A peak_dHPC_R peak_vHPC_A peak_vHPC_R
        clear maps_dHPC_A maps_dHPC_R maps_vHPC_A maps_vHPC_R
        clear aversiveTS aversiveTS_run rewardTS rewardTS_run baselineTS
        clear behavior camara Cell_type_classification cellulartype cluster
        clear group_dHPC group_vHPC K Kinfo leftvalve rightvalve movement R
        clear Rewards_filt Shocks_filt spks_dHPC spks_vHPC r rZ
        clear cluster cluster_dHPC cluster_vHPC coordinated coordinatedV coordinatedV_refined
        clear camaraA count dX dX_int dY dY_int gc grps i ii MUA segments tmp WAKE REM NREM
        disp(['-- Finished folder #' , num2str(t) , ' from rat #' , num2str(tt) , ' ---'])
        disp('  ')
    end
    disp(['-------- Finished rat#' , num2str(tt) , ' --------'])
    disp(' ')
end

% Plots 
figure(1);clf;hold on;
% sgtitle('Population vector')
h = cdfplot(pv_corr.dhpc_bet);hold on;
h.Color = 'black';
h.LineWidth = 2;

h=cdfplot(pv_corr.dhpc_within_ave);hold on;
h.Color = [1 0.2 0.2];
h.LineWidth = 2;

h=cdfplot(pv_corr.dhpc_within_rew);hold on;
h.Color = [0.5 0.3 1];
h.LineWidth = 2;

ylabel('Cumulative freq');
xlabel('Spatial bin corr. coef.(r)')
title('Mean population vector - dHPC')
legend('Between','Within ave', 'Within rew')
% vHPC
figure(2);clf;hold on;
% sgtitle('Population vector')
h = cdfplot(pv_corr.vhpc_bet);hold on;
h.Color = 'black';
h.LineWidth = 2;

h=cdfplot(pv_corr.vhpc_within_ave);hold on;
h.Color = [1 0.2 0.2];
h.LineWidth = 2;

h=cdfplot(pv_corr.vhpc_within_rew);hold on;
h.Color = [0.5 0.3 1];
h.LineWidth = 2;

ylabel('Cumulative freq');
xlabel('Spatial bin corr. coef.(r)')
title('Mean population vector - dHPC')


%vHPC
figure(2);clf;hold on;
% sgtitle('Population vector')
h = cdfplot(pv_corr.vhpc_bet);hold on;
h.Color = 'black';
h.LineWidth = 2;

h=cdfplot(pv_corr.vhpc_within_ave);hold on;
h.Color = [1 0.2 0.2];
h.LineWidth = 2;

h=cdfplot(pv_corr.vhpc_within_rew);hold on;
h.Color = [0.5 0.3 1];
h.LineWidth = 2;

ylabel('Cumulative freq');
xlabel('Spatial bin corr. coef.(r)')
title('Mean population vector - dHPC')
legend('Between','Within ave', 'Within rew')
% vHPC
figure(2);clf;hold on;
% sgtitle('Population vector')
h = cdfplot(pv_corr.vhpc_bet);hold on;
h.Color = 'black';
h.LineWidth = 2;

h=cdfplot(pv_corr.vhpc_within_ave);hold on;
h.Color = [1 0.2 0.2];
h.LineWidth = 2;

h=cdfplot(pv_corr.vhpc_within_rew);hold on;
h.Color = [0.5 0.3 1];
h.LineWidth = 2;

ylabel('Cumulative freq');
xlabel('Spatial bin corr. coef.(r)')
title('Mean population vector - vHPC')