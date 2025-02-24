clear
clc
close all
%% PARAMETERS - run this section before run any other
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'; ...
    '\\Maryjackson\e\Rat126\Ephys\in_Pyr'};%List of folders from the path
% path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat103\usable';'E:\Rat132\recordings\in_pyr'};%List of folders from the path
% for SU
criteria_fr = 0.01; %criteria to include or not a SU into the analysis
criteria_n = 6; % minimal number of neurons from each structure
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
n_SU_V = [];
n_SU_D = [];

binSize = 0.001; % bin size for replay events detection
bin_size = 1; % to bin pos ans spks in between/within  
% Behavior
minimal_speed = 2.5;% minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods
%Place cells map
%Parameters for FiringCurve.m
min_time = 0.15; %min time in a given bin to be consider
min_peak = 0.2; %
min_size = 4;% min number of bins for the place fields
Xedges = 60; %number of bins for RateMap construction - 2.5cm bin (because of the extremes cutting)
sigma = 2;%round(15/(180/Xedges)); %defined for gauss kernel of 15cm

%% VELOCITY / FR CORRELATION 
%Spesific parameters
 bin_size = 1; %sec
%Output
dhpc_r_all= []; 
vhpc_r_all= []; 

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
        %Rewards selection -c1:open valve c2: close valve
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
        
        %%%%%%%%% Binned velocity 
        
        %Aversive
        start = behavior.speed.aversive(1,1);   stop = behavior.speed.aversive(end,1);
        bin_limits = [start:bin_size:stop]; 
        [mov_ave] = InIntervals(bin_limits,behavior.quiet.aversive); 
        mov_ave=~mov_ave;
        bin_velocity_ave= nan(size(bin_limits,2),1); 
        for v = 1:size(bin_velocity_ave,1)
            
          if v == size(bin_velocity_ave,1) % for the last bin
              ini = bin_limits(v); fin = stop; 
              [status,~,~] = InIntervals(behavior.speed.aversive(:,1),[ini,fin]);
              bin_velocity_ave(v)=nanmean(behavior.speed.aversive(status,2));
          else 
              
          ini = bin_limits(v); fin = bin_limits(v+1); 
          [status,~,~] = InIntervals(behavior.speed.aversive(:,1),[ini,fin]);
          bin_velocity_ave(v)=nanmean(behavior.speed.aversive(status,2));
          end
         end
        
        bin_velocity_ave=[bin_limits', bin_velocity_ave];
        mean_velocity_ave = nanmean(behavior.speed.aversive(:,2));
        
        %Reward
        start = behavior.speed.reward(1,1);   stop = behavior.speed.reward(end,1);
        bin_limits = [start:bin_size:stop]; 
        [mov_rew] = InIntervals(bin_limits,behavior.quiet.reward); 
        mov_rew=~mov_rew;
        bin_velocity_rew= nan(size(bin_limits,2),1); 
        for v = 1:size(bin_velocity_rew,1)
            
          if v == size(bin_velocity_rew,1)
              ini = bin_limits(v); fin = stop; 
              [status,~,~] = InIntervals(behavior.speed.reward(:,1),[ini,fin]);
              bin_velocity_rew(v)=nanmean(behavior.speed.reward(status,2));
          else 
          ini = bin_limits(v); fin = bin_limits(v+1); 
          [status,~,~] = InIntervals(behavior.speed.reward(:,1),[ini,fin]);
          bin_velocity_rew(v)=nanmean(behavior.speed.reward(status,2));
          end
        end
        
        bin_velocity_rew=[bin_limits', bin_velocity_rew];
        mean_velocity_rew = nanmean(behavior.speed.reward(:,2));
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

        %% Fr/velocity correlation 
        
        %%%%%%% PC
        try
            load('dHPC_pc_skaggs_circular.mat');
            dhpc_pc = [];
            for d=1:size(dHPC,2)
                pc = dHPC{d}; 
                dhpc_pc = [dhpc_pc;pc.id]; 
            end
            
            group_dHPC = group_dHPC(ismember(group_dHPC(:,1),dhpc_pc),:); % keep only pc 
             
            for d=1:size(group_dHPC,1)
                
                cluster = group_dHPC(d,1);
                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster
                
                %Binned spikes
                [dN,~]=binspikes(spks,1/bin_size,[behavior.speed.aversive(1,1) behavior.speed.aversive(end,1)]); 
                binned_fr = dN./bin_size; binned_fr =binned_fr(1:end-1);
                
%                 [R_ave,~] = corrcoef(bin_velocity_ave(mov_ave,2),binned_fr(mov_ave)); R_ave=R_ave(1,2);
%                 
%                 figure(2);clf; hold on; scatter(bin_velocity_ave(mov_ave,2),binned_fr(mov_ave), 'filled'); xlabel('Velocity'); ylabel('Fr')
                dhpc_table = array2table([bin_velocity_ave(mov_ave,2),binned_fr(mov_ave)],'VariableNames',{'bin velocity','bin fr'}); 
                mdl = fitlm(dhpc_table);R_ave= mdl.Rsquared.Ordinary; 
%                 figure(1);plot(mdl)

                %Binned spikes
                [dN,~]=binspikes(spks,1/bin_size,[behavior.speed.reward(1,1) behavior.speed.reward(end,1)]); 
                binned_fr = dN./bin_size; binned_fr =binned_fr(1:end-1);
                
%                 [R_rew,~] = corrcoef(bin_velocity_rew(mov_rew,2),binned_fr(mov_rew));R_rew= R_rew(1,2);
%               figure(); hold on; scatter(bin_velocity_rew(mov_rew,2),binned_fr(mov_rew), 'filled'); xlabel('Velocity'); ylabel('Fr');
                

                dhpc_table = array2table([bin_velocity_rew(mov_rew,2),binned_fr(mov_rew)],'VariableNames',{'bin velocity','bin fr'}); 
                mdl = fitlm(dhpc_table);R_rew= mdl.Rsquared.Ordinary;
%                 figure; plot(mdl)
                
                dhpc_r_all= [dhpc_r_all;R_ave, R_rew]; 
                
            end

        catch 
            dHPC = []; disp(['No dHPC_pc.mat file in ',session]);
        end 
        
        try
            load('vHPC_pc_skaggs_circular.mat');
            vhpc_pc = [];
            for d=1:size(vHPC,2)
               pc = vHPC{d}; 
               vhpc_pc = [vhpc_pc;pc.id];
            end
            group_vHPC = group_vHPC(ismember(group_vHPC(:,1),vhpc_pc),:); % keep only pc 
            for d=1:size(group_vHPC,1)
                
                cluster = group_vHPC(d,1);
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2); % select tspk from cluster
                
                %Binned spikes
                [dN,~]=binspikes(spks,1/bin_size,[behavior.speed.aversive(1,1) behavior.speed.aversive(end,1)]); 
                binned_fr = dN./bin_size; binned_fr =binned_fr(1:end-1);
                
%                 [R_ave,~] = corrcoef(bin_velocity_ave(mov_ave,2),binned_fr(mov_ave)); R_ave=R_ave(1,2);
      
%                 figure; hold on; scatter(bin_velocity',binned_fr(1:end-1) , 'filled'); xlabel('Velocity'); ylabel('Fr')
               
                table = array2table([bin_velocity_ave(mov_ave,2),binned_fr(mov_ave)],'VariableNames',{'bin velocity','bin fr'}); 
                mdl = fitlm(table);R_ave= mdl.Rsquared.Ordinary; 
                
                %Binned spikes
                [dN,~]=binspikes(spks,1/bin_size,[behavior.speed.reward(1,1) behavior.speed.reward(end,1)]); 
                binned_fr = dN./bin_size;binned_fr =binned_fr(1:end-1);
                
%                 [R_rew,~] = corrcoef(bin_velocity_rew(mov_rew,2),binned_fr(mov_rew)); R_rew=R_rew(1,2);
%                 figure; hold on; scatter(bin_velocity_rew',binned_fr(1:end-1) , 'filled'); xlabel('Velocity'); ylabel('Fr')

                table = array2table([bin_velocity_rew(mov_rew,2),binned_fr(mov_rew)],'VariableNames',{'bin velocity','bin fr'}); 
                mdl = fitlm(table);R_rew= mdl.Rsquared.Ordinary;
                
                vhpc_r_all= [vhpc_r_all;R_ave,R_rew]; 
                
            end
            
        catch
            vHPC= []; disp(['No vHPC_pc.mat file in ',session]);
        end


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
        clear ripple_bursts ripplesD ripplesV x w position_shocks posx posy ejeX ejeY PC replay Replay limits
        disp(['-- Finished folder #' , num2str(t) , ' from rat #' , num2str(tt) , ' ---'])
        disp('  ')
    end
    disp(['-------- Finished rat#' , num2str(tt) , ' --------'])
    disp(' ')
end


%%%%%Plot

%Cumilative 
figure(6);clf;hold on; cdfplot([dhpc_r_all(:,1);dhpc_r_all(:,2)]); 
cdfplot([vhpc_r_all(:,1);vhpc_r_all(:,2)]); 

%Correlational plot 
figure(1);clf; hold on;title('dhpc-green vhpc-blue');
scatter(dhpc_r_all(:,1),dhpc_r_all(:,2),[],[0 1 0],'filled','g'); xlabel('Rave'); ylabel('Rrew'); ylim([0 0.8]);xlim([0 0.6]);
scatter(vhpc_r_all(:,1),vhpc_r_all(:,2),[],[0 0 1],'filled','b'); xlabel('Rave'); ylabel('Rrew'); ylim([0 0.8]);xlim([0 0.6]);

% plot([-0.2 0.6],[-0.2 0.8]);

figure(2);clf; hold on; 
subplot(1,2,1);mdl_dhpc = fitlm(array2table([dhpc_r_all(:,1),dhpc_r_all(:,2)])); plot(mdl_dhpc); title('dHPC');xlabel('Rave'); ylabel('Rrew');
subplot(1,2,2);mdl_vhpc = fitlm(array2table([vhpc_r_all(:,1),vhpc_r_all(:,2)]));plot(mdl_vhpc); title('vHPC');xlabel('Rave'); ylabel('Rrew');

%%%%% Stats
[h,p] = kstest([dhpc_r_all(:,1);dhpc_r_all(:,2)]) %  no normal data 
[h,p] = kstest([vhpc_r_all(:,1);vhpc_r_all(:,2)]) %  no normal data 

[p,h] = signrank(dhpc_r_all(:,1),dhpc_r_all(:,2))
[p, h] = signrank(dhpc_r_all(:,1)); % h=1, p sig, reject null hypothesis -> data differ from 0
[p, h] = signrank(dhpc_r_all(:,2)); % h=1, p sig, reject null hypothesis -> data differ from 0

[p,h] = signrank(vhpc_r_all(:,1),vhpc_r_all(:,2))
[p, h] = signrank(vhpc_r_all(:,1)); % h=0,p no sig, data do not differ from 0
[p, h] = signrank(vhpc_r_all(:,2)); % h=1, p  sig, data  differ from 0

% Scatter plots 
figure(2);clf;hold on
sgtitle('Firing rate - velocity')
subplot(1,2,1); title('dHPC');
x = [ones(size(dhpc_r_all,1),1);ones(size(dhpc_r_all,1),1)*2];
y= [dhpc_r_all(:,1);dhpc_r_all(:,2)];

scatter(x,y,30,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0.5 2.5]); ylim([-0.2 0.8]);
ylabel('r2');
xticks([1 2]); xticklabels({'Ave', 'Rew'});
hold on
s1=scatter(1,nanmedian(dhpc_r_all(:,1)), "filled");s1.MarkerFaceColor = [1 0 0];
s2=scatter(2,nanmedian(dhpc_r_all(:,2)), "filled");s2.MarkerFaceColor = [0 0 1];
 
subplot(1,2,2); title('vHPC');
x = [ones(size(vhpc_r_all,1),1);ones(size(vhpc_r_all,1),1)*2];
y= [vhpc_r_all(:,1);vhpc_r_all(:,2)];

scatter(x,y,30,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0.5 2.5]); ylim([-0.2 0.8]);
ylabel('r2');
xticks([1 2]); xticklabels({'Ave', 'Rew'});
hold on
s1=scatter(1,nanmedian(vhpc_r_all(:,1)), "filled");s1.MarkerFaceColor = [1 0 0];
s2=scatter(2,nanmedian(vhpc_r_all(:,2)), "filled");s2.MarkerFaceColor = [0 0 1];

%%%%%Ratio r2 ave/ r2 rew
figure(3);clf;hold on
sgtitle('R2 ratio')
x = [ones(size(dhpc_r_all,1),1);ones(size(vhpc_r_all,1),1)*2];
y= [dhpc_r_all(:,1)./dhpc_r_all(:,2);vhpc_r_all(:,1)./vhpc_r_all(:,2)];

scatter(x,y,30,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0.5 2.5]); ylim([-100 100]);
ylabel('r2 ratio');
xticks([1 2]); xticklabels({'dHPC', 'vHPC'});
hold on
s1=scatter(1,nanmedian(dhpc_r_all(:,1)./dhpc_r_all(:,2)), "filled");s1.MarkerFaceColor = [1 0 0];
s2=scatter(2,nanmedian(vhpc_r_all(:,1)./vhpc_r_all(:,2)), "filled");s2.MarkerFaceColor = [0 0 1];
yline(0); 

%stats 
[p,h] = ranksum(dhpc_r_all(:,1)./dhpc_r_all(:,2),vhpc_r_all(:,1)./vhpc_r_all(:,2))
[p, h] = signrank(dhpc_r_all(:,1)./dhpc_r_all(:,2)); % h=1, p sig, reject null hypothesis -> data differ from 0
[p, h] = signrank(vhpc_r_all(:,1)./vhpc_r_all(:,2)); % h=1, p sig, reject null hypothesis -> data differ from 0

%% SUBSTRACT NEURONS WITH HIGH VELOCITY/FR CORRELATION 
% using output of previous section
 figure(1);clf; 
 sgtitle('dHPC')
 Q = quantile(dhpc_r_all(:,1),[0.025 0.25 0.5 0.75 0.975]); 
 subplot(1,2,1);h=histogram(dhpc_r_all(:,1));h.EdgeColor = 'none';h.FaceColor=[.3,.3,.3];xlim([0 0.4]);ylim([0 450]);
 xl = xline(Q(4),'r');
 xlabel('r2 aversive'); ylabel('Counts');legend(xl,'Q 0.75')
 subplot(1,2,2);h=histogram(dhpc_r_all(:,2));xlim([0 0.4]);ylim([0 450]);h.EdgeColor = 'none';h.FaceColor=[.3,.3,.3]; 
 xline(Q(4),'r'); xlabel('r2 reward'); ylabel('Counts');
 tresh_dhpc = Q(4); 
 
 figure(1);clf; 
 sgtitle('vHPC')
 subplot(1,2,1);h= histogram(vhpc_r_all(:,1));h.EdgeColor = 'none';h.FaceColor=[.3,.3,.3];xlim([0 0.2]);ylim([0 130])
 Q = quantile(vhpc_r_all(:,1),[0.025 0.25 0.5 0.75 0.975]); 
 xl = xline(Q(4),'r'); xlabel('r2 aversive'); ylabel('Counts');legend(xl,'Q 0.75')
 subplot(1,2,2);h= histogram(vhpc_r_all(:,2));h.EdgeColor = 'none';h.FaceColor=[.3,.3,.3];xlim([0 0.2]);ylim([0 130]);
 xline(Q(4),'r'); xlabel('r2 reward'); ylabel('Counts');
 tresh_vhpc = Q(4); 
 

%Spesific parameters
 bin_size = 1; %sec
%Output
dhpc_params= []; 
vhpc_params= []; 

% not the fastest way of doing this 
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
        %Rewards selection -c1:open valve c2: close valve
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
        
        %%%%%%%%% Binned velocity 
        
        %Aversive
        start = behavior.speed.aversive(1,1);   stop = behavior.speed.aversive(end,1);
        bin_limits = [start:bin_size:stop]; 
        [mov_ave] = InIntervals(bin_limits,behavior.quiet.aversive); 
        
        bin_velocity_ave= nan(size(bin_limits,2),1); 
        for v = 1:size(bin_velocity_ave,1)
            
          if v == size(bin_velocity_ave,1)
              ini = bin_limits(v); fin = stop; 
              [status,~,~] = InIntervals(behavior.speed.aversive(:,1),[ini,fin]);
              bin_velocity_ave(v)=nanmean(behavior.speed.aversive(status,2));
          else 
          ini = bin_limits(v); fin = bin_limits(v+1); 
          [status,~,~] = InIntervals(behavior.speed.aversive(:,1),[ini,fin]);
          bin_velocity_ave(v)=nanmean(behavior.speed.aversive(status,2));
          end
         end
        
        bin_velocity_ave=[bin_limits', bin_velocity_ave];
        mean_velocity_ave = nanmean(behavior.speed.aversive(:,2));
        
        %Reward
        start = behavior.speed.reward(1,1);   stop = behavior.speed.reward(end,1);
        bin_limits = [start:bin_size:stop]; 
        [mov_rew] = InIntervals(bin_limits,behavior.quiet.reward); 
        bin_velocity_rew= nan(size(bin_limits,2),1); 
        for v = 1:size(bin_velocity_rew,1)
            
          if v == size(bin_velocity_rew,1)
              ini = bin_limits(v); fin = stop; 
              [status,~,~] = InIntervals(behavior.speed.reward(:,1),[ini,fin]);
              bin_velocity_rew(v)=nanmean(behavior.speed.reward(status,2));
          else 
          ini = bin_limits(v); fin = bin_limits(v+1); 
          [status,~,~] = InIntervals(behavior.speed.reward(:,1),[ini,fin]);
          bin_velocity_rew(v)=nanmean(behavior.speed.reward(status,2));
          end
        end
        
        bin_velocity_rew=[bin_limits', bin_velocity_rew];
        mean_velocity_rew = nanmean(behavior.speed.reward(:,2));
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

        %% Fr/velocity correlation 
        
        %%%%%%% PC
        try
            load('dHPC_pc_skaggs_circular.mat');
            dhpc_pc = [];
            for d=1:size(dHPC,2)
                pc = dHPC{d}; 
                dhpc_pc = [dhpc_pc;pc.id]; 
            end
            
            group_dHPC = group_dHPC(ismember(group_dHPC(:,1),dhpc_pc),:); % keep only pc 
             
            for d=1:size(group_dHPC,1)
                
                cluster = group_dHPC(d,1);
                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster
                
                %Binned spikes
                [dN,~]=binspikes(spks,1/bin_size,[behavior.speed.aversive(1,1) behavior.speed.aversive(end,1)]); 
                binned_fr = dN./bin_size; binned_fr =binned_fr(1:end-1);
                
                [R_ave,~] = corrcoef(bin_velocity_ave(mov_ave,2),binned_fr(mov_ave)); 
      
%                 figure(4); hold on; scatter(bin_velocity_ave(mov_ave,2),binned_fr(mov_ave), 'filled'); xlabel('Velocity'); ylabel('Fr')
               
                %Binned spikes
                [dN,~]=binspikes(spks,1/bin_size,[behavior.speed.reward(1,1) behavior.speed.reward(end,1)]); 
                binned_fr = dN./bin_size; binned_fr =binned_fr(1:end-1);
                
                [R_rew,~] = corrcoef(bin_velocity_rew(mov_rew,2),binned_fr(mov_rew)); 
%                 figure; hold on; scatter(bin_velocity_rew',binned_fr(1:end-1) , 'filled'); xlabel('Velocity'); ylabel('Fr')
                
                if R_ave(1,2) && R_rew(1,2) < tresh_dhpc
                    dhpc_params= [dhpc_params;[dHPC{d}.between,1;dHPC{d}.within_ave,2;dHPC{d}.within_rew,3]]; 
                end 
            end

        catch 
            dHPC = []; disp(['No dHPC_pc.mat file in ',session]);
        end 
        
        try
            load('vHPC_pc_skaggs_circular.mat');
            vhpc_pc = [];
            for d=1:size(vHPC,2)
               pc = vHPC{d}; 
               vhpc_pc = [vhpc_pc;pc.id];
            end
            group_vHPC = group_vHPC(ismember(group_vHPC(:,1),vhpc_pc),:); % keep only pc 
            for d=1:size(group_vHPC,1)
                
                cluster = group_vHPC(d,1);
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2); % select tspk from cluster
                
                %Binned spikes
                [dN,~]=binspikes(spks,1/bin_size,[behavior.speed.aversive(1,1) behavior.speed.aversive(end,1)]); 
                binned_fr = dN./bin_size; binned_fr =binned_fr(1:end-1);
                
                [R_ave,~] = corrcoef(bin_velocity_ave(mov_ave,2),binned_fr(mov_ave));  
      
%                 figure; hold on; scatter(bin_velocity',binned_fr(1:end-1) , 'filled'); xlabel('Velocity'); ylabel('Fr')
               
                %Binned spikes
                [dN,~]=binspikes(spks,1/bin_size,[behavior.speed.reward(1,1) behavior.speed.reward(end,1)]); 
                binned_fr = dN./bin_size;binned_fr =binned_fr(1:end-1);
                
                [R_rew,~] = corrcoef(bin_velocity_rew(mov_rew,2),binned_fr(mov_rew)); 
%                 figure; hold on; scatter(bin_velocity_rew',binned_fr(1:end-1) , 'filled'); xlabel('Velocity'); ylabel('Fr')
                
               if R_ave(1,2) && R_rew(1,2) < tresh_vhpc
                    vhpc_params= [vhpc_params;[vHPC{d}.between,1;vHPC{d}.within_ave,2;vHPC{d}.within_rew,3]]; 
                end  
                
            end
            
        catch
            vHPC= []; disp(['No vHPC_pc.mat file in ',session]);
        end


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
        clear ripple_bursts ripplesD ripplesV x w position_shocks posx posy ejeX ejeY PC replay Replay limits
        disp(['-- Finished folder #' , num2str(t) , ' from rat #' , num2str(tt) , ' ---'])
        disp('  ')
    end
    disp(['-------- Finished rat#' , num2str(tt) , ' --------'])
    disp(' ')
end
 

%Only 3 plots
ylabels ={'Spatial correlation', 'Overlap', 'Pf shift'};
idx = [1,3,4]; 
xlabels = {'Between', 'Within Ave', 'Within Rew'}; 
figure(1);clf;hold on, 
sgtitle('dHPC')
for a=1:size(ylabels,2)
    p = idx(a);
    subplot(1,3,a); hold on; 
    x = dhpc_params(:,6);
    y= dhpc_params(:,p);%Change the column number (1-4) to choose which variable to plot 
    c = [.3,.3,.3];
    scatter(x,y,30,c,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0 4]);
    ylabel(ylabels{a});
    xticks([1 2 3])
    xticklabels({'Bet', 'WA', 'WR'});
    hold on
    s1=scatter(1,nanmedian(dhpc_params(dhpc_params(:,6)==1,p)), "filled");s1.MarkerFaceColor = [0 0 0];
    s2=scatter(2,nanmedian(dhpc_params(dhpc_params(:,6)==2,p)), "filled");s2.MarkerFaceColor = [1 0.1 0.2];
    s3=scatter(3,nanmedian(dhpc_params(dhpc_params(:,6)==3,p)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];
end

%Stats
[P,ANOVATAB,STATS] = kruskalwallis(dhpc_sub(:,5),dhpc_sub(:,6));
c = multcompare(STATS)

tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

% Paired non-parametric anova - Friedman test
%Prepare data to test
p = 1; % 1=sp corr 2=fr change 3=overlap 4=pf shift
data=[dhpc_params(dhpc_params(:,6)==1,p), dhpc_params(dhpc_params(:,6)==2,p), ...
    dhpc_params(dhpc_params(:,6)==3,p), dhpc_params(dhpc_params(:,6)==4,p),dhpc_params(dhpc_params(:,6)==5,p)];
data(any(isnan(data), 2), :) = [];

[p,~,stats] =  friedman(data,1); 
c = multcompare(stats);
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

%Only 3 plots
ylabels ={'Spatial correlation', 'Overlap', 'Pf shift'};
idx = [1,3,4]; 
xlabels = {'Between', 'Within Ave', 'Within Rew'}; 
figure(1);clf;hold on, 
sgtitle('vHPC')
for a=1:size(ylabels,2)
    p = idx(a);
    subplot(1,3,a); hold on; 
    x = vhpc_params(:,6);
    y= vhpc_params(:,p);%Change the column number (1-4) to choose which variable to plot 
    c = [.3,.3,.3];
    scatter(x,y,30,c,"filled",'jitter','on', 'jitterAmount',0.3); xlim([0 4]);
    ylabel(ylabels{a});
    xticks([1 2 3])
    xticklabels({'Bet', 'WA', 'WR'});
    hold on
    s1=scatter(1,nanmedian(vhpc_params(vhpc_params(:,6)==1,p)), "filled");s1.MarkerFaceColor = [0 0 0];
    s2=scatter(2,nanmedian(vhpc_params(vhpc_params(:,6)==2,p)), "filled");s2.MarkerFaceColor = [1 0.1 0.2];
    s3=scatter(3,nanmedian(vhpc_params(vhpc_params(:,6)==3,p)), "filled");s3.MarkerFaceColor = [0.1 0.3 1];
end

%Stats
[P,ANOVATAB,STATS] = kruskalwallis(dhpc_sub(:,5),dhpc_sub(:,6));
c = multcompare(STATS)

tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

% Paired non-parametric anova - Friedman test
%Prepare data to test
p = 1; % 1=sp corr 2=fr change 3=overlap 4=pf shift
data=[vhpc_params(vhpc_params(:,6)==1,p), vhpc_params(vhpc_params(:,6)==2,p), ...
    vhpc_params(vhpc_params(:,6)==3,p), vhpc_params(vhpc_params(:,6)==4,p),vhpc_params(vhpc_params(:,6)==5,p)];
data(any(isnan(data), 2), :) = [];

[p,~,stats] =  friedman(data,1); 
c = multcompare(stats);
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

%% VELOCITY - REMAPPING CORRELATION
%Spesific parameters
 bin_size = 1; %sec

%Output

dhpc_params= []; 
vhpc_params= []; 


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
        %Rewards selection -c1:open valve c2: close valve
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
        
        %%%%%%%%% Binned velocity 
        
        %Aversive
        start = behavior.speed.aversive(1,1);   stop = behavior.speed.aversive(end,1);
        bin_limits = [start:bin_size:stop]; 
        [mov_ave] = InIntervals(bin_limits,behavior.quiet.aversive); 
        
        bin_velocity_ave= nan(size(bin_limits,2),1); 
        for v = 1:size(bin_velocity_ave,1)
            
          if v == size(bin_velocity_ave,1)
              ini = bin_limits(v); fin = stop; 
              [status,~,~] = InIntervals(behavior.speed.aversive(:,1),[ini,fin]);
              bin_velocity_ave(v)=nanmean(behavior.speed.aversive(status,2));
          else 
          ini = bin_limits(v); fin = bin_limits(v+1); 
          [status,~,~] = InIntervals(behavior.speed.aversive(:,1),[ini,fin]);
          bin_velocity_ave(v)=nanmean(behavior.speed.aversive(status,2));
          end
         end
        
        bin_velocity_ave=[bin_limits', bin_velocity_ave];
        mean_velocity_ave = nanmean(behavior.speed.aversive(:,2));
        
        %Reward
        start = behavior.speed.reward(1,1);   stop = behavior.speed.reward(end,1);
        bin_limits = [start:bin_size:stop]; 
        [mov_rew] = InIntervals(bin_limits,behavior.quiet.reward); 
        bin_velocity_rew= nan(size(bin_limits,2),1); 
        for v = 1:size(bin_velocity_rew,1)
            
          if v == size(bin_velocity_rew,1)
              ini = bin_limits(v); fin = stop; 
              [status,~,~] = InIntervals(behavior.speed.reward(:,1),[ini,fin]);
              bin_velocity_rew(v)=nanmean(behavior.speed.reward(status,2));
          else 
          ini = bin_limits(v); fin = bin_limits(v+1); 
          [status,~,~] = InIntervals(behavior.speed.reward(:,1),[ini,fin]);
          bin_velocity_rew(v)=nanmean(behavior.speed.reward(status,2));
          end
        end
        
        bin_velocity_rew=[bin_limits', bin_velocity_rew];
        mean_velocity_rew = nanmean(behavior.speed.reward(:,2));
        %% Saving remapping parameters 
        cd 'Spikesorting'
      
        try
            load('dHPC_pc_skaggs_circular.mat');

             
            for d=1:size(dHPC,2)
                dhpc_params= [dhpc_params;[dHPC{d}.between,mean_velocity_ave, mean_velocity_rew]]; 
            end

        catch 
            dHPC = []; disp(['No dHPC_pc.mat file in ',session]);
        end 
        
        try
            load('vHPC_pc_skaggs_circular.mat');
            
            for d=1:size(vHPC,2)
                vhpc_params= [vhpc_params;[vHPC{d}.between,mean_velocity_ave, mean_velocity_rew]];
  
            end
            
        catch
            vHPC= []; disp(['No vHPC_pc.mat file in ',session]);
        end


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
        clear ripple_bursts ripplesD ripplesV x w position_shocks posx posy ejeX ejeY PC replay Replay limits
        disp(['-- Finished folder #' , num2str(t) , ' from rat #' , num2str(tt) , ' ---'])
        disp('  ')
    end
    disp(['-------- Finished rat#' , num2str(tt) , ' --------'])
    disp(' ')
end

dhpc_table = array2table([(dhpc_params(:,7)-dhpc_params(:,6))./(dhpc_params(:,7)+dhpc_params(:,6)),dhpc_params(:,3)],'VariableNames',{'Velocity change','overlap'}); 
mdl = fitlm(dhpc_table)

vhpc_table = array2table([(vhpc_params(:,7)-vhpc_params(:,6))./(vhpc_params(:,7)+vhpc_params(:,6)),vhpc_params(:,1)],'VariableNames',{'Velocity change','overlap'}); 
mdl = fitlm(vhpc_table)

plot(mdl)

%% FR MAPS NORMALIZED BY VELOCITY 
%Spesific parameters
 bin_size = 1; %sec
 
dHPC_all = [];
vHPC_all = [];

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
        %Rewards selection -c1:open valve c2: close valve
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
        
        %%%%%%%%% Time binned velocity 
        
        %Aversive
        start = behavior.speed.aversive(1,1);   stop = behavior.speed.aversive(end,1);
        bin_limits = [start:bin_size:stop]; 
        [mov_ave] = InIntervals(bin_limits,behavior.quiet.aversive); 
        
        bin_velocity_ave= nan(size(bin_limits,2),1); 
        for v = 1:size(bin_velocity_ave,1)
            
          if v == size(bin_velocity_ave,1)
              ini = bin_limits(v); fin = stop; 
              [status,~,~] = InIntervals(behavior.speed.aversive(:,1),[ini,fin]);
              bin_velocity_ave(v)=nanmean(behavior.speed.aversive(status,2));
          else 
          ini = bin_limits(v); fin = bin_limits(v+1); 
          [status,~,~] = InIntervals(behavior.speed.aversive(:,1),[ini,fin]);
          bin_velocity_ave(v)=nanmean(behavior.speed.aversive(status,2));
          end
         end
        
        bin_velocity_ave=[bin_limits', bin_velocity_ave];
       
        
        %Reward
        start = behavior.speed.reward(1,1);   stop = behavior.speed.reward(end,1);
        bin_limits = [start:bin_size:stop]; 
        [mov_rew] = InIntervals(bin_limits,behavior.quiet.reward); 
        bin_velocity_rew= nan(size(bin_limits,2),1); 
        for v = 1:size(bin_velocity_rew,1)
            
          if v == size(bin_velocity_rew,1)
              ini = bin_limits(v); fin = stop; 
              [status,~,~] = InIntervals(behavior.speed.reward(:,1),[ini,fin]);
              bin_velocity_rew(v)=nanmean(behavior.speed.reward(status,2));
          else 
          ini = bin_limits(v); fin = bin_limits(v+1); 
          [status,~,~] = InIntervals(behavior.speed.reward(:,1),[ini,fin]);
          bin_velocity_rew(v)=nanmean(behavior.speed.reward(status,2));
          end
        end
        
        bin_velocity_rew=[bin_limits', bin_velocity_rew];
        clear strt stop bin_limits mov_rew mov_ave bin_velocity_rew bin_velocity_ave ini fin status
        %% Spatial binned velocity 
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
        
        pos_ave= [pos_ave, behavior.speed.aversive(:,2)]; 
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
        pos_rew= [pos_rew, behavior.speed.reward(:,2)]; 
        pos_rew = Restrict(pos_rew,in_lapR); pos_rew = Restrict(pos_rew , movement.reward);
        
        %%%%%%%%% Spatial binned velocity 
        spbin_velocity_ave =  spatial_bin_velocity(pos_ave, Xedges); 
        spbin_velocity_rew =  spatial_bin_velocity(pos_rew, Xedges); 
         
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

        %% Fr/velocity correlation 
        
        %%%%%%% PC
        try
            load('dHPC_pc_skaggs_circular.mat');
            dhpc_pc = [];
            for d=1:size(dHPC,2)
                pc = dHPC{d}; 
                dhpc_pc = [dhpc_pc;pc.id]; 
            end
            
            group_dHPC = group_dHPC(ismember(group_dHPC(:,1),dhpc_pc),:); % keep only pc 
             
            for d=1:size(group_dHPC,1)
                
                cluster = group_dHPC(d,1);
                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster
                
                %%%Aversive 
                %Comuted r fr vs velocity 
                [dN,~]=binspikes(spks,1/bin_size,[behavior.speed.aversive(1,1) behavior.speed.aversive(end,1)]); % binned spk
                binned_fr = dN./bin_size; binned_fr =binned_fr(1:end-1);

                dhpc_table = array2table([bin_velocity_ave(mov_ave,2),binned_fr(mov_ave)],'VariableNames',{'bin velocity','bin fr'}); 
                mdl = fitlm(dhpc_table);R_ave= mdl.Rsquared.Ordinary; 
                clear dN binned_fr dhpc_table  mdl
                
                %Normalization factor 
                spbin_velocity_ave = spbin_velocity_ave.*R_ave; 
                
                %%% Reward
                %Comuted r fr vs velocity 
                [dN,~]=binspikes(spks,1/bin_size,[behavior.speed.reward(1,1) behavior.speed.reward(end,1)]); 
                binned_fr = dN./bin_size; binned_fr =binned_fr(1:end-1);

                dhpc_table = array2table([bin_velocity_rew(mov_rew,2),binned_fr(mov_rew)],'VariableNames',{'bin velocity','bin fr'}); 
                mdl = fitlm(dhpc_table);R_rew= mdl.Rsquared.Ordinary;
                clear dN binned_fr dhpc_table  mdl
                
                 %Normalization factor 
                spbin_velocity_rew = spbin_velocity_rew.*R_rew; 
                
                %Compute the firing curve 
                spks_ave = spks; 
                %Restrict spk and position to laps and movement periods
                spks_ave = Restrict(spks_ave,in_lapA);spks_ave = Restrict(spks_ave , movement.aversive); % Restrict spikes to movement periods and laps

                % --- Reward ---
                spks_rew = spks; 
                %Restrict spk and position to laps and movement periods
                spks_rew = Restrict(spks_rew,in_lapR); spks_rew = Restrict(spks_rew, movement.reward); % Restrict to movement periods
            
                %Within aversive
                [withinA] = Within_pc(pos_ave,spks_ave,1,sigma,Xedges,'velocity_norm',spbin_velocity_ave);
                %Within reward 
                [withinR] = Within_pc(pos_rew,spks_rew,1,sigma,Xedges,'velocity_norm',spbin_velocity_rew);
                %Between
                [between] = Between_pc(pos_ave,spks_ave,pos_rew,spks_rew,bin_size,sigma,Xedges,'velocity_norm_ave',spbin_velocity_ave,'velocity_norm_rew',spbin_velocity_rew);
                
                dHPC_all = [dHPC_all; between,1;withinA,2;withinR,3];
                
            end

        catch 
            dHPC = []; disp(['No dHPC_pc.mat file in ',session]);
        end 
        
        try
            load('vHPC_pc_skaggs_circular.mat');
            vhpc_pc = [];
            for d=1:size(vHPC,2)
               pc = vHPC{d}; 
               vhpc_pc = [vhpc_pc;pc.id];
            end
            group_vHPC = group_vHPC(ismember(group_vHPC(:,1),vhpc_pc),:); % keep only pc 
            for d=1:size(group_vHPC,1)
                
                cluster = group_vHPC(d,1);
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2); % select tspk from cluster
                
                %Binned spikes
                [dN,~]=binspikes(spks,1/bin_size,[behavior.speed.aversive(1,1) behavior.speed.aversive(end,1)]); 
                binned_fr = dN./bin_size; binned_fr =binned_fr(1:end-1);
                
%                 [R_ave,~] = corrcoef(bin_velocity_ave(mov_ave,2),binned_fr(mov_ave)); R_ave=R_ave(1,2);
      
%                 figure; hold on; scatter(bin_velocity',binned_fr(1:end-1) , 'filled'); xlabel('Velocity'); ylabel('Fr')
               
                table = array2table([bin_velocity_ave(mov_ave,2),binned_fr(mov_ave)],'VariableNames',{'bin velocity','bin fr'}); 
                mdl = fitlm(table);R_ave= mdl.Rsquared.Ordinary; 
                
                %Binned spikes
                [dN,~]=binspikes(spks,1/bin_size,[behavior.speed.reward(1,1) behavior.speed.reward(end,1)]); 
                binned_fr = dN./bin_size;binned_fr =binned_fr(1:end-1);
                
%                 [R_rew,~] = corrcoef(bin_velocity_rew(mov_rew,2),binned_fr(mov_rew)); R_rew=R_rew(1,2);
%                 figure; hold on; scatter(bin_velocity_rew',binned_fr(1:end-1) , 'filled'); xlabel('Velocity'); ylabel('Fr')

                table = array2table([bin_velocity_rew(mov_rew,2),binned_fr(mov_rew)],'VariableNames',{'bin velocity','bin fr'}); 
                mdl = fitlm(table);R_rew= mdl.Rsquared.Ordinary;
                
                vhpc_r_all= [vhpc_r_all;R_ave,R_rew]; 
                
            end
            
        catch
            vHPC= []; disp(['No vHPC_pc.mat file in ',session]);
        end


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
        clear ripple_bursts ripplesD ripplesV x w position_shocks posx posy ejeX ejeY PC replay Replay limits
        disp(['-- Finished folder #' , num2str(t) , ' from rat #' , num2str(tt) , ' ---'])
        disp('  ')
    end
    disp(['-------- Finished rat#' , num2str(tt) , ' --------'])
    disp(' ')
end
