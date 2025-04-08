function [Pred_shock , Pred_member , Pred_nomember] = GLM_ShockSU_one_out2(path)
% This function calculates the Pre and Post sleep assemblies rate.
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% OUTPUT
% Pre, Post: Structure, it store the SU response sourrounding the ripple.
%
% T: vector, it contains time vector for visualization
%
% I: Structure, it contains the tag to determine if the SU is or not
%    shock responssive
%
% curve: Structure, it contains the response of the neuron sourrounding the
%        shock. It is zscored.
%
% Morci Juan Facundo 08/2024

% variables to use in the script
criteria_fr = 0;
criteria_type = 0; % criteria for celltype (0:pyr, 1:int, 2:all)
sm = 0;
bin = 0.025;

Pred_shock.dHPC = [];      Pred_shock.vHPC = [];
Pred_member.dHPC = [];     Pred_member.vHPC = [];
Pred_nomember.dHPC = [];   Pred_nomember.vHPC = [];

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
        load('behavioral_data.mat', 'movement')
        
        %% load sleep states
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);
        clear x states
        
        baselineTS = baselineTS./1000;                  aversiveTS = aversiveTS./1000;                    rewardTS = rewardTS./1000;
        aversiveTS_run = aversiveTS_run./1000;          rewardTS_run = rewardTS_run./1000;
        
        NREM.baseline = Restrict(NREM.all,baselineTS);   NREM.aversive = Restrict(NREM.all,aversiveTS);   NREM.reward = Restrict(NREM.all,rewardTS);
        REM.baseline = Restrict(REM.all,baselineTS);     REM.aversive = Restrict(REM.all,aversiveTS);     REM.reward = Restrict(REM.all,rewardTS);
        
        %% Load ripples
        if exist('ripplesD_customized2.csv')
            ripplesD = table2array(readtable('ripplesD_customized2.csv'));
            RD = true;
        else
            RD = false;
        end
        
        if exist('ripplesV_customized2.csv')
            ripplesV = table2array(readtable('ripplesV_customized2.csv'));
            RV = true;
        else
            RV = false;
        end
        
        if and(RD,RV)
            RB = true;
            % coordination
            coordinated = [];
            coordinatedV = [];
            cooridnated_event = [];
            for i = 1:length(ripplesD)
                r = ripplesD(i,:);
                tmp = sum(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1));
                if tmp>0
                    z = ripplesV(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1),:);
                    coordinatedV = [coordinatedV ; z];
                    [p,indice] = min(abs(r(2)-z(:,2)));
                    coordinated = [coordinated ; r];
                    
                    %                     if r(2)<z(2) % keep only when dorsal happen first
                    % %                         cooridnated_event = [cooridnated_event ; r];
                    %                         coordinated = [coordinated ; r];
                    %                         coordinatedV = [coordinatedV ; z];
                    %                     end
                    
                    peak = min(min(z(indice,2),r(2)))+abs(z(indice,2)-r(2))/2;
                    low = min([r(1) , z(indice,1)]);
                    up = max([r(3) , z(indice,3)]);
                    cooridnated_event = [cooridnated_event ; low , peak , up];
                    
                    clear tmp2 tmp1 p indice z peak low up
                end
                clear r
            end
            clear x tmp i
            
            [C,IA,IC] = unique(coordinatedV(:,1));
            coordinatedV  = coordinatedV(IA,:); clear C IA IC
            
            % Store events time stamps
            % dRipples
            ripples.dHPC.baseline = Restrict(ripplesD , NREM.baseline);
            ripples.dHPC.reward = Restrict(ripplesD , NREM.reward);
            ripples.dHPC.aversive = Restrict(ripplesD , NREM.aversive);
            % vRipples
            ripples.vHPC.baseline = Restrict(ripplesV , NREM.baseline);
            ripples.vHPC.reward = Restrict(ripplesV , NREM.reward);
            ripples.vHPC.aversive = Restrict(ripplesV , NREM.aversive);
            % coordinated dRipples
            ripples.dHPC.coordinated.all = coordinated;
            ripples.dHPC.uncoordinated.all = ripplesD(not(ismember(ripplesD(:,2) , coordinated(:,2))),:);
            ripples.dHPC.coordinated.baseline = Restrict(coordinated , NREM.baseline);
            ripples.dHPC.coordinated.reward = Restrict(coordinated , NREM.reward);
            ripples.dHPC.coordinated.aversive = Restrict(coordinated , NREM.aversive);
            % coordinated vRipples
            ripples.vHPC.coordinated.all = coordinatedV;
            ripples.vHPC.uncoordinated.all = ripplesV(not(ismember(ripplesV(:,2) , coordinatedV(:,2))),:);
            ripples.vHPC.coordinated.baseline = Restrict(coordinatedV , NREM.baseline);
            ripples.vHPC.coordinated.reward = Restrict(coordinatedV , NREM.reward);
            ripples.vHPC.coordinated.aversive = Restrict(coordinatedV , NREM.aversive);
            %coordinated event
            cooridnated_event((cooridnated_event(:,3)-cooridnated_event(:,1)<0.04),:) = [];
            ripple_event.baseline = Restrict(cooridnated_event,baselineTS);
            ripple_event.reward = Restrict(cooridnated_event,rewardTS);
            ripple_event.aversive = Restrict(cooridnated_event,aversiveTS);
            ripple_event.all = cooridnated_event;
        else
            RB = false;
        end
        
        %% Spikes
        % Load Units
        disp('Uploading Spiking activity')
        cd 'Spikesorting'
        [clusters , numberD , numberV , spks , spks_dHPC , spks_vHPC , cellulartype] = load_SU_FM(cd,criteria_type,criteria_fr,aversiveTS_run,rewardTS_run);
        
        if isfile('dorsalventral_assemblies_aversiveVF.mat')
            iterator = true;
            disp('Loading Aversive template')
            [patterns , cond , Thresholded] = load_assemblies(cd , 'dorsalventral_assemblies_aversiveVF.mat', clusters, numberD , 'aversive');
            th = sum(Thresholded.aversive(:,cond.both),2)>=1;
            patterns = patterns.aversive(:,cond.both);
            Thresholded = Thresholded.aversive(:,cond.both);
        else
            th = zeros(size(clusters.all,1),1);
            iterator = false;
        end
        
        if exist('dHPC_responsivness_all.mat')
            load('dHPC_responsivness_all.mat')
            Responssive.dHPC = [];
            for i = 1 : size(clusters.dHPC,1)
                Responssive.dHPC = [Responssive.dHPC ; dHPC_resp.resp_ave(dHPC_resp.id == clusters.dHPC(i))==1];
            end
        end
        
        if exist('vHPC_responsivness_all.mat')
            load('vHPC_responsivness_all.mat')
            Responssive.vHPC = [];
            for i = 1 : size(clusters.vHPC,1)
                Responssive.vHPC = [Responssive.vHPC ; vHPC_resp.resp_ave(vHPC_resp.id == clusters.vHPC(i))==1];
            end
        end
        
        if and(RB , iterator)
            disp('Running GLM')
            % SpikeTrains construction
            limits = [0 segments.Var1(end)/1000];
            events = [];
            [Spikes , bins , Clusters] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, 0.025, limits, events, false, true);
%             clear limits events
            
            a = assembly_activity(patterns , Spikes');
            
            for i = 1 : size(patterns,2)
                % Ids construction
                template = logical([Thresholded(:,i) [Responssive.dHPC ;  Responssive.vHPC]]);
                
                if sum(and(template(:,1) , template(:,2)))>0
                    % Model
                    [output_metrics] = train_test_poisson_glm_leave_one_out2([bins' , Spikes], template , [bins' a(i,:)'], movement.aversive, [ripple_event.aversive(:,1) ripple_event.aversive(:,3)], [numberD , numberD+numberV]);
                    % Saving data  
                    % dHPC
                    if not(isempty(output_metrics.dHPC.shock_member))
                        Pred_shock.dHPC = [Pred_shock.dHPC ; vertcat(output_metrics.dHPC.shock_member.mse_delta) , vertcat(output_metrics.dHPC.shock_member.pearson_delta) , vertcat(output_metrics.dHPC.shock_member.spearman_delta)];
                    end
                    if not(isempty(output_metrics.dHPC.non_member))
                        Pred_nomember.dHPC = [Pred_nomember.dHPC ; vertcat(output_metrics.dHPC.non_member.mse_delta) , vertcat(output_metrics.dHPC.non_member.pearson_delta) , vertcat(output_metrics.dHPC.non_member.spearman_delta)];
                    end                        
                    if not(isempty(output_metrics.dHPC.member_only))
                        Pred_member.dHPC = [Pred_member.dHPC ; vertcat(output_metrics.dHPC.member_only.mse_delta) , vertcat(output_metrics.dHPC.member_only.pearson_delta) , vertcat(output_metrics.dHPC.member_only.spearman_delta)];
                    end
                    % vHPC
                    if not(isempty(output_metrics.vHPC.shock_member))
                        Pred_shock.vHPC = [Pred_shock.vHPC ; vertcat(output_metrics.vHPC.shock_member.mse_delta) , vertcat(output_metrics.vHPC.shock_member.pearson_delta) , vertcat(output_metrics.vHPC.shock_member.spearman_delta)];
                    end
                    if not(isempty(output_metrics.vHPC.non_member))
                        Pred_nomember.vHPC = [Pred_nomember.vHPC ; vertcat(output_metrics.vHPC.non_member.mse_delta) , vertcat(output_metrics.vHPC.non_member.pearson_delta) , vertcat(output_metrics.vHPC.non_member.spearman_delta)];
                    end                        
                    if not(isempty(output_metrics.vHPC.member_only))
                        Pred_member.vHPC = [Pred_member.vHPC ; vertcat(output_metrics.vHPC.member_only.mse_delta) , vertcat(output_metrics.vHPC.member_only.pearson_delta) , vertcat(output_metrics.vHPC.member_only.spearman_delta)];
                    end
                end
                
                clear output_metrics template
            end
        end
        disp(' ')
        clear A aversiveTS aversiveTS_run baselineTS rewardTS rewardTS_run
        clear behavior bins Cell_type_classification cellulartype cond
        clear is K Kinfo group_dHPC group_vHPC matrixC matrixD matrixV
        clear NREM REM WAKE segmentation tmp cond config
        clear spiketrains_dHPC spiketrains_vHPC opts MUA
        clear patterns Thresholded i  ii numberD numberV movement cross crossN
        clear Spikes bins Clusters Shocks_filt Rewards_filt config n_SU_D n_SU_V
        clear clusters coordinated coordinated_ripple_bursts coordinatedV
        clear cooridnated_event coordinatedV_refined coordinatedV_refined
        clear ripple_bursts ripple_event ripplesD ripplesV
        clear spks spks_dHPC spks_vHPC ripples cooridnated_event
        clear cooridnated_eventDV cooridnated_eventVD segments movement RD RV RB
    end
% 
% 
% --- dHPC ---
y1 = Pred_nomember.dHPC(:,2); y1 = y1(~isnan(y1));
y2 = Pred_member.dHPC(:,2);   y2 = y2(~isnan(y2));
y3 = Pred_shock.dHPC(:,2);    y3 = y3(~isnan(y3));

% Outlier removal (IQR method)
y1 = y1(~isoutlier(y1));
y2 = y2(~isoutlier(y2));
y3 = y3(~isoutlier(y3));

% Determine smallest group size after outlier removal
refs = min([length(y1), length(y2), length(y3)]);

% Random downsampling
y1 = y1(randperm(length(y1), refs));
y2 = y2(randperm(length(y2), refs));
y3 = y3(randperm(length(y3), refs));

matrix1 = [ ...
    [ones(refs,1), y1]; ...
    [2*ones(refs,1), y2]; ...
    [3*ones(refs,1), y3] ...
    ];

% --- vHPC ---
y1 = Pred_nomember.vHPC(:,2); y1 = y1(~isnan(y1));
y2 = Pred_member.vHPC(:,2);   y2 = y2(~isnan(y2));
y3 = Pred_shock.vHPC(:,2);    y3 = y3(~isnan(y3));

% Outlier removal
y1 = y1(~isoutlier(y1));
y2 = y2(~isoutlier(y2));
y3 = y3(~isoutlier(y3));

% Determine new minimum
refs = min([length(y1), length(y2), length(y3)]);

% Random downsampling
y1 = y1(randperm(length(y1), refs));
y2 = y2(randperm(length(y2), refs));
y3 = y3(randperm(length(y3), refs));

matrix2 = [ ...
    [ones(refs,1), y1]; ...
    [2*ones(refs,1), y2]; ...
    [3*ones(refs,1), y3] ...
    ];
two_way_anova_with_posthoc(matrix1, matrix2)
yline(0,'--')
% 
% %% -------%%
% y_data1 = [Pred_nomember.dHPC(:,2) ; Pred_member.dHPC(:,2) ; Pred_shock.dHPC(:,2)];
% x_data1 = [ones(length(Pred_nomember.dHPC(:,1)),1) ; ones(length(Pred_member.dHPC(:,1)),1)*2 ; ones(length(Pred_shock.dHPC(:,1)),1)*3];
% data1 = [x_data1 , y_data1];
% 
% y = Pred_member.vHPC(randperm(length(Pred_member.vHPC)),2);
% y = y(1:66,1);
% y_data2 = [Pred_nomember.vHPC(:,2) ; y ; Pred_shock.vHPC(:,2)];
% x_data2 = [ones(length(Pred_nomember.vHPC(:,1)),1) ; ones(length(y(:,1)),1)*2 ; ones(length(Pred_shock.vHPC(:,1)),1)*3];
% data2 = [x_data2 , y_data2];
% 
% two_way_anova_with_posthoc(data1, data2)
% 
% % 
end