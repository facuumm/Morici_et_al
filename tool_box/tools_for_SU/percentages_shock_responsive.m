function [dHPC vHPC] = percentages_shock_responsive(path)
% This function calculates the Pre and Post sleep assemblies rate.
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% OUTPUT
% It will plot the percentage of dHPC and vHPC shock responssive SU.
% Down, it will show you the percentage discriminated by their
% participation in assemblies.
%
% Morci Juan Facundo 04/2024

criteria_fr = 0;
criteria_type = 0; % criteria for celltype (0:pyr, 1:int, 2:all)

dHPC = [];
vHPC = [];

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
        load('behavioral_data.mat')
        baselineTS = baselineTS./1000;
        aversiveTS = aversiveTS./1000;
        rewardTS = rewardTS./1000;
        aversiveTS_run = aversiveTS_run./1000;
        rewardTS_run = rewardTS_run./1000;
        
        %% Spikes
        % Load Units
        disp('Uploading Spiking activity')
        cd 'Spikesorting'
        [clusters , numberD , numberV , spks , spks_dHPC , spks_vHPC , cellulartype] = load_SU_FM(cd,criteria_type,criteria_fr,aversiveTS_run,rewardTS_run);
        
        %% Assemblies detection
%         if and(numberD > 3 , numberV > 3)
            
            % Load data regarding the shock
            load('vHPC_responsivness_all.mat')
            load('dHPC_responsivness_all.mat')
            
            % --- Aversive ---
            disp('Lets go for the assemblies')
            if isfile('dorsalventral_assemblies_aversiveVF.mat')
                disp('Loading Aversive template')
                load('dorsalventral_assemblies_aversiveVF.mat')
                
                Thresholded.aversive = Th;    patterns.aversive = pat;     clear cond Th pat
                patterns.aversive = patterns.aversive .* Thresholded.aversive;
                
                cond  = classification_of_asselblies(Thresholded.aversive,clusters.dHPC);% Detection of members
                
                i = sum(Thresholded.aversive(1:numberD,cond.dHPC),2)>0;
                ii = sum(Thresholded.aversive(1:numberD,cond.both),2)>0;
                iii = not(or(sum(Thresholded.aversive(1:numberD,cond.dHPC),2)>0 , sum(Thresholded.aversive(1:numberD,cond.both),2)>0));
                dHPC = [dHPC ; dHPC_resp.resp_ave , i , ii , iii]; clear i ii iii
                
                i = sum(Thresholded.aversive(numberD+1:end,cond.vHPC),2)>0;
                ii = sum(Thresholded.aversive(numberD+1:end,cond.both),2)>0;
                iii = not(or(sum(Thresholded.aversive(numberD+1:end,cond.vHPC),2)>0 , sum(Thresholded.aversive(numberD+1:end,cond.both),2)>0));
                vHPC = [vHPC ; vHPC_resp.resp_ave , i , ii , iii]; clear i ii iii
            end
            
%         end
        
        
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
    clear cooridnated_eventDV cooridnated_eventVD segments movement
    clear vHPC_resp dHPC_resp
end

% Percentage calculation for dHPC SUs
percentage.dHPC.responsive = [(sum(dHPC(:,1)==1)/length(dHPC(:,1)))*100];
percentage.vHPC.responsive = [(sum(vHPC(:,1)==1)/length(vHPC(:,1)))*100];

i    = sum(and(dHPC(:,1)==1 , dHPC(:,2)==1));
ii   = sum(and(dHPC(:,1)==1 , dHPC(:,3)==1));
iii  = sum(and(dHPC(:,1)==1 , dHPC(:,4)==1));
iiii = i+ii+iii;
percentages.dHPC.dHPC = (i/iiii)*100;
percentages.dHPC.joint = (ii/iiii)*100;
percentages.dHPC.non = (iii/iiii)*100;

i    = sum(and(vHPC(:,1)==1 , vHPC(:,2)==1));
ii   = sum(and(vHPC(:,1)==1 , vHPC(:,3)==1));
iii  = sum(and(vHPC(:,1)==1 , vHPC(:,4)==1));
iiii = i+ii+iii;
percentages.vHPC.vHPC = (i/iiii)*100;
percentages.vHPC.joint = (ii/iiii)*100;
percentages.vHPC.non = (iii/iiii)*100;

figure,
subplot(221),pie([percentage.dHPC.responsive , 100-percentage.dHPC.responsive] , {'Responssive' , 'No-Responssive'})
subplot(222),pie([percentage.vHPC.responsive , 100-percentage.vHPC.responsive] , {'Responssive' , 'No-Responssive'})
subplot(223),pie([percentages.dHPC.dHPC percentages.dHPC.joint percentages.dHPC.non],{'dHPC' , 'Joint' , 'NoMember'})
subplot(224),pie([percentages.vHPC.vHPC percentages.vHPC.joint percentages.vHPC.non],{'vHPC' , 'Joint' , 'NoMember'})
end
