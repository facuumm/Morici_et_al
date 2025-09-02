function [dorsal , ventral , extra] = ccg_shock_quiet(path)
% This function iterates inside the path uploading the Neurons, their
% responsivnes to shocks and lock their activity to transitions from
% quiet to movement.
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% OUTPUT
% dorsal, ventral: strucutre that contains the id, shock_responsivnes, curve
%
% extra: structure containing time vector for plotting curves, average
%        speed sourrounding quiet-to-active transitions, and # of
%        transitions.
%
% Morci Juan Facundo 09/2025

% Variables to use in the script
criteria_fr = 0;
criteria_type = 0; % criteria for celltype (0:pyr, 1:int, 2:all)

b = 0.1;
d = 10;
sm = 1;

% Storage variables
dorsal.id = [];    ventral.id = [];
dorsal.curve = [];    ventral.curve = [];
extra.time = [];
extra.speed = [];
extra.transitions = [];

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
        
        %% Load Session organization
        disp('Loading TimeStamps Session Structure')
        load('session_organization.mat')
        
        %% Load Behavioral data
        disp('Loading Behvioral Data')
        load('behavioral_dataVF.mat')
        
        %% speed sourrounding the shocks
        means = meanInGroups(behavior.speed.aversive(:,2), 3);
        downsampled_t = downsampleTimeVector(behavior.speed.aversive(:,1), 1/30, 0.1);
        
        if size(means,2) < size(downsampled_t,2)
            means = [downsampled_t(1:size(means,2))' , means']; clear downsampled_t
        elseif size(means,2) > size(downsampled_t,2)
            means = [downsampled_t' , means(1:size(downsampled_t,2))']; clear downsampled_t
        else
            means = [downsampled_t' , means']; clear downsampled_t
        end
        
        tmp = [];
        for i = 1 : length(behavior.quiet.aversive(:,2))
            [~ , ii] = min(abs(means(:,1)-behavior.quiet.aversive(i,2)));
            if and(ii+50 < length(means) , ii-50>=1)
                tmp = [tmp , means(ii-50 : ii+50 , 2)];
            end
        end
        Mean = nanmean(tmp'); clear tmp means
        extra.speed = [extra.speed , Mean']; clear Mean      
        extra.transitions = [extra.transitions ; length(behavior.quiet.aversive(:,2))];
        
        %% Spikes
        cd 'Spikesorting'
        % Load Units
        disp('Uploading Spiking activity')        
        [clusters , numberD , numberV , spks , spks_dHPC , spks_vHPC , cellulartype] = load_SU_FM(cd,criteria_type,criteria_fr,aversiveTS_run./1000,rewardTS_run./1000);
        
        %% Shock
        if exist('dHPC_Shock_VF.mat')
            load('dHPC_Shock_VF.mat')
            for i = 1 : size(dHPC_resp.id,1)
                cluster = dHPC_resp.id(i);
                resp = dHPC_resp.resp_ave(i);
                
                Times1 = behavior.quiet.aversive(:,2);
                Times2 = spks(spks(:,1)==cluster,2);
                
                [s,ids,groups] = CCGParameters(Times1,ones(length(Times1),1),Times2,ones(length(Times2),1)*2);
                [ccg,T] = CCG(s,ids,'binSize',b,'duration',d,'smooth',sm,'mode','ccg');
                ccg = (ccg(:,1,2)./length(Times1))./b;
                
                dorsal.id = [dorsal.id ; cluster , resp];
                dorsal.curve = [dorsal.curve , zscore(ccg)];
                extra.time = T;
                
                clear ccg s ids groups T Times1 Times2
            end
        end
        
        if exist('vHPC_Shock_VF.mat')
            load('vHPC_Shock_VF.mat')
            for i = 1 : size(vHPC_resp.id,1)
                cluster = vHPC_resp.id(i);
                resp = vHPC_resp.resp_ave(i);
                
                Times1 = behavior.quiet.aversive(:,2);
                Times2 = spks(spks(:,1)==cluster,2);
                
                [s,ids,groups] = CCGParameters(Times1,ones(length(Times1),1),Times2,ones(length(Times2),1)*2);
                [ccg,T] = CCG(s,ids,'binSize',b,'duration',d,'smooth',sm,'mode','ccg');
                ccg = (ccg(:,1,2)./length(Times1))./b;
                
                ventral.id = [ventral.id ; cluster , resp];
                ventral.curve = [ventral.curve , zscore(ccg)];
                extra.time = T;
                
                clear ccg s ids groups T Times1 Times2                
            end
        end
    end
end

figure,
i =dorsal.id(:,2)==1;
M = nanmean(zscore(dorsal.curve,0)');
S = nansem(zscore(dorsal.curve,0)');

plot(extra.time,M,'k'),hold on
ciplot(M-S , M+S , extra.time,'k'), alpha 0.2


i =ventral.id(:,2)==1;
M = nanmean(zscore(ventral.curve,0)');
S = nansem(zscore(ventral.curve,0)');
plot(extra.time,M,'b'),hold on
ciplot(M-S , M+S , extra.time,'b'), alpha 0.2

figure
plot(extra.time,nanmean(extra.speed')),
hold on
ciplot(nanmean(extra.speed') - nansem(extra.speed') , nanmean(extra.speed') + nansem(extra.speed') , extra.time), alpha 0.2
end