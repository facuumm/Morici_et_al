function [Pre , Post , T , I , curve] = Shock_responsive_ripples(path) 
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
dur = 1;
bin = 0.01;

% storage variables
Pre.aversive.dHPC = [];      Post.aversive.dHPC = [];
Pre.reward.dHPC = [];        Post.reward.dHPC = [];


Pre.aversive.vHPC = [];      Post.aversive.vHPC = [];
Pre.reward.vHPC = [];        Post.reward.vHPC = [];

I.dHPC = [];                 I.vHPC = [];

curve.dHPC = [];             curve.vHPC = [];
curve.id.dHPC = [];          curve.id.vHPC = [];
curve.speed.data = [];       curve.speed.id = [];

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
        for i = 1 : length(Shocks_filt)
            [~ , ii] = min(abs(means(:,1)-Shocks_filt(i)));
            if ii+20 < length(means)
                tmp = [tmp , means(ii-20 : ii+20 , 2)];
            end
        end
        Mean = nanmean(tmp'); clear tmp means
        curve.speed.data = [curve.speed.data , Mean']; clear Mean
        curve.speed.id = [curve.speed.id ; tt t length(Shocks_filt)];
        
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
            
            load([cd , '\coordinated_ripple_bursts.mat'])
            
            % Separation of bursts across emotional condition
            bursts.coordinated.DV(bursts.coordinated.DV(:,2)-bursts.coordinated.DV(:,1)<0.06, : ) = [];
            %             bursts.coordinated.DV(bursts.coordinated.DV(:,2)-bursts.coordinated.DV(:,1)>0.3, : ) = [];
            
            bursts.baseline = Restrict(bursts.coordinated.DV,NREM.baseline);
            bursts.aversive = Restrict(bursts.coordinated.DV,NREM.aversive);
            bursts.reward = Restrict(bursts.coordinated.DV,NREM.reward);
        else
            RB = false;
        end
        
        %% Spikes
        % Load Units
        disp('Uploading Spiking activity')
        cd 'Spikesorting'
        [clusters , numberD , numberV , spks , spks_dHPC , spks_vHPC , cellulartype] = load_SU_FM(cd,criteria_type,criteria_fr,aversiveTS_run,rewardTS_run);
        
        if isfile('dorsalventral_assemblies_aversiveVF.mat')
            disp('Loading Aversive template')
            [patterns , cond , Thresholded] = load_assemblies(cd , 'dorsalventral_assemblies_aversiveVF.mat', clusters, numberD , 'aversive');
            th = sum(Thresholded.aversive(:,cond.both),2)>=1;
        else
            th = zeros(size(clusters.all,1),1);
        end
        
        if exist('dHPC_responsivness_all_VF.mat')
            load('dHPC_responsivness_all_VF.mat')
            curve.dHPC = [curve.dHPC ; dHPC_resp.curve_ave];
            x = ones(size(dHPC_resp.id,1),1);
%             curve.id.dHPC = [curve.id.dHPC ; dHPC_resp.id , dHPC_resp.resp_ave x*tt x*t dHPC_resp.Speedcorrelation.all(:,2)];
            curve.id.dHPC = [curve.id.dHPC ; dHPC_resp.id , dHPC_resp.resp_ave x*tt x*t];

            if RB
                for i = 1 : 2
                    if i == 1
                        if aversiveTS_run(1)<rewardTS_run(1)
                            TS.pre = NREM.baseline;       TS.post = NREM.aversive;
                        else
                            TS.pre = NREM.reward;         TS.post = NREM.aversive;
                        end
                    else
                        if aversiveTS_run(1)<rewardTS_run(1)
                            TS.pre = NREM.aversive;        TS.post = NREM.reward;                    
                        else
                            TS.pre = NREM.baseline;        TS.post = NREM.reward;
                        end
                    end
                    
                    for ii = 1 : numberD
                        % Pre
                        y = Restrict(spks(spks(:,1)==clusters.dHPC(ii),2),TS.pre);
                        x = Restrict(ripples.dHPC.uncoordinated.all(:,1:3),TS.pre);
%                         x = Restrict(ripplesD(:,1:3),TS.pre);
                        
                        if and(size(x,1)>5 , size(y,1)>5)
                            x = Restrict(x(:,2),[y(1,1) y(end,1)]);
                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                            [ccg1,T] = CCG(s,ids,'binSize',bin,'duration',dur,'smooth',sm,'mode','ccg');
                            ccg1 = ccg1(:,1,2)./size(x,1);    ccg1 = ccg1./bin; %ccg1 = ccg1./m;
                        else
                            ccg1 = [];
                        end
                        clear x y s ids groups m
                        
                        % Post
                        y = Restrict(spks(spks(:,1)==clusters.dHPC(ii),2),TS.post);
                        x = Restrict(ripples.dHPC.uncoordinated.all(:,1:3),TS.post);

                        if and(size(x,1)>5 , size(y,1)>5)
                            x = Restrict(x(:,2),[y(1,1) y(end,1)]);
                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                            [ccg2,T] = CCG(s,ids,'binSize',bin,'duration',dur,'smooth',sm,'mode','ccg');
                            ccg2 = ccg2(:,1,2)./size(x,1);    ccg2 = ccg2./bin; %ccg1 = ccg1./m;
                        else
                            ccg2 = [];
                        end
                        clear x y s ids groups m
                        
                        if and(not(isempty(ccg1)) ,  not(isempty(ccg2)))
                            if i == 1
                                Pre.aversive.dHPC = [Pre.aversive.dHPC , ccg1];
                                Post.aversive.dHPC = [Post.aversive.dHPC , ccg2];
                                
                                I.dHPC = [I.dHPC ; dHPC_resp.resp_ave(dHPC_resp.id == clusters.dHPC(ii)) clusters.dHPC(ii) t tt th(ii)]; clear p
                            else
                                Pre.reward.dHPC = [Pre.reward.dHPC , ccg1];
                                Post.reward.dHPC = [Post.reward.dHPC , ccg2];
                            end
                        end
                        clear ccg1 ccg2
                    end
                    clear TS
                end
                clear i ii
            end
        end
        
        if exist('vHPC_responsivness_all_VF.mat')
            load('vHPC_responsivness_all_VF.mat')
            curve.vHPC = [curve.vHPC ; vHPC_resp.curve_ave];
            x = ones(size(vHPC_resp.id,1),1);
%             curve.id.vHPC = [curve.id.vHPC ; vHPC_resp.id , vHPC_resp.resp_ave x*tt x*t vHPC_resp.Speedcorrelation.all(:,2)];
            curve.id.vHPC = [curve.id.vHPC ; vHPC_resp.id , vHPC_resp.resp_ave x*tt x*t];

            if RB
                for i = 1 : 2
                    if i == 1
                        if aversiveTS_run(1)<rewardTS_run(1)
                            TS.pre = NREM.baseline;       TS.post = NREM.aversive;
                        else
                            TS.pre = NREM.reward;         TS.post = NREM.aversive;
                        end
                    else
                        if aversiveTS_run(1)<rewardTS_run(1)
                            TS.pre = NREM.aversive;        TS.post = NREM.reward;                       
                        else
                            TS.pre = NREM.baseline;        TS.post = NREM.reward;
                        end
                    end
                    
                    for ii = 1 : numberV
                        % Pre
                        y = Restrict(spks(spks(:,1)==clusters.vHPC(ii),2),TS.pre);
                        x = Restrict(ripples.vHPC.uncoordinated.all(:,1:3),TS.pre);

                        if and(size(x,1)>5 , size(y,1)>5)
                            x = Restrict(x(:,2),[y(1,1) y(end,1)]);
                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                            [ccg1,T] = CCG(s,ids,'binSize',bin,'duration',dur,'smooth',sm,'mode','ccg');
                            ccg1 = ccg1(:,1,2)./size(x,1);    ccg1 = ccg1./bin; %ccg1 = ccg1./m;
                        else
                            ccg1 = [];
                        end
                        clear x y s ids groups
                        
                        % Post
                        y = Restrict(spks(spks(:,1)==clusters.vHPC(ii),2),TS.post);
                        x = Restrict(ripples.vHPC.uncoordinated.all(:,1:3),TS.post);
                        
                        if and(size(x,1)>5 , size(y,1)>5)
                            x = Restrict(x(:,2),[y(1,1) y(end,1)]);
                            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
                            [ccg2,T] = CCG(s,ids,'binSize',bin,'duration',dur,'smooth',sm,'mode','ccg');
                            ccg2 = ccg2(:,1,2)./size(x,1);    ccg2 = ccg2./bin; %ccg1 = ccg1./m;
                        else
                            ccg2 = [];
                        end
                        clear x y s ids groups m
                        
                        if and(not(isempty(ccg1)) ,  not(isempty(ccg2)))
                            if i == 1
                                Pre.aversive.vHPC = [Pre.aversive.vHPC , ccg1];
                                Post.aversive.vHPC = [Post.aversive.vHPC , ccg2];
                                
                                I.vHPC = [I.vHPC ; vHPC_resp.resp_ave(vHPC_resp.id == clusters.vHPC(ii)) clusters.vHPC(ii) t tt th(ii+numberD)]; clear p
                            else
                                Pre.reward.vHPC = [Pre.reward.vHPC , ccg1];
                                Post.reward.vHPC = [Post.reward.vHPC , ccg2];
                            end
                        end
                        clear ccg1 ccg2
                    end
                    clear TS
                end
                clear i ii
            end
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
end
% %% To plot data as it is in the paper
% %% Sorting first shock responssive
% t = [-2:0.1:2];
% % indices to calculate the mean and sort
% [~ , low] = min(abs(t-0));
% [~ , up] = min(abs(t-1));
% 
% % --- dHPC ---
% % Separate the row with IDX == 1 from the rest of the matrix
% IDX = curve.id.dHPC(:,2);
% row_to_move = curve.dHPC(IDX == 1, :);  % Row where IDX == 1
% remaining_rows = curve.dHPC(IDX == -1, :);  % The remaining rows where IDX == 0
% 
% % Sort both the row with IDX == 1 and the remaining rows by the mean of each row
% [~, sort_idx_remaining] = sort(nanmean(remaining_rows(:,low:up)'));  % Sort remaining rows by mean
% sorted_remaining_rows = remaining_rows(sort_idx_remaining, :);  % Sorted remaining rows
% 
% [~, sort_idx_move] = sort(nanmean(row_to_move(:,low:up)'));  % Sort row with IDX == 1 by mean (although only one row)
% sorted_row_to_move = row_to_move(sort_idx_move, :);  % Sorted row with IDX == 1
% 
% % Reconstruct the matrix: first the sorted row with IDX == 1, then the sorted remaining rows
% reordered_matrix = [sorted_remaining_rows;sorted_row_to_move];
% IDX = logical([zeros(1,length(sort_idx_remaining)),ones(1,length(sort_idx_move))]);
% 
% figure
% subplot(121)
% imagesc(t,[1:size(reordered_matrix,1)],reordered_matrix),colormap 'gray',clim([0 1])
% set(gca, 'YDir','normal');
% hold on, xline(0), xline(1)
% 
% subplot(122)
% hold on;  % Keep the current plot
% % Mark the lines (rows in this case) according to the logical vector
% for i = 1:length(IDX)
%     if IDX(i)  % Check if the logical vector value is true
%         plot([1 size(curve.dHPC,2)], [i i], 'r', 'LineWidth', 0.0001);  % Plot a red line on the row
%     end
% end
% ylim([1, length(IDX)])
% hold off;  % Release the plot hold
% 
% zeros(1,length(sort_idx_remaining))
% 
% 
% 
% % --- vHPC ---
% % indices to calculate the mean and sort
% [~ , low] = min(abs(t-0));
% [~ , up] = min(abs(t-1));
% 
% % Separate the row with IDX == 1 from the rest of the matrix
% IDX = curve.id.vHPC(:,2);
% row_to_move = curve.vHPC(IDX == 1, :);  % Row where IDX == 1
% remaining_rows = curve.vHPC(IDX == -1, :);  % The remaining rows where IDX == 0
% 
% % Sort both the row with IDX == 1 and the remaining rows by the mean of each row
% [~, sort_idx_remaining] = sort(nanmean(remaining_rows(:,low:up)'));  % Sort remaining rows by mean
% sorted_remaining_rows = remaining_rows(sort_idx_remaining, :);  % Sorted remaining rows
% 
% [~, sort_idx_move] = sort(nanmean(row_to_move(:,low:up)'));  % Sort row with IDX == 1 by mean (although only one row)
% sorted_row_to_move = row_to_move(sort_idx_move, :);  % Sorted row with IDX == 1
% 
% % Reconstruct the matrix: first the sorted row with IDX == 1, then the sorted remaining rows
% reordered_matrix = [sorted_remaining_rows;sorted_row_to_move];
% IDX = logical([zeros(1,length(sort_idx_remaining)),ones(1,length(sort_idx_move))]);
% 
% figure
% subplot(121)
% imagesc(t,[1:size(reordered_matrix,1)],reordered_matrix),colormap 'gray',clim([0 1])
% set(gca, 'YDir','normal');
% hold on, xline(0), xline(1)
% 
% subplot(122)
% hold on;  % Keep the current plot
% % Mark the lines (rows in this case) according to the logical vector
% for i = 1:length(IDX)
%     if IDX(i)  % Check if the logical vector value is true
%         plot([1 size(curve.dHPC,2)], [i i], 'r', 'LineWidth', 0.0001);  % Plot a red line on the row
%     end
% end
% ylim([1, length(IDX)])
% hold off;  % Release the plot hold
% 
% zeros(1,length(sort_idx_remaining))