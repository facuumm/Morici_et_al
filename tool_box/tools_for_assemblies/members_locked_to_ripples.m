function [Pre Post] = members_locked_to_ripples(path)
% This function calculates the Pre and Post sleep assemblies rate.
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% OUTPUT
% Pre, Post: Structure, it store the Assemblies Rate for each type
%            of assemblies.
%
%               Architecture of each output:
%                   Pre.aversive
%                      .reward
%                   Post.aversive
%                       .reward
%
%
% Morci Juan Facundo 07/2024

% variables to use in the script
criteria_fr = 0;
criteria_type = 0; % criteria for celltype (0:pyr, 1:int, 2:all)
th = 5; % threshold for detecting peak assemblies
% 3 SD funciona
% 5 SD funciona
sm = 2;
dur = 0.6;

% storage variables
Pre.aversive.cooridnated.dRipples.dHPC = [];
Pre.aversive.cooridnated.dRipples.vHPC = [];
Pre.aversive.cooridnated.vRipples.dHPC = [];
Pre.aversive.cooridnated.vRipples.vHPC = [];
Pre.aversive.uncooridnated.dRipples.dHPC = [];
Pre.aversive.uncooridnated.dRipples.vHPC = [];
Pre.aversive.uncooridnated.vRipples.dHPC = [];
Pre.aversive.uncooridnated.vRipples.vHPC = [];

Post.aversive.cooridnated.dRipples.dHPC = [];
Post.aversive.cooridnated.dRipples.vHPC = [];
Post.aversive.cooridnated.vRipples.dHPC = [];
Post.aversive.cooridnated.vRipples.vHPC = [];
Post.aversive.uncooridnated.dRipples.dHPC = [];
Post.aversive.uncooridnated.dRipples.vHPC = [];
Post.aversive.uncooridnated.vRipples.dHPC = [];
Post.aversive.uncooridnated.vRipples.vHPC = [];


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
        
        %% load sleep states
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);
        clear x states
        
        baselineTS = baselineTS./1000;
        aversiveTS = aversiveTS./1000;
        rewardTS = rewardTS./1000;
        aversiveTS_run = aversiveTS_run./1000;
        rewardTS_run = rewardTS_run./1000;
        
        NREM.baseline = Restrict(NREM.all,baselineTS);
        NREM.aversive = Restrict(NREM.all,aversiveTS);
        NREM.reward = Restrict(NREM.all,rewardTS);
        
        REM.baseline = Restrict(REM.all,baselineTS);
        REM.aversive = Restrict(REM.all,aversiveTS);
        REM.reward = Restrict(REM.all,rewardTS);
        
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
                    
                    cooridnated_event = [cooridnated_event ; min([r(1) , z(indice,1)]) , min(min(z(indice,2),r(2)))+abs(z(indice,2)-r(2))/2 , max([r(3) , z(indice,3)])];
                    
                    clear tmp2 tmp1 p indice z
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
            bursts.baseline = Restrict(bursts.coordinated,NREM.baseline);
            bursts.aversive = Restrict(bursts.coordinated,NREM.aversive);
            bursts.reward = Restrict(bursts.coordinated,NREM.reward);
        end
        
        %% Spikes
        % Load Units
        disp('Uploading Spiking activity')
        cd 'Spikesorting'
        spks = double([readNPY('spike_clusters.npy') readNPY('spike_times.npy')]);
        K = tdfread('cluster_group.tsv'); % import clusters ID and groups (MUA,Noise,Good)
        Kinfo = tdfread('cluster_info.tsv'); % import information of clusters
        K = [K.cluster_id(K.group(:,1) == 'g') , Kinfo.ch(K.group(:,1) == 'g'),]; % defining good clusters
        % Load neuronal classification
        load('Cell_type_classification')
        K = [K , Cell_type_classification(:,6:8)];
        group_dHPC = K(K(:,2) > 63,:);
        group_vHPC = K(K(:,2) <= 63,:);
        
        %Loop to select dorsal or ventral LFP and SU
        % z=1 --> dorsal
        % z=2 --> ventral
        for z = 1:2
            if z == 1
                spks_dHPC = spks(ismember(spks(:,1),group_dHPC(:,1)),:); %keep spks from good clusters
            else
                spks_vHPC = spks(ismember(spks(:,1),group_vHPC(:,1)),:); %keep spks from good clusters
            end
        end
        clear z
        spks_vHPC(:,2) = double(spks_vHPC(:,2))./20000;
        spks_dHPC(:,2) = double(spks_dHPC(:,2))./20000;
        spks(:,2) = double(spks(:,2))./20000;
        
        % Selection of celltype to analyze
        if criteria_type == 0 %pyr
            cellulartype = [K(:,1) , K(:,4)];
        elseif criteria_type == 1 % int
            cellulartype = [K(:,1) , not(K(:,4))];
        elseif criteria_type == 2 % all
            cellulartype = [K(:,1) , ones(length(K),1)];
        end
        
        %% Counting the Number f SU
        numberD = 0;
        clusters.all = [];
        clusters.dHPC = [];
        for ii=1:size(group_dHPC,1)
            cluster = group_dHPC(ii,1);
            Cellulartype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if Cellulartype
                a = length(Restrict(spks(spks(:,1)==cluster,2),aversiveTS_run)) / ((aversiveTS_run(2)-aversiveTS_run(1)));
                r = length(Restrict(spks(spks(:,1)==cluster,2),rewardTS_run)) / ((rewardTS_run(2)-rewardTS_run(1)));
                if or(a > criteria_fr ,  r > criteria_fr)
                    numberD = numberD+1;
                    clusters.all = [clusters.all ; cluster];
                    clusters.dHPC = [clusters.dHPC ; cluster];
                end
            end
            clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5 r a
        end
        
        numberV = 0;
        clusters.vHPC = [];
        for ii=1:size(group_vHPC,1)
            cluster = group_vHPC(ii,1);
            Cellulartype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if Cellulartype
                a = length(Restrict(spks(spks(:,1)==cluster,2),aversiveTS_run)) / ((aversiveTS_run(2)-aversiveTS_run(1)));
                r = length(Restrict(spks(spks(:,1)==cluster,2),rewardTS_run)) / ((rewardTS_run(2)-rewardTS_run(1)));
                if or(a > criteria_fr ,  r > criteria_fr)
                    numberV = numberV+1;
                    clusters.all = [clusters.all ; cluster];
                    clusters.vHPC = [clusters.vHPC ; cluster];
                end
            end
            clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5 r a
        end
        clear freq limits
        clear camara shock rightvalve leftvalve
        clear ejeX ejeY dX dY dX_int dY_int
        
        %% Assemblies detection
        if and(numberD > 3 , numberV > 3)
            
            limits = [0 segments.Var1(end)/1000];
            events = [];
            [Spikes , bins , Clusters] = spike_train_construction([spks_dHPC;spks_vHPC], clusters.all, cellulartype, 0.025, limits, events, false, true);
            clear limits events
            
            % --- Aversive ---
            disp('Lets go for the assemblies')
            if isfile('dorsalventral_assemblies_aversiveVF.mat')
                disp('Loading Aversive template')
                load('dorsalventral_assemblies_aversiveVF.mat')
                
                Thresholded.aversive = Th;
                patterns.aversive = pat;
                clear cond Th pat
                
                patterns.aversive = patterns.aversive .* Thresholded.aversive;
                
                % Detection of members
                if not(isempty(patterns.aversive))
                    if numberD>0
                        cond1 =  sum(Thresholded.aversive(1:size(clusters.dHPC,1),:),1)>0; %checking of dHPC SU
                        cond2 =  sum(Thresholded.aversive(size(clusters.dHPC,1)+1:end,:),1)>0; %checking of vHPC SU
                        cond.both = and(cond1 , cond2); clear cond1 cond2
                    else
                        cond1 =  logical(zeros(1,size(Thresholded.aversive,2))); %checking of dHPC SU
                        cond2 =  logical(ones(1,size(Thresholded.aversive,2))); %checking of vHPC SU
                        cond.both = and(cond1 , cond2); clear cond1 cond2
                    end
                else
                    cond1 =  false; %checking of dHPC SU
                    cond2 =  logical(0); %checking of vHPC SU
                    cond.both = and(cond1 , cond2); clear cond1 cond2
                end
                
                if sum(cond.both)>0
                    
                    members = Thresholded.aversive(:,cond.both);
                    
                    if aversiveTS_run(1)<rewardTS_run(1)
                        TS.pre = NREM.baseline;
                        TS.post = NREM.aversive;
                    else
                        TS.pre = NREM.reward;
                        TS.post = NREM.aversive;
                    end
                    
                    for i = 1 : sum(cond.both)
                        % coordinated dRipples
                        % dMembers
                        for ii = 1 : numberD
                            if members(ii,i)
                                member = clusters.dHPC(ii);
                                x = spks_dHPC(spks_dHPC(:,1) == member , 2);
                                [pInc1 pDec1 surp1] = RippleModulation_assemblies(ripples.dHPC.coordinated.all,x,TS.pre);
                                [pInc2 pDec2 surp2] = RippleModulation_assemblies(ripples.dHPC.coordinated.all,x,TS.post);
                                Pre.aversive.cooridnated.dRipples.dHPC = [Pre.aversive.cooridnated.dRipples.dHPC ; pInc1 pDec1];
                                Post.aversive.cooridnated.dRipples.dHPC = [Post.aversive.cooridnated.dRipples.dHPC ; pInc2 pDec2];
                                clear pInc1 pDec1 surp1 pInc2 pDec2 surp2
                            end
                        end
                        % vMembers
                        for ii = 1 : numberV
                            if members(ii+numberD,i)
                                member = clusters.vHPC(ii);
                                x = spks_vHPC(spks_vHPC(:,1) == member , 2);
                                [pInc1 pDec1 surp1] = RippleModulation_assemblies(ripples.dHPC.coordinated.all,x,TS.pre);
                                [pInc2 pDec2 surp2] = RippleModulation_assemblies(ripples.dHPC.coordinated.all,x,TS.post);
                                Pre.aversive.cooridnated.dRipples.vHPC = [Pre.aversive.cooridnated.dRipples.vHPC ; pInc1 pDec1];
                                Post.aversive.cooridnated.dRipples.vHPC = [Post.aversive.cooridnated.dRipples.vHPC ; pInc2 pDec2];
                                clear pInc1 pDec1 surp1 pInc2 pDec2 surp2
                            end
                        end
                        
                        % uncoordinated dRipples
                        % dMembers
                        for ii = 1 : numberD
                            if members(ii,i)
                                member = clusters.dHPC(ii);
                                x = spks_dHPC(spks_dHPC(:,1) == member , 2);
                                [pInc1 pDec1 surp1] = RippleModulation_assemblies(ripples.dHPC.uncoordinated.all,x,TS.pre);
                                [pInc2 pDec2 surp2] = RippleModulation_assemblies(ripples.dHPC.uncoordinated.all,x,TS.post);
                                Pre.aversive.uncooridnated.dRipples.dHPC = [Pre.aversive.uncooridnated.dRipples.dHPC ; pInc1 pDec1];
                                Post.aversive.uncooridnated.dRipples.dHPC = [Post.aversive.uncooridnated.dRipples.dHPC ; pInc2 pDec2];
                                clear pInc1 pDec1 surp1 pInc2 pDec2 surp2
                            end
                        end
                        % vMembers
                        for ii = 1 : numberV
                            if members(ii+numberD,i)
                                member = clusters.vHPC(ii);
                                x = spks_vHPC(spks_vHPC(:,1) == member , 2);
                                [pInc1 pDec1 surp1] = RippleModulation_assemblies(ripples.dHPC.uncoordinated.all,x,TS.pre);
                                [pInc2 pDec2 surp2] = RippleModulation_assemblies(ripples.dHPC.uncoordinated.all,x,TS.post);
                                Pre.aversive.uncooridnated.dRipples.vHPC = [Pre.aversive.uncooridnated.dRipples.vHPC ; pInc1 pDec1];
                                Post.aversive.uncooridnated.dRipples.vHPC = [Post.aversive.uncooridnated.dRipples.vHPC ; pInc2 pDec2];
                                clear pInc1 pDec1 surp1 pInc2 pDec2 surp2
                            end
                        end
                        
                        % coordinated vRipples
                        % dMembers
                        for ii = 1 : numberD
                            if members(ii,i)
                                member = clusters.dHPC(ii);
                                x = spks_dHPC(spks_dHPC(:,1) == member , 2);
                                [pInc1 pDec1 surp1] = RippleModulation_assemblies(ripples.vHPC.coordinated.all,x,TS.pre);
                                [pInc2 pDec2 surp2] = RippleModulation_assemblies(ripples.vHPC.coordinated.all,x,TS.post);
                                Pre.aversive.cooridnated.vRipples.dHPC = [Pre.aversive.cooridnated.vRipples.dHPC ; pInc1 pDec1];
                                Post.aversive.cooridnated.vRipples.dHPC = [Post.aversive.cooridnated.vRipples.dHPC ; pInc2 pDec2];
                                clear pInc1 pDec1 surp1 pInc2 pDec2 surp2
                            end
                        end
                        % vMembers
                        for ii = 1 : numberV
                            if members(ii+numberD,i)
                                member = clusters.vHPC(ii);
                                x = spks_vHPC(spks_vHPC(:,1) == member , 2);
                                [pInc1 pDec1 surp1] = RippleModulation_assemblies(ripples.vHPC.coordinated.all,x,TS.pre);
                                [pInc2 pDec2 surp2] = RippleModulation_assemblies(ripples.vHPC.coordinated.all,x,TS.post);
                                Pre.aversive.cooridnated.vRipples.vHPC = [Pre.aversive.cooridnated.vRipples.vHPC ; pInc1 pDec1];
                                Post.aversive.cooridnated.vRipples.vHPC = [Post.aversive.cooridnated.vRipples.vHPC ; pInc2 pDec2];
                                clear pInc1 pDec1 surp1 pInc2 pDec2 surp2
                            end
                        end
                        
                        % uncoordinated vRipples
                        % dMembers
                        for ii = 1 : numberD
                            if members(ii,i)
                                member = clusters.dHPC(ii);
                                x = spks_dHPC(spks_dHPC(:,1) == member , 2);
                                [pInc1 pDec1 surp1] = RippleModulation_assemblies(ripples.vHPC.uncoordinated.all,x,TS.pre);
                                [pInc2 pDec2 surp2] = RippleModulation_assemblies(ripples.vHPC.uncoordinated.all,x,TS.post);
                                Pre.aversive.uncooridnated.vRipples.dHPC = [Pre.aversive.uncooridnated.vRipples.dHPC ; pInc1 pDec1];
                                Post.aversive.uncooridnated.vRipples.dHPC = [Post.aversive.uncooridnated.vRipples.dHPC ; pInc2 pDec2];
                                clear pInc1 pDec1 surp1 pInc2 pDec2 surp2
                            end
                        end
                        % vMembers
                        for ii = 1 : numberV
                            if members(ii+numberD,i)
                                member = clusters.vHPC(ii);
                                x = spks_vHPC(spks_vHPC(:,1) == member , 2);
                                [pInc1 pDec1 surp1] = RippleModulation_assemblies(ripples.vHPC.uncoordinated.all,x,TS.pre);
                                [pInc2 pDec2 surp2] = RippleModulation_assemblies(ripples.vHPC.uncoordinated.all,x,TS.post);
                                Pre.aversive.uncooridnated.vRipples.vHPC = [Pre.aversive.uncooridnated.vRipples.vHPC ; pInc1 pDec1];
                                Post.aversive.uncooridnated.vRipples.vHPC = [Post.aversive.uncooridnated.vRipples.vHPC ; pInc2 pDec2];
                                clear pInc1 pDec1 surp1 pInc2 pDec2 surp2
                            end
                        end
                    end
                    
                    
                end
            end
            % Put here reward
            
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
    clear cooridnated_eventDV cooridnated_eventVD segments movement
end

end
%
%
%
%% Calculation por percentage
p = 0.01; % pvalue to define significance

figure
subplot(2,11,1)
% coordinated dRipples vs dMembers
C1 = Pre.aversive.cooridnated.dRipples.dHPC(:,1)<p;
C1 = (sum(C1)/length(C1))*100;
pie([C1 100-C1])
title('Pre')
subplot(2,11,2)
C2 = Post.aversive.cooridnated.dRipples.dHPC(:,1)<p;
C2 = (sum(C2)/length(C2))*100;
pie([C2 100-C2])
title('Post')

% coordinated dRipples vs vMembers
subplot(2,11,4)
C3 = Pre.aversive.cooridnated.dRipples.vHPC(:,1)<p;
C3 = (sum(C3)/length(C3))*100;
pie([C3 100-C3])
title('Pre')
subplot(2,11,5)
C4 = Post.aversive.cooridnated.dRipples.vHPC(:,1)<p;
C4 = (sum(C4)/length(C4))*100;
pie([C4 100-C4])
title('Post')

% coordinated vRipples vs dMembers
subplot(2,11,7)
C5 = Pre.aversive.cooridnated.vRipples.dHPC(:,1)<p;
C5 = (sum(C5)/length(C5))*100;
pie([C5 100-C5])
title('Pre')
subplot(2,11,8)
C6 = Post.aversive.cooridnated.vRipples.dHPC(:,1)<p;
C6 = (sum(C6)/length(C6))*100;
pie([C6 100-C6])
title('Post')

% coordinated vRipples vs vMembers
subplot(2,11,10)
C7 = Pre.aversive.cooridnated.vRipples.vHPC(:,1)<p;
C7 = (sum(C7)/length(C7))*100;
pie([C7 100-C7])
title('Pre')
subplot(2,11,11)
C8 = Post.aversive.cooridnated.vRipples.vHPC(:,1)<p;
C8 = (sum(C8)/length(C8))*100;
pie([C8 100-C8])
title('Post')

% uncoordinated dRipples vs dMembers
subplot(2,11,12)
C9 = Pre.aversive.uncooridnated.dRipples.dHPC(:,1)<p;
C9 = (sum(C9)/length(C9))*100;
pie([C9 100-C9])
title('Pre')
subplot(2,11,13)
C10 = Post.aversive.uncooridnated.dRipples.dHPC(:,1)<p;
C10 = (sum(C10)/length(C10))*100;
pie([C10 100-C10])
title('Post')

% uncoordinated dRipples vs vMembers
subplot(2,11,15)
C11 =  Pre.aversive.uncooridnated.dRipples.vHPC(:,1)<p;
C11 = (sum(C11)/length(C11))*100;
pie([C11 100-C11])
title('Pre')
subplot(2,11,16)
C12 = Post.aversive.uncooridnated.dRipples.vHPC(:,1)<p;
C12 = (sum(C12)/length(C12))*100;
pie([C12 100-C12])
title('Post')

% uncoordinated vRipples vs dMembers
subplot(2,11,18)
C13 = Pre.aversive.uncooridnated.vRipples.dHPC(:,1)<p;
C13 = (sum(C13)/length(C13))*100;
pie([C13 100-C13])
title('Pre')
subplot(2,11,19)
C14 = Post.aversive.uncooridnated.vRipples.dHPC(:,1)<p;
C14 = (sum(C14)/length(C14))*100;
pie([C14 100-C14])
title('Post')

% uncoordinated vRipples vs vMembers
subplot(2,11,21)
C15 = Pre.aversive.uncooridnated.vRipples.vHPC(:,1)<p;
C15 = (sum(C15)/length(C15))*100;
pie([C15 100-C15])
title('Pre')
subplot(2,11,22)
C16 = Post.aversive.uncooridnated.vRipples.vHPC(:,1)<p;
C16 = (sum(C16)/length(C16))*100;
pie([C16 100-C16])
title('Post')
% 