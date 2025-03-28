clear
clc
% close all

%% Parameters
path = {'E:\Rat126\Ephys\in_Pyr';'E:\Rat103\usable';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path

%Sleep
time_criteria = 1; % minimal time to include a NREM epoch (in min)

% for SU
criteria_fr = 0; %criteria to include or not a SU into the analysis
criteria_n = [3 3]; % minimal number of neurons from each structure [vHPC dHPC]
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
binSize = 0.025;
n_SU_V = 0;
n_SU_D = 0;

% Behavior
minimal_speed = 7; % minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods

EV.aversive.dvHPCi = [];   EV.reward.dvHPCi = [];
EV.aversive.dvHPC = [];   EV.reward.dvHPC = [];

FiringRate.dHPC.baseline = [];
FiringRate.dHPC.reward = [];
FiringRate.dHPC.aversive = [];
FiringRate.vHPC.baseline = [];
FiringRate.vHPC.reward = [];
FiringRate.vHPC.aversive = [];

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
        %Loading TS of the sessions
        disp('Uploading session time stamps')

        
        load('behavioral-data.mat')
        load('session_organization.mat')
        
        %% load sleep states
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);    WAKE.all = ToIntervals(states==1);
        clear x states
        
%         NREM.all(NREM.all(:,2)-NREM.all(:,1)<60*time_criteria,:)=[];
        NREM.baseline = Restrict(NREM.all,baselineTS./1000);
        NREM.aversive = Restrict(NREM.all,aversiveTS./1000);
        NREM.reward = Restrict(NREM.all,rewardTS./1000);    
        
        REM.baseline = Restrict(REM.all,baselineTS./1000);
        REM.aversive = Restrict(REM.all,aversiveTS./1000);
        REM.reward = Restrict(REM.all,rewardTS./1000);
                
        
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
        
        if and(RV, RD)
            % coordination
            coordinated = [];
            coordinatedV = [];
            coordinatedV_refined = [];
            cooridnated_event = [];
            cooridnated_eventDV = [];
            cooridnated_eventVD = [];
            for i = 1:length(ripplesD)
                r = ripplesD(i,:);
                tmp = sum(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1));
                if tmp>0
                    z = ripplesV(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1),:);
                    coordinatedV = [coordinatedV ; z];
                    [p,indice] = min(abs(r(2)-z(:,2)));
                    coordinatedV_refined = [coordinatedV_refined ; z(indice,:)];
                    coordinated = [coordinated ; r];
                    
                    cooridnated_event = [cooridnated_event ; min([r(1) , z(indice,1)]) , min(min(z(indice,2),r(2)))+abs(z(indice,2)-r(2))/2 , max([r(3) , z(indice,3)])];
                    
                    if r(2)<z(indice,2)
                        cooridnated_eventDV = [cooridnated_eventDV ; min([r(1) , z(indice,1)]) , min(min(z(indice,2),r(2)))+abs(z(indice,2)-r(2))/2 , max([r(3) , z(indice,3)])];
                    else
                        cooridnated_eventVD = [cooridnated_eventVD ; min([r(1) , z(indice,1)]) , min(min(z(indice,2),r(2)))+abs(z(indice,2)-r(2))/2 , max([r(3) , z(indice,3)])];
                    end
                    
                    clear tmp2 tmp1 p indice z
                end
                clear r
            end
            clear x tmp i
            
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
            %eliminate ripples <50 ms
%             ri = (ripples.dHPC.uncoordinated.all(:,3) - ripples.dHPC.uncoordinated.all(:,1))<0.05;
%             ripples.dHPC.uncoordinated.all(ri,:) = []; clear ri
            ripples.dHPC.coordinated.baseline = Restrict(coordinated , NREM.baseline);
            ripples.dHPC.coordinated.reward = Restrict(coordinated , NREM.reward);
            ripples.dHPC.coordinated.aversive = Restrict(coordinated , NREM.aversive);
            % coordinated vRipples
            ripples.vHPC.coordinated.all = coordinatedV_refined;
            ripples.vHPC.uncoordinated.all = ripplesV(not(ismember(ripplesV(:,2) , coordinatedV_refined(:,2))),:);
            %eliminate ripples <50 ms
%             ri = (ripples.vHPC.uncoordinated.all(:,3) - ripples.vHPC.uncoordinated.all(:,1))<0.05;
%             ripples.vHPC.uncoordinated.all(ri,:) = []; clear ri
            ripples.vHPC.coordinated.baseline = Restrict(coordinatedV_refined , NREM.baseline);
            ripples.vHPC.coordinated.reward = Restrict(coordinatedV_refined , NREM.reward);
            ripples.vHPC.coordinated.aversive = Restrict(coordinatedV_refined , NREM.aversive);
            % ALL Uncoordinated
            % downsampling to have same amount of ripples
            if size(ripples.dHPC.uncoordinated.all,1)>size(ripples.dHPC.coordinated.all,1)
                ran = randperm(size(ripples.dHPC.uncoordinated.all,1));
                ran = ripples.dHPC.uncoordinated.all(ran,:);
                ran = ran(1 : size(ripples.dHPC.coordinated.all,1) , :);
            else
                ran = ripples.dHPC.uncoordinated.all;
            end
            
            if size(ripples.vHPC.uncoordinated.all,1)>size(ripples.vHPC.coordinated.all,1)
                ran1 = randperm(size(ripples.vHPC.uncoordinated.all,1));
                ran1 = ripples.vHPC.uncoordinated.all(ran1,:);
                ran1 = ran1(1 : size(ripples.vHPC.coordinated.all,1) , :);
            else
                ran1 = ripples.vHPC.uncoordinated.all;
            end
            
            u = [ran ; ran1]; clear ran ran1
            
%             u = [ripples.dHPC.uncoordinated.all ripples.dHPC.uncoordinated.all ; ripples.vHPC.uncoordinated.all ripples.vHPC.uncoordinated.all];
            
            [B I] = sort(u(:,1));
            u = ConsolidateIntervals([u(I,1) u(I,3)]);
            
            ripple_event.uncoordinated = u; clear u
            %coordinated event
            cooridnated_event((cooridnated_event(:,3)-cooridnated_event(:,1)<0.04),:) = [];
            ripple_event.baseline = Restrict(cooridnated_event,baselineTS./1000);
            ripple_event.reward = Restrict(cooridnated_event,rewardTS./1000);
            ripple_event.aversive = Restrict(cooridnated_event,aversiveTS./1000);
            ripple_event.all = cooridnated_event;
        elseif RD
            % Store events time stamps
            % dRipples
            ripples.dHPC.baseline = Restrict(ripplesD , NREM.baseline);
            ripples.dHPC.reward = Restrict(ripplesD , NREM.reward);
            ripples.dHPC.aversive = Restrict(ripplesD , NREM.aversive);
        elseif RV
            % vRipples
            ripples.vHPC.baseline = Restrict(ripplesV , NREM.baseline);
            ripples.vHPC.reward = Restrict(ripplesV , NREM.reward);
            ripples.vHPC.aversive = Restrict(ripplesV , NREM.aversive);
        end        
        
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
            cellulartype = [K(:,1) , ones(size(K,1),1)];
        end
        % Definition of the emotional transition
        if aversiveTS_run(1)<rewardTS_run(1)
            config = 1; % 1 if the order was Aversive -> Reward
        else
            config = 2;% 1 if the order was Reward -> Aversive
        end
        
        %% Constructing Spiketrains
        freq = 1/binSize;
        limits = [0 segments.Var1(end)/1000];
        spiketrains_dHPC.pyr = [];        spiketrains_dHPC.int = [];
        spiketrains_vHPC.pyr = [];        spiketrains_vHPC.int = [];
        for ii=1:size(group_dHPC,1)
            cluster = group_dHPC(ii,1);
            Cellulartype = logical(Cell_type_classification(Cell_type_classification(:,1) == cluster,7));
            if Cellulartype
                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2);
                [tmp,bins]=binspikes(spks,freq,limits);

                spiketrains_dHPC.pyr = [spiketrains_dHPC.pyr , (tmp)];
            end
            clear spks tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5
        end
        
        for ii=1:size(group_vHPC,1)
            cluster = group_vHPC(ii,1);
            Cellulartype = logical(Cell_type_classification(Cell_type_classification(:,1) == cluster,7));
            if Cellulartype
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2);
                [tmp,bins]=binspikes(spks,freq,limits);
                
                spiketrains_vHPC.pyr = [spiketrains_vHPC.pyr , (tmp)];
            end
            clear spks tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5
        end
        clear freq limits
        clear spks spks_dHPC spks_vHPC camara shock rightvalve leftvalve
        clear ejeX ejeY dX dY dX_int dY_int
        
        %% Explained variance calculation
        if and(size(spiketrains_vHPC.pyr,2) >= criteria_n(1),size(spiketrains_dHPC.pyr,2) >= criteria_n(2))
            disp('Lets go for the SUs')
            %Restricting bins inside each condition
            
%             is.baseline.sws = InIntervals(bins,NREM.baseline);
%             is.reward.sws = InIntervals(bins,NREM.reward);
%             is.aversive.sws = InIntervals(bins,NREM.aversive);
            p = Restrict([ripple_event.baseline(:,2)-0.2 ripple_event.baseline(:,2)+0.2],[(aversiveTS_run(1)/1000)-3600 (aversiveTS_run(1)/1000)]);
            is.baseline.sws = InIntervals(bins,p);clear p
            
            p = Restrict([ripple_event.reward(:,2)-0.2 ripple_event.reward(:,2)+0.2],[(aversiveTS_run(1)/1000)-3600 (aversiveTS_run(1)/1000)]);
            is.reward.sws = InIntervals(bins,p); clear p
            
            p = Restrict([ripple_event.aversive(:,2)-0.2 ripple_event.aversive(:,2)+0.2],[(aversiveTS_run(2)/1000) (aversiveTS_run(2)/1000)+3600]);
            is.aversive.sws = InIntervals(bins,p); clear p
            
            is.shock = InIntervals(bins,[Shocks_filt Shocks_filt+1]);

            %%
            if and(and(~isempty(NREM.aversive),~isempty(NREM.reward)),~isempty(NREM.baseline))
                if aversiveTS_run(1)>rewardTS_run(1)
                    %% Aversive
                    x = [spiketrains_dHPC.pyr(is.reward.sws,:)];
                    y = [spiketrains_vHPC.pyr(is.reward.sws,:)];
                    [S3 , p] = corr(x,y);
%                     S3 = S3 - diag(diag(S3));
                    clear x y p
                    
                    % Correlation Matrix Calculation
                    x = [spiketrains_dHPC.pyr(is.shock,:)];
                    y = [spiketrains_vHPC.pyr(is.shock,:)];
                    [S4 , p] = corr(x,y);
%                     S4 = S4 - diag(diag(S4));
                    clear x y p
                    
                    x = [spiketrains_dHPC.pyr(is.aversive.sws,:)];
                    y = [spiketrains_vHPC.pyr(is.aversive.sws,:)];
                    [S5 , p] = corr(x,y);
%                     S5 = S5 - diag(diag(S5));
                    clear x y p
                    
                    % EV and REV
                    Sx = corrcoef(S4,S5,'rows','complete');     Sx = Sx(1,2);
                    Sy = corrcoef(S4,S3,'rows','complete');     Sy = Sy(1,2);
                    Sz = corrcoef(S5,S3,'rows','complete');     Sz = Sz(1,2);
                    
                    ev = (Sx-Sy*Sz);
                    ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    rev = (Sy-Sx*Sz);
                    rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    EV.aversive.dvHPC = [EV.aversive.dvHPC ; rev*100 , ev*100 tt t config];
                    clear Sx Sy Sz ev rev
                    
                    clear rev1 ev1 rev2 ev2
                    
                else
                    %% Aversive
                    % Correlation Matrix Calculation
                    x = [spiketrains_dHPC.pyr(is.baseline.sws,:)];
                    y = [spiketrains_vHPC.pyr(is.baseline.sws,:)];
                    [S1 , p] = corr(x,y);
%                     S1 = S1 - diag(diag(S1));
                    clear x y p
                    
                    x = [spiketrains_dHPC.pyr(is.shock,:)];
                    y = [spiketrains_vHPC.pyr(is.shock,:)];
                    [S2 , p] = corr(x,y);
%                     S2 = S2 - diag(diag(S2));
                    clear x y p
                    
                    x = [spiketrains_dHPC.pyr(is.aversive.sws,:)];
                    y = [spiketrains_vHPC.pyr(is.aversive.sws,:)];
                    [S3 , p] = corr(x,y);
%                     S3 = S3 - diag(diag(S3));
                    clear x y p
                    
                    % EV and REV
                    Sx = corrcoef(S2,S3,'rows','complete');     Sx = Sx(1,2);
                    Sy = corrcoef(S2,S1,'rows','complete');     Sy = Sy(1,2);
                    Sz = corrcoef(S3,S1,'rows','complete');     Sz = Sz(1,2);
                    
                    ev = (Sx-Sy*Sz);
                    ev = (ev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    rev = (Sy-Sx*Sz);
                    rev = (rev/sqrt((1-(Sy^2))*(1-(Sz^2))))^2;
                    
                    EV.aversive.dvHPC = [EV.aversive.dvHPC ; rev*100 , ev*100 tt t config];
                    clear Sx Sy Sz ev rev
                end
            end
            
        end
        disp(['-- Analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' finished --'])
        disp(' ')
        clear spiketrains_dHPC_int spiketrains_dHPC_pyr spiketrains_vHPC_int spiketrains_vHPC_pyr
        clear aversiveTS aversiveTS_run rewardTS rewardTS_run behavior bins Cell_type_classification cellulartype
        clear group_dHPC group_vHPC is K Kinfo NREM REM WAKE ripplesD ripplesV segments
        clear spiketrains_dHPC spiketrains_vHPC ripple_bursts ripple_event baselineTS cond
        clear coordinated coordinatedV coordinatedV_refined movement
    end
end

% save([cd,'\Explained_Variance_REM.mat'] , 'EV')

figure
boxplot([EV.aversive.dvHPC(:,1) , EV.aversive.dvHPC(:,2)]) , ylim([0 4]), [h p] =ranksum(EV.aversive.dvHPC(:,1) , EV.aversive.dvHPC(:,2),'tail','left'),hold on
x = [[EV.aversive.dvHPC(:,1) ; EV.aversive.dvHPC(:,2)] , [ones(length(EV.aversive.dvHPC),1) ; ones(length(EV.aversive.dvHPC),1)*2]];
scatter(x(:,2),x(:,1),"filled",'jitter','on', 'jitterAmount',0.1)

% subplot(1,4,3),boxplot([EV.reward.dvHPCi(:,1) , EV.reward.dvHPCi(:,2)]) , ylim([0 3]), [h p] =ranksum(EV.reward.dvHPCi(:,1) , EV.reward.dvHPCi(:,2),'tail','left')
% subplot(1,4,4),boxplot([EV.aversive.dvHPCi(:,1) , EV.aversive.dvHPCi(:,2)]) , ylim([0 3]), [h p] =ranksum(EV.aversive.dvHPCi(:,1) , EV.aversive.dvHPCi(:,2),'tail','left')

% figure
% subplot(1,2,1),scatter([ones(length(EV.reward.dvHPC(:,1)),1) ; ones(length(EV.reward.dvHPC(:,2)),1)*2],[EV.reward.dvHPC(:,1) ; EV.reward.dvHPC(:,2)],"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 3]),hold on
% scatter([1 2] , [nanmean(EV.reward.dvHPC(:,1)) nanmean(EV.reward.dvHPC(:,2))],'filled')
% subplot(1,2,2),scatter([ones(length(EV.aversive.dvHPC(:,1)),1) ; ones(length(EV.aversive.dvHPC(:,2)),1)*2],[EV.aversive.dvHPC(:,1) ; EV.aversive.dvHPC(:,2)],"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 3]),hold on
% scatter([1 2] , [nanmean(EV.aversive.dvHPC(:,1)) nanmean(EV.aversive.dvHPC(:,2))],'filled')


x = [EV.reward.dvHPC(:,1) ; EV.reward.dvHPC(:,2) ; EV.aversive.dvHPC(:,1) ; EV.aversive.dvHPC(:,2)];
y = [ones(length([EV.reward.dvHPC(:,1) ; EV.reward.dvHPC(:,2)]),1) ; ones(length([EV.aversive.dvHPC(:,1) ; EV.aversive.dvHPC(:,2)]),1)*2];
y = [y , [ones(length(EV.reward.dvHPC(:,1)),1) ; ones(length(EV.reward.dvHPC(:,2)),1)*2 ; ones(length(EV.aversive.dvHPC(:,1)),1) ; ones(length(EV.aversive.dvHPC(:,2)),1)*2]]

[~,~,stats] = anovan(x,y,'model','interaction','varnames',{'condition','tipo'})
[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);

x = [EV.reward.dvHPC(:,1) ; EV.reward.dvHPC(:,2) ; EV.aversive.dvHPC(:,1) ; EV.aversive.dvHPC(:,2)];
y = {'REV';'REV';'REV';'REV';'REV';'REV';'REV';'REV';'REV';'REV';'REV';'REV';'REV';'REV';'EV';'EV';'EV';'EV';'EV';'EV';'EV';'EV';'EV';'EV';'EV';'EV';'EV';'EV';'REV';'REV';'REV';'REV';'REV';'REV';'REV';'REV';'EV';'EV';'EV';'EV';'EV';'EV';'EV';'EV'};
z = {'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'R';'A';'A';'A';'A';'A';'A';'A';'A';'A';'A';'A';'A';'A';'A';'A';'A'};
[~,~,stats] = anovan(x,{y,z},'interaction')
[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);


x = EV.aversive.reward.dvHPC .* 100;
y = EV.reward.aversive.dvHPC .* 100;

figure
subplot(121),boxplot(x),ylim([0 8]), ranksum(x(:,3),x(:,4),'tail','left')
subplot(122),boxplot(y),ylim([0 8]), ranksum(y(:,3),y(:,4),'tail','left')

