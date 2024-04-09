clear
clc
% close all

%% Parameters
path = {'E:\Rat126\Ephys\in_Pyr';'E:\Rat103\usable';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path

%Sleep
time_criteria = 600; %time criteria to define the maximal time of sleep to include

% Ripples
shuf = false; % to shuffle ripples to see if the coordination level change
downSamp = false; % to downsample dRipples
q = 0.25; %quantile to restrict above it ripples according to their peak amplitude
ripples_coordinated_numbers = []; %for storing the number of cooridnated events all conditions pooled
ripples_coordinated_numbers_shuffled = [];
ripples_coordinated_numbers_downsampled = [];
rate.V.all = []; rate.D.all = []; % rate
rate.V.B = []; rate.D.B = []; % rate baseline
rate.V.R = []; rate.D.R = []; % rate reward
rate.V.A = []; rate.D.A = []; % rate aversive

IRI.D.all = []; IRI.D.B = []; IRI.D.R = []; IRI.D.A = [];
IRI.V.all = []; IRI.V.B = []; IRI.V.R = []; IRI.V.A = [];

% for ripples in general
BurstIndex.D.B = []; BurstIndex.D.R = []; BurstIndex.D.A = [];
BurstIndex.V.B = []; BurstIndex.V.R = []; BurstIndex.V.A = [];
BurstIndex.D.all = []; BurstIndex.V.all = [];

%for cooridnated and uncooridnated events
BurstIndex.D.coordinated.B = []; BurstIndex.D.coordinated.R = []; BurstIndex.D.coordinated.A = [];
BurstIndex.V.coordinated.B = []; BurstIndex.V.coordinated.R = []; BurstIndex.V.coordinated.A = [];
BurstIndex.D.uncooridnated.B = []; BurstIndex.D.uncooridnated.R = []; BurstIndex.D.uncooridnated.A = [];
BurstIndex.V.uncoordinated.B = []; BurstIndex.V.uncoordinated.R = []; BurstIndex.V.uncoordinated.A = [];
BurstIndex.D.coordinated.all = []; BurstIndex.V.coordinated.all = [];
BurstIndex.D.uncoordinated.all = []; BurstIndex.V.uncoordinated.all = [];


durations.D.all = []; durations.D.B = []; durations.D.R = []; durations.D.A = [];
durations.V.all = []; durations.V.B = []; durations.V.R = []; durations.V.A = [];

amplitude.D.all = []; amplitude.D.B = []; amplitude.D.R = []; amplitude.D.A = [];
amplitude.V.all = []; amplitude.V.B = []; amplitude.V.R = []; amplitude.V.A = [];
amplitude.D.coordinated = []; amplitude.V.coordinated = []; 
amplitude.D.uncoordinated = []; amplitude.V.uncoordinated = []; 

coordinatedD_IRI_B = []; coordinatedD_IRI_R = []; coordinatedD_IRI_A = [];
coordinatedV_IRI_B = []; coordinatedV_IRI_R = []; coordinatedV_IRI_A = [];

dRipples_coordinated_single = [];   vRipples_coordinated_single = [];
dRipples_uncoordinated_single = []; vRipples_uncoordinated_single = [];
dRipples_coordinated_single_ds = [];   vRipples_coordinated_single_ds = [];
dRipples_uncoordinated_single_ds = []; vRipples_uncoordinated_single_ds = [];

Auto.dRipples.all = [] ; Auto.dRipples.aversive = [] ; Auto.dRipples.reward = [] ; Auto.dRipples.baseline = [] ; 
Auto.vRipples.all = [] ; Auto.vRipples.aversive = [] ; Auto.vRipples.reward = [] ; Auto.vRipples.baseline = [] ; 
Auto.coordinated.dRipples = []; Auto.coordinated.vRipples = []; Auto.uncoordinated.dRipples = []; Auto.uncoordinated.vRipples = [];
Auto.coordinated.subsample.dRipples = []; Auto.coordinated.subsample.vRipples = []; Auto.uncoordinated.subsample.dRipples = []; Auto.uncoordinated.subsample.vRipples = [];
Cross.DV.all = []; Cross.DV.B = []; Cross.DV.R = []; Cross.DV.A = []; 

Cross.CoorDorsal_allVentral = []; Cross.CoorVentral_allDorsal = []; Cross.CoorDorsal_allDorsal = []; Cross.CoorVentral_allVentral = [];
Cross.bursts = [];

PrefPos.dHPC = []; PrefPos.vHPC = [];
T = [];

bursts_durations.dHPC = [];   bursts_durations.vHPC = [];
Cross.t2vst1 = [];
BurstOverlapping = [];

%% Main tloop to iterate across sessions
for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    for t = 1 : length(subFolders)-2
        disp(' ')
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)
        
        disp('Uploading session time stamps')
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
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);    WAKE.all = ToIntervals(states==1);
        clear x states
        
        %keep only WAKE in HomeCage
        %         WAKE = Restrict(WAKE, [aversiveTS ; rewardTS ; baselineTS] ./1000);
        
        % NREM events restriction according conditions
        NREM.B = Restrict(NREM.all,baselineTS./1000);
        NREM.R = Restrict(NREM.all,rewardTS./1000);
        NREM.A = Restrict(NREM.all,aversiveTS./1000);
        % REM events restriction according conditions
        REM.B = Restrict(REM.all,baselineTS./1000);
        REM.R = Restrict(REM.all,rewardTS./1000);
        REM.A = Restrict(REM.all,aversiveTS./1000);
        
        %% Upload ripples
        if isfile('ripplesD_customized2.csv')
            disp('Uploading dRipples')
            ripplesD = table2array(readtable('ripplesD_customized2.csv'));
            D = true;
        else
            disp('No data related to dRipples')
            D = false;
        end
        
        if isfile('ripplesV_customized2.csv')
            disp('Uploading vRipples')
            ripplesV = table2array(readtable('ripplesV_customized2.csv'));
            V = true;
        else
            disp('No data related to vRipples')
            V = false;
        end
        
        %% ---- dRipples ----
        if D 
            disp('Lets analyze dRipples')
            rate.D.all = [rate.D.all ; length(ripplesD)/sum(NREM.all(:,2)-NREM.all(:,1))];
            durations.D.all = [durations.D.all ; ripplesD(:,3)-ripplesD(:,1)];
            amplitude.D.all = [amplitude.D.all ; ripplesD(:,4)];
            
            % IRI
            % all
            tmp = diff(ripplesD(:,2)); %pool of inter-dRipples-interval
            IRI.D.all = [IRI.D.all ; tmp];
            clear tmp
            %Aversive
            tmp = Restrict(ripplesD(:,2),NREM.A);
            tmp = diff(tmp); %pool of inter-dRipples-interval
            IRI.D.A = [IRI.D.A ; tmp];
            clear tmp
            %Reward
            tmp = Restrict(ripplesD(:,2),NREM.R);
            tmp = diff(tmp); %pool of inter-dRipples-interval
            IRI.D.R = [IRI.D.R ; tmp];
            clear tmp
            %Baseline
            tmp = Restrict(ripplesD(:,2),NREM.B);
            tmp = diff(tmp); %pool of inter-dRipples-interval
            IRI.D.B = [IRI.D.B ; tmp];
            clear tmp
            
            %Durations per condition
            tmp = Restrict(ripplesD,NREM.B);
            durations.D.B = [durations.D.B ; tmp(:,3)-tmp(:,1)];
            tmp1 = Restrict(ripplesD,NREM.R);
            durations.D.R = [durations.D.R ; tmp1(:,3)-tmp1(:,1)];
            tmp2 = Restrict(ripplesD,NREM.A);
            durations.D.A = [durations.D.A ; tmp2(:,3)-tmp2(:,1)];

            %Amplitude per conditions
            amplitude.D.B = [amplitude.D.B ; tmp(:,4)];
            amplitude.D.R = [amplitude.D.R ; tmp1(:,4)];
            amplitude.D.A = [amplitude.D.A ; tmp2(:,4)];
            rate.D.B = [rate.D.B ; length(tmp)/sum(NREM.B(:,2)-NREM.B(:,1))];
            rate.D.R = [rate.D.R ; length(tmp1)/sum(NREM.R(:,2)-NREM.R(:,1))];
            rate.D.A = [rate.D.A ; length(tmp2)/sum(NREM.A(:,2)-NREM.A(:,1))];
            clear tmp tmp1 tmp2
            
            %dIRI per condtion
            % Baselmine
            tmp = diff(Restrict(ripplesD(:,2),NREM.B)); %pool of inter-dRipples-interval
            tmp1 = 0;
            for p = 1:length(tmp)
                if and(tmp(p) <= 0.20 , tmp(p)>=0.00)
                    tmp1 = tmp1 + 1;
                end
            end
            BurstIndex.D.B = [BurstIndex.D.B ; tmp1/length(tmp)];
            clear tmp N EDGES ii p iii tmp1
            
            % Reward
            tmp = diff(Restrict(ripplesD(:,2),NREM.R)); %pool of inter-dRipples-interval
            tmp1 = 0;
            for p = 1:length(tmp)
                if and(tmp(p) <= 0.20 , tmp(p)>=0.00)
                    tmp1 = tmp1 + 1;
                end
            end
            BurstIndex.D.R = [BurstIndex.D.R ; tmp1/length(tmp)];
            clear tmp N EDGES ii p iii tmp1
            
            % Aversive
            tmp = diff(Restrict(ripplesD(:,2),NREM.A)); %pool of inter-dRipples-interval
            tmp1 = 0;
            for p = 1:length(tmp)
                if and(tmp(p) <= 0.20 , tmp(p)>=0.00)
                    tmp1 = tmp1 + 1;
                end
            end
            BurstIndex.D.A = [BurstIndex.D.A ; tmp1/length(tmp)];
            clear tmp N EDGES ii p iii tmp1
            
            
            % Auto-correlogram dRipples
            % All
            x = Restrict(ripplesD(:,2),NREM.all);
            y = Restrict(ripplesD(:,2),NREM.all);
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',1,'mode','ccg');
            
            Auto.dRipples.all = [Auto.dRipples.all , ccg(:,1,1)];
            
            % Aversive
            x = Restrict(ripplesD(:,2),NREM.A);
            y = Restrict(ripplesD(:,2),NREM.A);
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',1,'mode','ccg');
            
            Auto.dRipples.aversive = [Auto.dRipples.aversive , ccg(:,1,1)];
            
            % Reward
            x = Restrict(ripplesD(:,2),NREM.R);
            y = Restrict(ripplesD(:,2),NREM.R);
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',1,'mode','ccg');
            
            Auto.dRipples.reward = [Auto.dRipples.reward , ccg(:,1,1)];
            
            % Baseline
            x = Restrict(ripplesD(:,2),NREM.B);
            y = Restrict(ripplesD(:,2),NREM.B);
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',1,'mode','ccg');
            
            Auto.dRipples.baseline = [Auto.dRipples.baseline , ccg(:,1,1)];    
            
            % Burst Index
            count = 0;
            y = ripplesD(:,2); %exclude coordinated + associated ripples
            for s = 1 : length(ripplesD(:,2))
                seed = ripplesD(s,2);
                tmp = sum(and(y > seed-0.1 , y < seed+0.1));
                if tmp>1
                    count = count + 1;
                end
            end
            BurstIndex.D.all = [BurstIndex.D.all ; count/length(x)];
            clear ccg time x y count
            
        end
        
        %% ---- vRipples ----
        if V
            disp('Lets analyze vRipples')
            rate.V.all = [rate.V.all ; length(ripplesV)/sum(NREM.all(:,2)-NREM.all(:,1))];
            durations.V.all = [durations.V.all ; ripplesV(:,3)-ripplesV(:,1)];
            amplitude.V.all = [amplitude.V.all ; ripplesV(:,4)];
            
            % IRI
            % all
            tmp = diff(ripplesV(:,2)); %pool of inter-dRipples-interval
            IRI.V.all = [IRI.V.all ; tmp];
            clear tmp
            %Aversive
            tmp = Restrict(ripplesV(:,2),NREM.A);
            tmp = diff(tmp); %pool of inter-dRipples-interval
            IRI.V.A = [IRI.V.A ; tmp];
            clear tmp
            %Reward
            tmp = Restrict(ripplesV(:,2),NREM.R);
            tmp = diff(tmp); %pool of inter-dRipples-interval
            IRI.V.R = [IRI.V.R ; tmp];
            clear tmp
            %Baseline
            tmp = Restrict(ripplesV(:,2),NREM.B);
            tmp = diff(tmp); %pool of inter-dRipples-interval
            IRI.V.B = [IRI.V.B ; tmp];
            clear tmp
            
            %Durations per condition
            tmp = Restrict(ripplesV,NREM.B);
            durations.V.B = [durations.V.B ; tmp(:,3)-tmp(:,1)];
            tmp1 = Restrict(ripplesV,NREM.R);
            durations.V.R = [durations.V.R ; tmp1(:,3)-tmp1(:,1)];
            tmp2 = Restrict(ripplesV,NREM.A);
            durations.V.A = [durations.V.A ; tmp2(:,3)-tmp2(:,1)];

            %Amplitude per conditions
            amplitude.V.B = [amplitude.V.B ; tmp(:,4)];
            amplitude.V.R = [amplitude.V.R ; tmp1(:,4)];
            amplitude.V.A = [amplitude.V.A ; tmp2(:,4)];
            rate.V.B = [rate.V.B ; length(tmp)/sum(NREM.B(:,2)-NREM.B(:,1))];
            rate.V.R = [rate.V.R ; length(tmp1)/sum(NREM.R(:,2)-NREM.R(:,1))];
            rate.V.A = [rate.V.A ; length(tmp2)/sum(NREM.A(:,2)-NREM.A(:,1))];
            clear tmp tmp1 tmp2
            
            %dIRI per condtion
            % Baselmine
            tmp = diff(Restrict(ripplesV(:,2),NREM.B)); %pool of inter-dRipples-interval
            tmp1 = 0;
            for p = 1:length(tmp)
                if and(tmp(p) <= 0.20 , tmp(p)>=0.00)
                    tmp1 = tmp1 + 1;
                end
            end
            BurstIndex.V.B = [BurstIndex.V.B ; tmp1/length(tmp)];
            clear tmp N EDGES ii p iii tmp1
            
            % Reward
            tmp = diff(Restrict(ripplesV(:,2),NREM.R)); %pool of inter-dRipples-interval
            tmp1 = 0;
            for p = 1:length(tmp)
                if and(tmp(p) <= 0.20 , tmp(p)>=0.00)
                    tmp1 = tmp1 + 1;
                end
            end
            BurstIndex.V.R = [BurstIndex.V.R ; tmp1/length(tmp)];
            clear tmp N EDGES ii p iii tmp1
            
            % Aversive
            tmp = diff(Restrict(ripplesV(:,2),NREM.A)); %pool of inter-dRipples-interval
            tmp1 = 0;
            for p = 1:length(tmp)
                if and(tmp(p) <= 0.20 , tmp(p)>=0.00)
                    tmp1 = tmp1 + 1;
                end
            end
            BurstIndex.V.A = [BurstIndex.V.A ; tmp1/length(tmp)];
            clear tmp N EDGES ii p iii tmp1
            
            
            % Auto-correlogram vRipples
            % All
            x = Restrict(ripplesV(:,2),NREM.all);
            y = Restrict(ripplesV(:,2),NREM.all);
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',1,'mode','ccg');
            
            Auto.vRipples.all = [Auto.vRipples.all , ccg(:,1,1)];
            
            % Aversive
            x = Restrict(ripplesV(:,2),NREM.A);
            y = Restrict(ripplesV(:,2),NREM.A);
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',1,'mode','ccg');
            
            Auto.vRipples.aversive = [Auto.vRipples.aversive , ccg(:,1,1)];
            
            % Reward
            x = Restrict(ripplesV(:,2),NREM.R);
            y = Restrict(ripplesV(:,2),NREM.R);
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',1,'mode','ccg');
            
            Auto.vRipples.reward = [Auto.vRipples.reward , ccg(:,1,1)];
            
            % Baseline
            x = Restrict(ripplesV(:,2),NREM.B);
            y = Restrict(ripplesV(:,2),NREM.B);
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',1,'mode','ccg');
            
            Auto.vRipples.baseline = [Auto.vRipples.baseline , ccg(:,1,1)];
            
            % Burdt Index
            count = 0;
            y = ripplesV(:,2); %exclude coordinated + associated ripples
            for s = 1 : length(ripplesV(:,2))
                seed = ripplesV(s,2);
                tmp = sum(and(y > seed-0.1 , y < seed+0.1));
                if tmp>1
                    count = count + 1;
                end
            end
            BurstIndex.V.all = [BurstIndex.V.all ; count/length(x)];
            clear ccg time x y count
        end
        
        %% Coordinated dHPC ripples
        if and(D,V)
            coordinated = [];
            coordinatedV = [];
            coordinatedV_refined = [];
            tmpB_D = [];        tmpR_D = [];        tmpA_D = [];
            tmpB_V = [];        tmpR_V = [];        tmpA_V = [];
            for i = 1:length(ripplesD)
                r = ripplesD(i,:);
                tmp = sum(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1));
                if tmp>0
                    z = ripplesV(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1),:);
                    coordinatedV = [coordinatedV ; z];
                    [p,indice] = min(abs(r(2)-z(:,2)));
                    coordinatedV_refined = [coordinatedV_refined ; z(indice,:)];
                    coordinated = [coordinated ; r];
                    
                    if and(r(1,2)>baselineTS(1)/1000 , r(1,2)<baselineTS(2)/1000)
                        iri_t = ripplesD(and(ripplesD(:,2)>r(1,2)-0.4 , ripplesD(:,2)<r(1,2)+0.4),2);
                        tmpB_D = [tmpB_D ; iri_t]; clear iri_t
                        
                        iri_t = ripplesV(and(ripplesV(:,2)>z(indice,2)-0.4 , ripplesV(:,2)<z(indice,2)+0.4),2);
                        tmpB_V = [tmpB_V ; z(indice,2)]; clear iri_t
                    elseif and(r(1,2)>rewardTS(1)/1000 , r(1,2)<rewardTS(2)/1000)
                        iri_t = ripplesD(and(ripplesD(:,2)>r(1,2)-0.4 , ripplesD(:,2)<r(1,2)+0.4),2);
                        tmpR_D = [tmpR_D ; iri_t]; clear iri_t
                        
                        iri_t = ripplesV(and(ripplesV(:,2)>z(indice,2)-0.4 , ripplesV(:,2)<z(indice,2)+0.4),2);
                        tmpR_V = [tmpR_V ; z(indice,2)]; clear iri_t
                    elseif and(r(1,2)>aversiveTS(1)/1000 , r(1,2)<aversiveTS(2)/1000)
                        iri_t = ripplesD(and(ripplesD(:,2)>r(1,2)-0.4 , ripplesD(:,2)<r(1,2)+0.4),2);
                        tmpA_D = [tmpA_D ; iri_t]; clear iri_t
                        
                        iri_t = ripplesV(and(ripplesV(:,2)>z(indice,2)-0.4 , ripplesV(:,2)<z(indice,2)+0.4),2);
                        tmpA_V = [tmpA_V ; z(indice,2)]; clear iri_t
                    end
                    
                    clear tmp2 tmp1 p indice z
                end
                clear r
            end
            clear x tmp i

            [C,IA,IC] = unique(coordinatedV(:,1));
            coordinatedV  = coordinatedV(IA,:); clear C IA IC
            coordinatedB = Restrict(coordinated,NREM.B);    coordinatedA = Restrict(coordinated,NREM.A);    coordinatedR = Restrict(coordinated,NREM.R);
            coordinatedB_V = Restrict(coordinatedV,NREM.B);    coordinatedR_V = Restrict(coordinatedV,NREM.R);    coordinatedA_V = Restrict(coordinatedV_refined,NREM.A);
            coordinatedB_V_non_refined = Restrict(coordinatedV,NREM.B);    coordinatedR_V_non_refined  = Restrict(coordinatedV,NREM.R);    coordinatedA_V_non_refined  = Restrict(coordinatedV,NREM.A);
            %         coordinatedB_V = Restrict(coordinatedV,NREM_B);    coordinatedR_V = Restrict(coordinatedV,NREM_R);    coordinatedA_V = Restrict(coordinatedV,NREM_A);
            
            % Detection of uncoordinated ripples
            uncoordinated = ripplesD(~ismember(ripplesD(:,1),coordinated(:,1)),:);
            uncoordinatedV = ripplesV(~ismember(ripplesV(:,1),unique(coordinatedV(:,1))),:);
            
            % Save amplitude of events
            mean_D = nanmean(ripplesD(:,4));            SD_D = nanstd(ripplesD(:,4));
            mean_V = nanmean(ripplesV(:,4));            SD_V = nanstd(ripplesV(:,4));
            
            amplitude.D.coordinated = [amplitude.D.coordinated ; ((coordinated(:,4)-mean_D)./SD_D)];
            amplitude.V.coordinated = [amplitude.V.coordinated ; ((coordinatedV(:,4)-mean_V)./SD_V)];
            amplitude.D.uncoordinated = [amplitude.D.uncoordinated ; ((uncoordinated(:,4)-mean_D)./SD_D)];
            amplitude.V.uncoordinated = [amplitude.V.uncoordinated ; ((uncoordinatedV(:,4)-mean_V)./SD_V)];             
            clear mean_D mean_V SD_D SD_V
                    
            %% xcorr to study burst
            % Coordinated dRipples vs all vRipples
            x = coordinated(:,2);
            y = ripplesV(:,2);
            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
            [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
            ccg = ccg(:,1,2);
            Cross.CoorDorsal_allVentral = [Cross.CoorDorsal_allVentral , ccg]; clear ccg time s ids groups x y
            
            % Coordinated vRipples vs all dRipples
            x = unique(coordinatedV(:,2));
            y = ripplesD(:,2);
            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
            [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
            ccg = ccg(:,1,2);
            Cross.CoorVentral_allDorsal = [Cross.CoorVentral_allDorsal , ccg]; clear ccg time s ids groups x y

%             % Coordinated dRipples vs all dRipples
%             x = coordinated(:,2);
%             y = uncoordinated(:,2);
%             [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%             [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
%             ccg = ccg(:,1,2);
%             Cross.CoorDorsal_allDorsal = [Cross.CoorDorsal_allDorsal , ccg]; clear ccg time s ids groups x y
%             
%             % Coordinated vRipples vs all dRipples
%             x = unique(coordinatedV(:,2));
%             y = unique(uncoordinatedV(:,2));
%             [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%             [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
%             ccg = ccg(:,1,2);
%             Cross.CoorVentral_allVentral = [Cross.CoorVentral_allVentral , ccg]; clear ccg time s ids groups x y            

            %% Including all ripples bu eliminating the bin 0
            % Coordinated dRipples vs all dRipples
            x = coordinated(:,2);
            y = ripplesD(:,2);
            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
            [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',0,'mode','ccg');
            [i ii] = min(abs(time-0));
            ccg = ccg(:,1,2);
            ccg(ii) = 0; clear i ii
            Cross.CoorDorsal_allDorsal = [Cross.CoorDorsal_allDorsal , ccg]; clear ccg time s ids groups x y
            
            % Coordinated vRipples vs all dRipples
            x = unique(coordinatedV(:,2));
            y =ripplesV(:,2);
            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
            [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',0,'mode','ccg');
            [~, ii] = min(abs(time-0));
            ccg = ccg(:,1,2);
            ccg(ii) = 0; clear i ii
            Cross.CoorVentral_allVentral = [Cross.CoorVentral_allVentral , ccg]; clear ccg time s ids groups x y            
            
            %% For storing number of coordinated events
            ripples_coordinated_numbers = [ripples_coordinated_numbers ; length(coordinated) , length(ripplesD) , length(coordinatedV) , length(ripplesV) , length(Restrict(coordinated,NREM.B)) , length(Restrict(ripplesD,NREM.B)) , length(Restrict(coordinated,NREM.R)) , length(Restrict(ripplesD,NREM.R)) , length(Restrict(coordinated,NREM.A)) , length(Restrict(ripplesD,NREM.A)) , length(Restrict(unique(coordinatedV(:,1)),NREM.B)) , length(Restrict(ripplesV,NREM.B)) , length(Restrict(unique(coordinatedV(:,1)),NREM.R)) , length(Restrict(ripplesV,NREM.R)) , length(Restrict(unique(coordinatedV(:,1)),NREM.A)) , length(Restrict(ripplesV,NREM.A))];
            
            %% shuffle dRipples to see if that change %coordinated
            if shuf
                disp('Shuffling dRipples and calculating coordination')
                meanCoordination = [];
                for ii = 1 : 100
                    CD = [];
                    CV = [];
                    Rt = ShuffleSpks(ripplesD(:,2));
                    for i = 1:length(Rt)
                        r = Rt(i);
                        tmp = sum(and(ripplesV(:,2)>= r-0.1, ripplesV(:,2)<= r+0.1));
                        if tmp>0
                            z = ripplesV(and(ripplesV(:,2)>= r-0.1, ripplesV(:,2)<= r+0.1),:);
                            [p,indice] = min(abs(r-z(:,2)));
                            CV = [CV ; z(indice,2)];
                            CD = [CD ; r];
                            clear tmp2 tmp1 p indice z
                        end
                        clear r
                    end
                    clear x tmp i
                    
                    % Detection of uncoordinated ripples
                    UD = ripplesD(~ismember(Rt,CD));
                    UV = ripplesV(~ismember(ripplesV(:,2),CV));
                    
                    % For storing number of coordinated events
                    meanCoordination = [meanCoordination ; length(CD) , length(Rt) , length(CV) , length(ripplesV(:,2))];
                    clear Rt CV CD
                end
                ripples_coordinated_numbers_shuffled = [ripples_coordinated_numbers_shuffled ; round(mean(meanCoordination))]; clear meanCoordination
            end
            
            %% shuffle dRipples to see if that change %coordinated
            if downSamp
                disp('Downsampling Ripples and calculating coordination')
                meanCoordination = [];
                for ii = 1 : 100
                    CD = [];
                    CV = [];
                    if length(ripplesD(:,2)) > length(ripplesV(:,2))
                        Rt = ripplesD(randperm(length(ripplesD(:,2))),2);
                        Rt = sort(Rt(1:length(ripplesV(:,2))));
                        Rtt = ripplesV(:,2);
                    
                    else
                        Rtt = ripplesV(randperm(length(ripplesV(:,2))),2);
                        Rtt = sort(Rtt(1:length(ripplesD(:,2))));
                        Rt = ripplesD(:,2);
                    end
                    
                    
                    for i = 1:length(Rt)
                        r = Rt(i);
                        tmp = sum(and(Rtt>= r-0.1, Rtt<= r+0.1));
                        if tmp>0
                            z = Rtt(and(Rtt>= r-0.1, Rtt<= r+0.1),:);
                            [p,indice] = min(abs(r-z));
                            CV = [CV ; z(indice)];
                            CD = [CD ; r];
                            clear tmp2 tmp1 p indice z
                        end
                        clear r
                    end
                    clear x tmp i
                    
                    % Detection of uncoordinated ripples
                    UD = Rt(~ismember(Rt,CD));
                    UV = Rtt(~ismember(Rtt,CV));
                    
                    % For storing number of coordinated events
                    meanCoordination = [meanCoordination ; length(CD) , length(Rt) , length(CV) , length(Rtt)];
                    clear Rt CV CD Rtt UD UV
                end
                ripples_coordinated_numbers_downsampled = [ripples_coordinated_numbers_downsampled ; round(mean(meanCoordination))]; clear meanCoordination
            end
            
            %% Autoccorelogram coordionated and uncooridnated without subsampling
            % Burst Index calculation as we discuss with GG
            % coordinated dRipples
            x = coordinated(:,2);
            y = coordinated(:,2);
            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
            [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
            ccg = ccg(:,1,1);
            Auto.coordinated.dRipples = [Auto.coordinated.dRipples , ccg];
            
            count = 0;
            y = ripplesD;
            bursts.dHPC.members = [];
            bursts.dHPC.events = [];
            temporal = 0;
            for s = 1 : length(x)
                seed = x(s);
                tmp = and(y(:,2) > seed-0.1 , y(:,2) < seed+0.1);
                if sum(and(y(:,2) > seed-0.1 , y(:,2) < seed+0.1))>1
                    tmp = y(tmp,:);
                    [i ii] = min(abs(tmp(:,2) - seed));
                    if not(tmp(ii,2) == temporal)
%                         PrefPos.dHPC = [PrefPos.dHPC ; ii];
                    end
                    bursts.dHPC.members = [bursts.dHPC.members ; tmp];
                    bursts.dHPC.events = [bursts.dHPC.events ; min(tmp(:,1)) ((max(tmp(:,3))-min(tmp(:,1)))/2)+min(tmp(:,1)) max(tmp(:,3))];
                    count = count + 1;
                    temporal = tmp(ii,2); clear i ii
                end
                clear seed tmp
            end
            
            [i ii iii] = unique(bursts.dHPC.members(:,2));
            bursts.dHPC.members = bursts.dHPC.members(ii,:); clear i ii iii temporal
            
            [i ii iii] = unique(bursts.dHPC.events(:,2));
            bursts.dHPC.events = bursts.dHPC.events(ii,:); clear i ii iii
            
            BurstIndex.D.coordinated.all = [BurstIndex.D.coordinated.all ; count/length(x)];
            clear ccg time x y count
            
            %% coordinated vRipples
            x = unique(coordinatedV(:,2));
            y = unique(coordinatedV(:,2));
            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
            [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
            ccg = ccg(:,1,1);
            Auto.coordinated.vRipples = [Auto.coordinated.vRipples , ccg];
            
            count = 0;
            y = ripplesV;
            bursts.vHPC.members = [];
            bursts.vHPC.events = [];
            temporal = 0;
            for s = 1 : length(x)
                seed = x(s);
                tmp = and(y(:,2) > seed-0.1 , y(:,2) < seed+0.1);
                if sum(and(y(:,2) > seed-0.1 , y(:,2) < seed+0.1))>1
                    tmp = y(tmp,:);
                    [i ii] = min(abs(tmp(:,2) - seed));
                    if not(tmp(ii,2) == temporal)
%                         PrefPos.vHPC = [PrefPos.vHPC ; ii];
                    end
                    bursts.vHPC.members = [bursts.vHPC.members ; tmp];
                    bursts.vHPC.events = [bursts.vHPC.events ; min(tmp(:,1)) ((max(tmp(:,3))-min(tmp(:,1)))/2)+min(tmp(:,1)) max(tmp(:,3))];
                    count = count + 1;
                    temporal = tmp(ii,2); clear i ii
                end
                clear seed tmp
            end
            
            [i ii iii] = unique(bursts.vHPC.members(:,2));
            bursts.vHPC.members = bursts.vHPC.members(ii,:); clear i ii iii
            
            [i ii iii] = unique(bursts.vHPC.events(:,2));
            bursts.vHPC.events = bursts.vHPC.events(ii,:); clear i ii iii
            
            BurstIndex.V.coordinated.all = [BurstIndex.V.coordinated.all ; count/length(x)];
            clear ccg time x y count
            
            %% uncoordinated dRipples
            x = uncoordinated(:,2);
            y = uncoordinated(:,2);
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
            ccg = ccg(:,1,1);
            Auto.uncoordinated.dRipples = [Auto.uncoordinated.dRipples  , ccg];
            
            count = 0;
            y = ripplesD; %exclude coordinated + associated ripples
            bursts.uncoordinated.dHPC = [];
            for s = 1 : length(x)
                seed = x(s);
                tmp = and(y(:,2) > seed-0.1 , y(:,2) < seed+0.1);
                if sum(and(y(:,2) > seed-0.1 , y(:,2) < seed+0.1))>1
                    tmp = y(tmp,:);
                    bursts.uncoordinated.dHPC = [bursts.uncoordinated.dHPC ; tmp(1,1) tmp(end,3)];
                    count = count + 1;
                end
                clear seed tmp
            end
            
            [i ii iii] = unique(bursts.uncoordinated.dHPC(:,2));
            bursts.uncoordinated.dHPC = bursts.uncoordinated.dHPC(ii,:); clear i ii iii
            
            BurstIndex.D.uncoordinated.all = [BurstIndex.D.uncoordinated.all ; count/length(x)];
            clear ccg time x y count
            
            %% uncoordinated vRipples
            x = unique(uncoordinatedV(:,2));
            y = unique(uncoordinatedV(:,2));
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
            ccg = ccg(:,1,1);
            Auto.uncoordinated.vRipples = [Auto.uncoordinated.vRipples , ccg];
            
            count = 0;
            y = ripplesV;%exclude coordinated + associated ripples
            bursts.uncoordinated.vHPC = [];
            for s = 1 : length(x)
                seed = x(s);
                tmp = and(y(:,2) > seed-0.1 , y(:,2) < seed+0.1);
                if sum(and(y(:,2) > seed-0.1 , y(:,2) < seed+0.1))>1
                    tmp = y(tmp,:);
                    bursts.uncoordinated.vHPC = [bursts.uncoordinated.vHPC ;  tmp(1,1) tmp(end,3)];
                    count = count + 1;
                end
                clear seed tmp
            end
            
            [i ii iii] = unique(bursts.uncoordinated.vHPC(:,2));
            bursts.uncoordinated.vHPC = bursts.uncoordinated.vHPC(ii,:); clear i ii iii  
            
            BurstIndex.V.uncoordinated.all = [BurstIndex.V.uncoordinated.all ; count/length(x)];
            clear ccg time x y count
     
            %%  Detection of coordinated bursts
            % compilation of dHPC members
            y = [bursts.dHPC.members(:,1)  bursts.dHPC.members(:,3)];
            % compilation of dHPC coordinated events
            x = [coordinated(:,1) coordinated(:,3)];
            % integration of both
            z = [y;x];
            % keep unique ones
            [i ii iii] = unique(z(:,2));
            z = z(ii,:); 
            % sort events
            [i ii] = sort(z(:,1));
            z = z(ii,:);
            % Merge events that are separated by 100ms
            [z] = merge_events(z, 0.1); clear x xx i ii y
            
            % compilation of vHPC members
            y = [bursts.vHPC.members(:,1)  bursts.vHPC.members(:,3)];
            % compilation of vHPC coordinated events
            x = [coordinatedV(:,1) coordinatedV(:,3)];
            % integration of both
            z1 = [y;x];
            % keep unique ones
            [i ii iii] = unique(z1(:,2));
            z1 = z1(ii,:); 
            % sort events
            [i ii] = sort(z1(:,1));
            z1 = z1(ii,:);
            % Merge events that are separated by 100ms
            [z1] = merge_events(z1, 0.1); clear x xx i ii y 

            % Compilation of dHPC and vHPC events
            x = [z ;z1]; clear z z1
            %ordering of events
            [i ii] = sort(x(:,1));
            x = x(ii,:);
            
            % Merge overlapping events
            [x xx] = ConsolidateIntervals(x);
            % Merge events that are separated by 100ms
            [bursts.coordinated] = merge_events(x, 0.1); clear x xx i ii
            
            bursts.uncoordinated.vHPC = merge_events(bursts.uncoordinated.vHPC , 0.1);
            bursts.uncoordinated.dHPC = merge_events(bursts.uncoordinated.dHPC , 0.1);
            
%             i = InIntervals([0 : 1/1250 : segments.Var1(end)/1000],coordinatd_bursts);
            
%             %% Detection of coordinated bursts
%             x = [ones(size(bursts.dHPC.members,1),1) , bursts.dHPC.members];
%             x = [x ; ones(size(coordinated,1),1) , coordinated];
%             [i ii iii] = unique(x(:,3));
%             x = x(ii,:);
%             
%             y = [ones(size(bursts.vHPC.members,1),1) , bursts.vHPC.members];
%             y = [y ; ones(size(coordinatedV,1),1) , coordinatedV];
%             [i ii iii] = unique(y(:,3));
%             y = y(ii,:);
%             
%             xy = [x ; y];
%             [i ii] = sort(xy(:,3));
%             xy = xy(ii,:); clear x y i ii
%             xy
%             bursts.coordinated = [];
%             
%             for i = 1 : length(xy)-1
%                 tmp1 = xy(i,:);
%                 tmp2 = xy(i+1,:);
%                     if abs(tmp2(3) - tmp1(3))<0.1
%                         tmp3 = [tmp1(2) tmp1(4) ; tmp2(2) tmp2(4)];
%                         bursts.coordinated = [bursts.coordinated ; min(tmp3(:,1)) , max(tmp3(:,2))];
%                         clear tmp3
%                     end
%                     clear tmp1 tmp2
%                 end
%             end   

            %% Detection of bursts across session
            iri = diff(ripplesD(:,2))<=0.1;
            iri = ToIntervals([2:length(ripplesD)] , iri);
            iri(:,1) = iri(:,1)-1;
            B.dHPC = [];
            for i = 1 : size(iri,1)
                tmp = ripplesD(iri(i,1) : iri(i,2) , :);
                B.dHPC = [B.dHPC ; tmp(1,1) tmp(1,1)+((tmp(end,2)-tmp(1,1))/2) tmp(end,2)];
                ruler = [1:size(tmp,1)]';
                PrefPos.dHPC = [PrefPos.dHPC ; ruler(ismember(tmp(:,2) , coordinated(:,2)))];
            end
            clear iri i
            bursts_durations.dHPC = [bursts_durations.dHPC ; B.dHPC(:,2)-B.dHPC(:,1)];   

            
            % ventral
            iri = diff(ripplesV(:,2))<=0.1;
            iri = ToIntervals([2:length(ripplesV)] , iri);
            iri(:,1) = iri(:,1)-1;
            B.vHPC = [];
            for i = 1 : size(iri,1)
                tmp = ripplesV(iri(i,1) : iri(i,2) , :);
                B.vHPC = [B.vHPC ; tmp(1,1) tmp(1,1)+((tmp(end,2)-tmp(1,1))/2) tmp(end,2)];
                ruler = [1:size(tmp,1)]';
                PrefPos.vHPC = [PrefPos.vHPC ; ruler(ismember(tmp(:,2) , coordinatedV(:,2)))];
            end
            clear iri i
            bursts_durations.vHPC = [bursts_durations.vHPC ; B.vHPC(:,2)-B.vHPC(:,1)];
            
            % CCG dT2 vs vt1
            x = B.dHPC(:,3);
            y = B.vHPC(:,1);
            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
            [ccg,Time] = CCG(s,ids,'binSize',0.05,'duration',1,'smooth',1,'mode','ccg');
            ccg = ccg(:,1,2);
            Cross.t2vst1 = [Cross.t2vst1 , ccg];
            
            % both
            for i = 1 : size(B.dHPC,1)
                T = [T ; B.dHPC(i,3) - B.vHPC(:,1)];
            end
%             histogram(tmp , 'BinEdges' , [-0.5 : 0.05 : 0.5])
            
            B.both = [B.dHPC(:,1) ,  B.dHPC(:,3) ; B.vHPC(:,1) , B.vHPC(:,3)];
            [i ii] = sort(B.both(:,1));
            B.both = ConsolidateIntervals(B.both(ii,:)); clear i ii
            
%             %% CCG between bursts
%             x = bursts.dHPC.events(:,2);
%             y = bursts.vHPC.events(:,2);
%             [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%             [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',1,'smooth',1,'mode','ccg');
%             ccg = ccg(:,1,2);
%             Cross.bursts = [Cross.bursts , ccg]; clear s ids x y groups
            
            x = B.dHPC(:,2);
            y = B.vHPC(:,2);
            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
            [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',1,'smooth',1,'mode','ccg');
            ccg = ccg(:,1,2);
            Cross.bursts = [Cross.bursts , ccg]; clear s ids x y groups
            
            % Overlap beween burst
            x = [0 : 0.001 : segments.Var1(end)/1000];
            x1 = InIntervals(x,[B.dHPC(:,1) B.dHPC(:,3)]);
            x2 = InIntervals(x,[B.vHPC(:,1) B.vHPC(:,3)]);
            x3 = and(x1,x2); x3 = ToIntervals(x,x3);
            x3 = x3(:,2) - x3(:,1);
            BurstOverlapping = [BurstOverlapping ; x3]; clear x x1 x2 x3
            
            %% save bursts timestamps
            if (isfile('coordinated_ripple_bursts.mat'))
                delete([cd,'\coordinated_ripple_bursts.mat'])
            end
            save ([cd,'\coordinated_ripple_bursts.mat'],'bursts' , 'B')
            
            %% Autoccorelogram coordionated and uncooridnated with subsampling
            % Burst Index calculation as we discuss with GG
            % coordinated dRipples
            x = coordinated(:,2);
            y = coordinated(:,2);
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
            ccg = ccg(:,1,1);
            Auto.coordinated.subsample.dRipples = [Auto.coordinated.subsample.dRipples , ccg];
            clear s ids groups ccg
            
            % coordinated vRipples
            x = unique(coordinatedV(:,2));
            y = unique(coordinatedV(:,2));
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
            ccg = ccg(:,1,1);
            Auto.coordinated.subsample.vRipples = [Auto.coordinated.subsample.vRipples , ccg];
            clear s ids groups ccg

            % uncoordinated dRipples
            tmp = [];
            for i = 1:50
                x = uncoordinated(randperm(length(uncoordinated)),2);
                x = x(1:length(coordinated));
                x = sort(x);
                y = x;
                [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
                ccg = ccg(:,1,1);
                tmp = [tmp , ccg];
                clear x y ccg time s ids groups
            end
            Auto.uncoordinated.subsample.dRipples = [Auto.uncoordinated.subsample.dRipples , mean(tmp,2)];
            clear ccg time x y count tmp
            
            % uncoordinated vRipples
            if length(unique(coordinatedV(:,2))) < length(uncoordinatedV)
                tmp = [];
                for i = 1:50
                    x = uncoordinatedV(randperm(length(uncoordinatedV)),2);
                    x = x(1:length(unique(coordinatedV(:,2))));
                    x = sort(x);
                    y = x;
                    [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                    [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
                    ccg = ccg(:,1,1);
                    tmp = [tmp , ccg];
                    clear x y ccg time s ids groups
                end
            else
                tmp = [];
                for i = 1:50
                    x = uncoordinatedV(:,2);
                    y = x;
                    [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                    [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
                    ccg = ccg(:,1,1);
                    tmp = [tmp , ccg];
                    clear x y ccg time s ids groups
                end
            end
            
            Auto.uncoordinated.subsample.vRipples = [Auto.uncoordinated.subsample.vRipples , mean(tmp,2)];
            clear ccg time x y count tmp

            %% Crosscorrelogram dRipples - vRipples
            x = ripplesD(:,2);
            y = ripplesV(:,2);
            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
            [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
            ccg = ccg(:,1,2);
            Cross.DV.all = [Cross.DV.all , ccg]; clear x y s ids groups
        end
        clear ripplesD ripplesV
        clear aversiveTS aversiveTS_run rewardTS rewardTS_run baselineTS
        clear coordinated coordinatedA coordinatedB coordinatedR
        clear coordinatedV coordinatedA_V coordinatedB_V coordinatedR_V
        clear coordinatedA_V_non_refined coordinatedB_V_non_refined coordinatedR_V_non_refined
        clear uncoordinated uncoordinatedA uncoordinatedA_V uncoordinatedB uncoordinatedB_V
        clear uncoordinatedR uncoordinatedR_V uncoordinatedV
        clear REM REM_A REM_B REM_R NREM NREM_A NREM_B NREM_R WAKE D V A
        clear burstD burstV ccg coordinated_ripple_bursts q R tmp tmpA_D tmpA_V tmpB_D tmpB_V tmpR_D tmpR_V
    end
end
% save([cd,'\Ripples_Characterization_all_ratsVf.mat'])

%% Plotting Ripple rate
figure
x = [[rate.D.all ; rate.V.all] , [ones(length(rate.D.all),1) ; ones(length(rate.V.all),1)*2]];
scatter(x(:,2),x(:,1),"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 3]),hold on
scatter([1 2] , [nanmean(rate.D.all) nanmean(rate.V.all)],'filled'),ylim([0 2])
[h p] = ranksum(rate.D.all , rate.V.all)

figure
subplot(121)
x = [[rate.D.B ; rate.D.R ; rate.D.A] , [ones(length(rate.D.B),1) ; ones(length(rate.D.R),1)*2 ; ones(length(rate.D.A),1)*3]];
scatter(x(:,2),x(:,1),"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 4]),hold on
scatter([1 2 3] , [nanmean(rate.D.B) nanmean(rate.D.R) nanmean(rate.D.A)],'filled'),ylim([0 2])
% [p t stats] = kruskalwallis(x(:,1) , x(:,2))
% c = multcompare(stats)

subplot(122)
x = [[rate.V.B ; rate.V.R ; rate.V.A] , [ones(length(rate.V.B),1) ; ones(length(rate.V.R),1)*2 ; ones(length(rate.V.A),1)*3]];
scatter(x(:,2),x(:,1),"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 4]),hold on
scatter([1 2 3] , [nanmean(rate.V.B) nanmean(rate.V.R) nanmean(rate.V.A)],'filled'),ylim([0 2])
% [p t stats] = kruskalwallis(x(:,1) , x(:,2))
% c = multcompare(stats)


%% Plotting Crosscorrelogram
figure
plot(time,Cross.DV.all./sum(Cross.DV.all),'k'),hold on
plot(time,mean(Cross.DV.all,2)./sum(mean(Cross.DV.all,2)),'r')
xline(-0.2), xline(0.2)

%% Plotting Autocorrelograms
figure
subplot(131)
plot(tttt,mean(Auto.dRipples.all./sum(Auto.dRipples.all),2),'k'),hold on
plot(tttt,mean(Auto.vRipples.all./sum(Auto.vRipples.all),2),'r'),xlim([-1 1]),ylim([0 0.0055])
subplot(132)
plot(tttt,mean(Auto.dRipples.baseline./sum(Auto.dRipples.baseline),2),'k'),hold on
plot(tttt,mean(Auto.dRipples.aversive./sum(Auto.dRipples.aversive),2),'r'),xlim([-1 1]),hold on
plot(tttt,mean(Auto.dRipples.reward./sum(Auto.dRipples.reward),2),'b'),xlim([-1 1]),ylim([0 0.0055]),hold on
subplot(133)
plot(tttt,mean(Auto.vRipples.baseline./sum(Auto.vRipples.baseline),2),'k'),hold on
plot(tttt,mean(Auto.vRipples.aversive./sum(Auto.vRipples.aversive),2),'r'),xlim([-1 1]),hold on
plot(tttt,mean(Auto.vRipples.reward./sum(Auto.vRipples.reward),2),'b'),xlim([-1 1]),ylim([0 0.0055]),hold on

%% Plotting Inter-Ripple-Interval
figure
histogram(IRI.D.all,[0 : 0.01 : 1],'Normalization','probability')
hold on
histogram(IRI.V.all,[0 : 0.01 : 1],'Normalization','probability')
xlim([0 1])

%% Plotting BurstIndex
figure
x = [[BurstIndex.D.all ; BurstIndex.V.all] , [ones(length(BurstIndex.D.all),1) ; ones(length(BurstIndex.V.all),1)*2]];
scatter(x(:,2),x(:,1),"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 3]),hold on
scatter([1 2] , [nanmean(BurstIndex.D.all) nanmean(BurstIndex.V.all)],'filled'),ylim([0 2])
[h p] = ranksum(BurstIndex.D.all , BurstIndex.V.all)


%% Plotting Ripple durations
figure
histogram(durations.D.all,[0.015 : 0.005 : 0.1],'Normalization','probability')
hold on
histogram(durations.V.all,[0.015 : 0.005 : 0.1],'Normalization','probability')
xlim([0.015 0.1])


%% Plot Auto-Cooridnated
figure
subplot(121)
x = Auto.coordinated.dRipples ./ sum(Auto.coordinated.dRipples);
plot(time,mean(x,2)),hold on
x = Auto.uncoordinated.dRipples ./ sum(Auto.uncoordinated.dRipples);
plot(time,mean(x,2)),ylim([0 0.018])

subplot(122)
x = Auto.coordinated.vRipples ./ sum(Auto.coordinated.vRipples);
plot(time,mean(x,2)),hold on
x = Auto.uncoordinated.vRipples ./ sum(Auto.uncoordinated.vRipples);
plot(time,mean(x,2)),ylim([0 0.018])


% Burst Index
figure
subplot(121)
x = [[BurstIndex.D.uncoordinated.all ; BurstIndex.D.coordinated.all] , [ones(length(BurstIndex.D.uncoordinated.all),1) ; ones(length(BurstIndex.D.coordinated.all),1)*2]];
scatter(x(:,2),x(:,1),"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 3]),hold on
scatter([1 2] , [nanmean(BurstIndex.D.uncoordinated.all) nanmean(BurstIndex.D.coordinated.all)],'filled'),ylim([0 0.4])
[h p] = ranksum(BurstIndex.D.uncoordinated.all , BurstIndex.D.coordinated.all,'tail','left')

subplot(122)
x = [[BurstIndex.V.uncoordinated.all ; BurstIndex.V.coordinated.all] , [ones(length(BurstIndex.V.uncoordinated.all),1) ; ones(length(BurstIndex.V.coordinated.all),1)*2]];
scatter(x(:,2),x(:,1),"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 3]),hold on
scatter([1 2] , [nanmean(BurstIndex.V.uncoordinated.all) nanmean(BurstIndex.V.coordinated.all)],'filled'),ylim([0 0.4])
[h p] = ranksum(BurstIndex.V.uncoordinated.all , BurstIndex.V.coordinated.all,'tail','left')

%% Plot Auto-Cooridnated subsampled
subplot(121)
x = Auto.coordinated.subsample.dRipples ./ sum(Auto.coordinated.subsample.dRipples);
plot(time,mean(x,2)),hold on
x = Auto.uncoordinated.subsample.dRipples ./ sum(Auto.uncoordinated.subsample.dRipples);
plot(time,mean(x,2)),ylim([0 0.018])

subplot(122)
x = Auto.coordinated.subsample.vRipples ./ sum(Auto.coordinated.subsample.vRipples);
plot(time,mean(x,2)),hold on
x = Auto.uncoordinated.subsample.vRipples ./ sum(Auto.uncoordinated.subsample.vRipples);
plot(time,mean(x,2)),ylim([0 0.018])


%% plot percentage 
figure
x = [ (ripples_coordinated_numbers(:,1) ./ ripples_coordinated_numbers(:,2))*100 ;  (ripples_coordinated_numbers(:,3) ./ ripples_coordinated_numbers(:,4))*100];
x = [x , [ones(40,1) ; ones(40,1)*2]];
scatter(x(:,2),x(:,1),"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 3]),hold on
scatter([1 2] , [nanmean((ripples_coordinated_numbers(:,1) ./ ripples_coordinated_numbers(:,2))*100) nanmean((ripples_coordinated_numbers(:,3) ./ ripples_coordinated_numbers(:,4))*100)],'filled')%,ylim([0 0.4])
boxplot(x(:,1),x(:,2)),ylim([0 60])
% [h p] = ranksum((ripples_coordinated_numbers(:,1) ./ ripples_coordinated_numbers(:,2))*100 ,  (ripples_coordinated_numbers(:,3) ./ ripples_coordinated_numbers(:,4))*100)

figure,
subplot(121)
x = [ (ripples_coordinated_numbers(:,5) ./ ripples_coordinated_numbers(:,6))*100 ;  (ripples_coordinated_numbers(:,7) ./ ripples_coordinated_numbers(:,8))*100;  (ripples_coordinated_numbers(:,9) ./ ripples_coordinated_numbers(:,10))*100];
x = [x , [ones(40,1) ; ones(40,1)*2 ; ones(40,1)*3]];
scatter(x(:,2),x(:,1),"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 4]),hold on
scatter([1 2 3] , [nanmean((ripples_coordinated_numbers(:,5) ./ ripples_coordinated_numbers(:,6))*100) ;  nanmean((ripples_coordinated_numbers(:,7) ./ ripples_coordinated_numbers(:,8))*100) ;  nanmean((ripples_coordinated_numbers(:,9) ./ ripples_coordinated_numbers(:,10))*100)],'filled'),ylim([0 70])
boxplot(x(:,1),x(:,2)),ylim([0 60])
% [h p] = kruskalwallis(x(:,1) , x(:,2))

subplot(122)
x = [ (ripples_coordinated_numbers(:,11) ./ ripples_coordinated_numbers(:,12))*100 ;  (ripples_coordinated_numbers(:,13) ./ ripples_coordinated_numbers(:,14))*100;  (ripples_coordinated_numbers(:,15) ./ ripples_coordinated_numbers(:,16))*100];
x = [x , [ones(40,1) ; ones(40,1)*2 ; ones(40,1)*3]];
scatter(x(:,2),x(:,1),"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 4]),hold on
scatter([1 2 3] , [nanmean((ripples_coordinated_numbers(:,11) ./ ripples_coordinated_numbers(:,12))*100) ;  nanmean((ripples_coordinated_numbers(:,13) ./ ripples_coordinated_numbers(:,14))*100) ;  nanmean((ripples_coordinated_numbers(:,15) ./ ripples_coordinated_numbers(:,16))*100)],'filled'),ylim([0 70])
boxplot(x(:,1),x(:,2)),ylim([0 60])
% [h p] = kruskalwallis(x(:,1) , x(:,2))

%% plot percentage after shuffling dRipples
figure
x = [ (ripples_coordinated_numbers_shuffled(:,1) ./ ripples_coordinated_numbers_shuffled(:,2))*100 ;  (ripples_coordinated_numbers_shuffled(:,3) ./ ripples_coordinated_numbers_shuffled(:,4))*100];
x = [x , [ones(40,1) ; ones(40,1)*2]];
scatter(x(:,2),x(:,1),"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 3]),hold on
scatter([1 2] , [nanmean((ripples_coordinated_numbers_shuffled(:,1) ./ ripples_coordinated_numbers_shuffled(:,2))*100) nanmean((ripples_coordinated_numbers_shuffled(:,3) ./ ripples_coordinated_numbers_shuffled(:,4))*100)],'filled')%,ylim([0 0.4])
boxplot(x(:,1),x(:,2)),ylim([0 60])
[h p] = ranksum((ripples_coordinated_numbers_shuffled(:,1) ./ ripples_coordinated_numbers_shuffled(:,2))*100 ,  (ripples_coordinated_numbers_shuffled(:,3) ./ ripples_coordinated_numbers_shuffled(:,4))*100)

%% plot percentage after downsampling
figure
x = [ (ripples_coordinated_numbers_downsampled(:,1) ./ ripples_coordinated_numbers_downsampled(:,2))*100 ;  (ripples_coordinated_numbers_downsampled(:,3) ./ ripples_coordinated_numbers_downsampled(:,4))*100];
x = [x , [ones(40,1) ; ones(40,1)*2]];
scatter(x(:,2),x(:,1),"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 3]),hold on
scatter([1 2] , [nanmean((ripples_coordinated_numbers_downsampled(:,1) ./ ripples_coordinated_numbers_downsampled(:,2))*100) nanmean((ripples_coordinated_numbers_downsampled(:,3) ./ ripples_coordinated_numbers_downsampled(:,4))*100)],'filled')%,ylim([0 0.4])
boxplot(x(:,1),x(:,2)),ylim([0 60])
[h p] = ranksum((ripples_coordinated_numbers_downsampled(:,1) ./ ripples_coordinated_numbers_downsampled(:,2))*100 ,  (ripples_coordinated_numbers_downsampled(:,3) ./ ripples_coordinated_numbers_downsampled(:,4))*100)

%% Plot percentage all, downsampled, shuffled
x1 = (ripples_coordinated_numbers(:,1) ./ ripples_coordinated_numbers(:,2))*100;
x2 = (ripples_coordinated_numbers_downsampled(:,1) ./ ripples_coordinated_numbers_downsampled(:,2))*100;
x3 = (ripples_coordinated_numbers_shuffled(:,1) ./ ripples_coordinated_numbers_shuffled(:,2))*100;
x4 = (ripples_coordinated_numbers(:,3) ./ ripples_coordinated_numbers(:,4))*100;
x5 = (ripples_coordinated_numbers_downsampled(:,3) ./ ripples_coordinated_numbers_downsampled(:,4))*100;
x6 = (ripples_coordinated_numbers_shuffled(:,3) ./ ripples_coordinated_numbers_shuffled(:,4))*100;

x = [x1 ; x2 ; x3 ; x4 ; x5 ; x6];
x = [x , [ones(40,1) ; ones(40,1)*2 ; ones(40,1)*3 ; ones(40,1)*4 ; ones(40,1)*5 ; ones(40,1)*6]];
x = [x , [ones(40,1) ; ones(40,1) ; ones(40,1) ; ones(40,1)*2 ; ones(40,1)*2 ; ones(40,1)*2]];
x = [x , [ones(40,1) ; ones(40,1)*2 ; ones(40,1)*3 ; ones(40,1) ; ones(40,1)*2 ; ones(40,1)*3]];

figure
scatter(x(:,2),x(:,1),"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 7]),hold on
scatter([1 2 3 4 5 6] , [nanmean(x1) nanmean(x2) nanmean(x3) nanmean(x4) nanmean(x5) nanmean(x6)],'filled')%,ylim([0 0.4])
boxplot(x(:,1),x(:,2)),ylim([0 60])

[~,~,stats] = anovan(x(:,1) , x(:,3:4) , 'model','interaction','varnames',{'Structure','Condition'})
[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);

%% Plot CDF of amplitude across structures
% x = [amplitude.D.uncoordinated ; amplitude.D.coordinated];
% x = [x , [ones(length(amplitude.D.uncoordinated),1) ; ones(length(amplitude.D.coordinated),1)*2]];
% boxplot(x(:,1) , x(:,2))
% 
% x = [amplitude.V.uncoordinated ; amplitude.V.coordinated];
% x = [x , [ones(length(amplitude.V.uncoordinated),1) ; ones(length(amplitude.V.coordinated),1)*2]];
% boxplot(x(:,1) , x(:,2))

figure
subplot(121),cdfplot(amplitude.D.coordinated),hold on
cdfplot(amplitude.D.uncoordinated),xlim([-2 6])
[h p] = kstest2(amplitude.D.coordinated,amplitude.D.uncoordinated)
xline(3,'--')

subplot(122),cdfplot(amplitude.V.coordinated),hold on
cdfplot(amplitude.V.uncoordinated),xlim([-2 6])
xline(3,'--')
[h p] = kstest2(amplitude.V.coordinated,amplitude.V.uncoordinated)

figure, % zoom
subplot(121),cdfplot(amplitude.D.coordinated),hold on
cdfplot(amplitude.D.uncoordinated),xlim([0 3]),ylim([0.6 1])
[h p] = kstest2(amplitude.D.coordinated,amplitude.D.uncoordinated)

subplot(122),cdfplot(amplitude.V.coordinated),hold on
cdfplot(amplitude.V.uncoordinated),xlim([0 3]),ylim([0.6 1])
[h p] = kstest2(amplitude.V.coordinated,amplitude.V.uncoordinated)



figure
x = and(amplitude.D.coordinated>-2 , amplitude.D.coordinated<6);
x = amplitude.D.coordinated(x);
subplot(121),cdfplot(x),hold on
x = and(amplitude.D.uncoordinated>-2 , amplitude.D.uncoordinated<6);
x = amplitude.D.uncoordinated(x);
cdfplot(x),hold on
[h p] = kstest2(amplitude.D.coordinated,amplitude.D.uncoordinated)

x = and(amplitude.V.coordinated>-2 , amplitude.V.coordinated<6);
x = amplitude.V.coordinated(x);
subplot(122),cdfplot(x),hold on
x = and(amplitude.V.uncoordinated>-2 , amplitude.V.uncoordinated<6);
x = amplitude.V.uncoordinated(x);
cdfplot(x),hold on
[h p] = kstest2(amplitude.D.coordinated,amplitude.D.uncoordinated)

%% Plot Cross CoorDorsal_allDorsal and CoorVentral_allVentral
figure
subplot(221),plot(time , mean(Cross.CoorDorsal_allDorsal./sum(Cross.CoorDorsal_allDorsal),2)),xlim([-1 1]),ylim([0 0.012])
subplot(222),plot(time , mean(Cross.CoorVentral_allVentral./sum(Cross.CoorVentral_allVentral),2)),xlim([-1 1]),ylim([0 0.012])
subplot(223),plot(time , mean(Cross.CoorDorsal_allVentral./sum(Cross.CoorDorsal_allVentral),2)),xlim([-0.2 0.2]),ylim([0 0.025])
subplot(224),plot(time , mean(Cross.CoorVentral_allDorsal./sum(Cross.CoorVentral_allDorsal),2)),xlim([-0.2 0.2]),ylim([0 0.025])

%% x-coord bursts
figure,
plot([-0.5 : 0.01 : 0.5],mean(Cross.bursts./sum(Cross.bursts),2)),xlim([-0.5 0.5])

%%
figure
subplot(121),histogram(PrefPos.dHPC,'Normalization','probability'),xlim([0 10]),ylim([0 0.5])
subplot(122),histogram(PrefPos.vHPC,'Normalization','probability'),xlim([0 10]),ylim([0 0.5])

%% plot distribution of durations
figure
subplot(121), histogram(bursts_durations.dHPC,[0:0.01:0.4],'Normalization','probability'),ylim([0 0.4])
subplot(122), histogram(bursts_durations.vHPC,[0:0.01:0.4],'Normalization','probability'),ylim([0 0.4])

%% plot CCG dT2 vd vT1
figure
plot(Time,mean(Cross.t2vst1./sum(Cross.t2vst1),2))
