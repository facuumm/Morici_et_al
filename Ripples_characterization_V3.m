clear
clc
close all

%% Parameters
path = {'E:\Rat126\Ephys\in_Pyr';'E:\Rat103\usable';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path

%Sleep
time_criteria = 600; %time criteria to define the maximal time of sleep to include

% Ripples
q = 0.25; %quantile to restrict above it ripples according to their peak amplitude
ripples_coordinated_percentage = []; %for storing percentage of coordnated events across conditions
ripples_coordinated_numbers = []; %for storing the number of cooridnated events all conditions pooled
rate.V.all = []; rate.D.all = []; % rate
rate.V.B = []; rate.D.B = []; % rate baseline
rate.V.R = []; rate.D.R = []; % rate reward
rate.V.A = []; rate.D.A = []; % rate aversive

IRI.D.all = []; IRI.D.B = []; IRI.D.R = []; IRI.D.A = [];
IRI.V.all = []; IRI.V.B = []; IRI.V.R = []; IRI.V.A = [];

% for ripples in general
BurstIndex.D.B = []; BurstIndex.D.R = []; BurstIndex.D.A = [];
BurstIndex.V.B = []; BurstIndex.V.R = []; BurstIndex.V.A = [];

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

coordinatedD_IRI_B = []; coordinatedD_IRI_R = []; coordinatedD_IRI_A = [];
coordinatedV_IRI_B = []; coordinatedV_IRI_R = []; coordinatedV_IRI_A = [];

dRipples_coordinated_single = [];   vRipples_coordinated_single = [];
dRipples_uncoordinated_single = []; vRipples_uncoordinated_single = [];
dRipples_coordinated_single_ds = [];   vRipples_coordinated_single_ds = [];
dRipples_uncoordinated_single_ds = []; vRipples_uncoordinated_single_ds = [];

Auto.dRipples.all = [] ; Auto.dRipples.aversive = [] ; Auto.dRipples.reward = [] ; Auto.dRipples.baseline = [] ; 
Auto.vRipples.all = [] ; Auto.vRipples.aversive = [] ; Auto.vRipples.reward = [] ; Auto.vRipples.baseline = [] ; 
Auto.coordinated.dRipples = []; Auto.coordinated.vRipples = []; Auto.uncoordinated.dRipples = []; Auto.uncoordinated.vRipples = [];
Cross.DV.all = []; Cross.DV.B = []; Cross.DV.R = []; Cross.DV.A = []; 


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
            
            coordinatedB = Restrict(coordinated,NREM.B);    coordinatedA = Restrict(coordinated,NREM.A);    coordinatedR = Restrict(coordinated,NREM.R);
            coordinatedB_V = Restrict(coordinatedV_refined,NREM.B);    coordinatedR_V = Restrict(coordinatedV_refined,NREM.R);    coordinatedA_V = Restrict(coordinatedV_refined,NREM.A);
            coordinatedB_V_non_refined = Restrict(coordinatedV,NREM.B);    coordinatedR_V_non_refined  = Restrict(coordinatedV,NREM.R);    coordinatedA_V_non_refined  = Restrict(coordinatedV,NREM.A);
            %         coordinatedB_V = Restrict(coordinatedV,NREM_B);    coordinatedR_V = Restrict(coordinatedV,NREM_R);    coordinatedA_V = Restrict(coordinatedV,NREM_A);
            
            % Detection of uncoordinated ripples
            uncoordinated = ripplesD(~ismember(ripplesD(:,1),coordinated(:,1)),:);
            uncoordinatedV = ripplesV(~ismember(ripplesV(:,1),coordinatedV_refined(:,1)),:);
            
%             %% Detection of ripple burst
%             coordinated_ripple_bursts = [];
%             for i = 1:length(coordinated)
%                 [p,ii] = min(abs(coordinatedV_refined(:,2) - coordinated(i,2)));
%                 
%                 %Burst detection for dorsal ripples
%                 tmp = and(ripplesD(:,2) > coordinated(i)-0.2 , ripplesD(:,2) < coordinated(i)+0.2);
%                 if sum(tmp)>1
%                     burstD = [coordinated(i,:) ; ripplesD(tmp,:)];
%                 else
%                     burstD = [coordinated(i,:)];
%                 end
%                 clear tmp
%                 
%                 %Burst detection for ventral ripples
%                 tmp = and(ripplesV(:,2) > coordinatedV_refined(ii,2)-0.2 , ripplesV(:,2) < coordinatedV_refined(ii,2)+0.2);
%                 if sum(tmp)>1
%                     burstV = [coordinatedV_refined(ii,:) ; ripplesV(tmp,:)];
%                 else
%                     burstV = [coordinatedV_refined(ii,:)];
%                 end
%                 clear p ii tmp
%                 
%                 if or(size(burstD,1)>1 , size(burstV,1)>1)
%                     %                 %To keep the overlapping area
%                     %                 tmp = [min(burstD(:,1)) , max(burstD(:,3)) ; min(burstV(:,1)) , max(burstV(:,3))];
%                     %                 tmp = [max(tmp(:,1)) , max(tmp(:,1)) + ((min(tmp(:,2))-max(tmp(:,1)))/2) , min(tmp(:,2))];
%                     % To keep the entier burst
%                     tmp = [burstD;burstV];
%                     tmp = [min(tmp(:,1)) , ((max(tmp(:,3))-min(tmp(:,1)))/2)+min(tmp(:,1)), max(tmp(:,3))];
%                     
%                     coordinated_ripple_bursts = [coordinated_ripple_bursts ; tmp];
%                     clear tmp
%                 end
%             end
% %             save ([cd,'\coordinated_ripple_bursts.mat'],'coordinated_ripple_bursts')
            
            %% For storing number of coordinated events
            ripples_coordinated_numbers = [ripples_coordinated_numbers ; length(coordinated) , length(ripplesD) , length(coordinatedV) , length(ripplesV)];
            
            %% Autoccorelogram coordionated and uncooridnated without subsampling
            % Burst Index calculation as we discuss with GG
            % coordinated dRipples
            x = coordinated(:,2);
            y = coordinated(:,2);
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
            ccg = ccg(:,1,1);
            Auto.coordinated.dRipples = [Auto.coordinated.dRipples , ccg];
            
            count = 0;
            y = ripplesD(:,2);
            for s = 1 : length(x)
                seed = x(s);
                tmp = sum(and(y > seed-0.2 , y < seed+0.2));
                if tmp>1
                    count = count + 1;
                end
            end
            BurstIndex.D.coordinated.all = [BurstIndex.D.coordinated.all ; count/length(x)];
            clear ccg time x y count
            
            % coordinated vRipples
            x = unique(coordinatedV(:,2));
            y = unique(coordinatedV(:,2));
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
            ccg = ccg(:,1,1);
            Auto.coordinated.vRipples = [Auto.coordinated.vRipples , ccg];
            
            count = 0;
            y = ripplesV(:,2);
            for s = 1 : length(x)
                seed = x(s);
                tmp = sum(and(y > seed-0.2 , y < seed+0.2));
                if tmp>1
                    count = count + 1;
                end
            end
            BurstIndex.V.coordinated.all = [BurstIndex.V.coordinated.all ; count/length(x)];
            clear ccg time x y count
            
            % uncoordinated dRipples
            x = uncoordinated(:,2);
            y = uncoordinated(:,2);
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
            ccg = ccg(:,1,1);
            Auto.uncoordinated.dRipples = [Auto.uncoordinated.dRipples  , ccg];
            
            count = 0;
            y = ripplesD(:,2); %exclude coordinated + associated ripples
            for s = 1 : length(x)
                seed = x(s);
                tmp = sum(and(y > seed-0.2 , y < seed+0.2));
                if tmp>1
                    count = count + 1;
                end
            end
            BurstIndex.D.uncoordinated.all = [BurstIndex.D.uncoordinated.all ; count/length(x)];
            clear ccg time x y count
            
            % uncoordinated vRipples
            x = unique(uncoordinatedV(:,2));
            y = unique(uncoordinatedV(:,2));
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
            ccg = ccg(:,1,1);
            Auto.uncoordinated.vRipples = [Auto.uncoordinated.vRipples , ccg];
            
            count = 0;
            y = ripplesV(:,2);%exclude coordinated + associated ripples
            for s = 1 : length(x)
                seed = x(s);
                tmp = sum(and(y > seed-0.2 , y < seed+0.2));
                if tmp>1
                    count = count + 1;
                end
            end
            BurstIndex.V.uncoordinated.all = [BurstIndex.V.uncoordinated.all ; count/length(x)];
            clear ccg time x y count
            
            %% VER ESTOOOOOOOO
%             %% Autoccorelogram coordionated and uncooridnated with subsampling
%             % Burst Index calculation as we discuss with GG
%             % coordinated dRipples
%             x = coordinated(:,2);
%             y = coordinated(:,2);
%             [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%             [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
%             ccg = ccg(:,1,1);
%             dRipples_coordinated_single_ds = [dRipples_coordinated_single_ds , ccg];
%             clear ccg time x y count tmp
%             
%             % coordinated vRipples
%             x = unique(coordinatedV(:,2));
%             y = unique(coordinatedV(:,2));
%             [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%             [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
%             ccg = ccg(:,1,1);
%             vRipples_coordinated_single_ds = [vRipples_coordinated_single_ds , ccg];
%             clear ccg time x y count tmp
%             
%             % uncoordinated dRipples
%             tmp = [];
%             for i = 1:50
%                 x = uncoordinated(randperm(length(uncoordinated)),2);
%                 x = x(1:length(coordinated));
%                 y = x;
%                 [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%                 [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
%                 ccg = ccg(:,1,1);
%                 tmp = [tmp , ccg];
%                 clear x y ccg time s ids groups
%             end
%             dRipples_uncoordinated_single_ds = [dRipples_uncoordinated_single_ds , mean(tmp,2)];
%             clear ccg time x y count tmp
%             
%             % uncoordinated vRipples
%             if length(coordinatedV) < length(uncoordinatedV)
%                 tmp = [];
%                 for i = 1:50
%                     x = uncoordinatedV(randperm(length(uncoordinatedV)),2);
%                     x = x(1:length(coordinatedV));
%                     y = x;
%                     [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%                     [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
%                     ccg = ccg(:,1,1);
%                     tmp = [tmp , ccg];
%                     clear x y ccg time s ids groups
%                 end
%             else
%                 tmp = [];
%                 for i = 1:50
%                     x = uncoordinatedV(randperm(length(uncoordinatedV)),2);
%                     y = x;
%                     [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%                     [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
%                     ccg = ccg(:,1,1);
%                     tmp = [tmp , ccg];
%                     clear x y ccg time s ids groups
%                 end
%             end
%             
%             vRipples_uncoordinated_single_ds = [vRipples_uncoordinated_single_ds , mean(tmp,2)];
%             clear ccg time x y count tmp

            %% Crosscorrelogram dRipples - vRipples
            x = ripplesD(:,2);
            y = ripplesV(:,2);
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
            ccg = ccg(:,1,2);
            Cross.DV.all = [Cross.DV.all , ccg]; clear x y s ids groups
            
%             %% Crosscorrelogram coordionated/uncooridnated with all ripples
%             % Burst Index calculation as we discuss with GG
%             % coordinated dRipples
%             x = coordinated(:,2);
%             y = ripplesD(:,2);
%             [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%             [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
%             ccg = ccg(:,1,2);
%             dRipples_coordinated_single_cross_with_all = [dRipples_coordinated_single_cross_with_all , ccg];
%             clear ccg time x y count tmp
%             
%             % coordinated vRipples
%             x = unique(coordinatedV(:,2));
%             y = ripplesV(:,2);
%             [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%             [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
%             ccg = ccg(:,1,2);
%             vRipples_coordinated_single_cross_with_all = [vRipples_coordinated_single_cross_with_all , ccg];
%             clear ccg time x y count tmp
%             
%             % uncoordinated dRipples
%             x = uncoordinated(:,2);
%             y = ripplesD(:,2);
%             [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%             [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
%             ccg = ccg(:,1,2);
%             dRipples_uncoordinated_single_cross_with_all = [dRipples_uncoordinated_single_cross_with_all , ccg];
%             clear ccg time x y count tmp
%             
%             % uncoordinated dRipples
%             x = uncoordinatedV(:,2);
%             y = ripplesV(:,2);
%             [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%             [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
%             ccg = ccg(:,1,2);
%             vRipples_uncoordinated_single_cross_with_all = [vRipples_uncoordinated_single_cross_with_all , ccg];
%             clear ccg time x y count tmp
%             
%             %% Crosscorrelogram coordionated vs uncooridnated
%             % coordinated dRipples
%             x = coordinated(:,2);
%             y = uncoordinated(:,2);
%             [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%             [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
%             ccg = ccg(:,1,2);
%             dRipples_coordinated_single_cross_with_uncoordinated = [dRipples_coordinated_single_cross_with_uncoordinated , ccg];
%             clear ccg time x y count tmp
%             
%             % coordinated vRipples
%             x = unique(coordinatedV(:,2));
%             y = unique(uncoordinatedV(:,2));
%             y = y(~ismember(y,x));
%             [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%             [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
%             ccg = ccg(:,1,2);
%             vRipples_coordinated_single_cross_with_uncoordinated = [vRipples_coordinated_single_cross_with_uncoordinated , ccg];
%             clear ccg time x y count tmp
%             
%             
%             %% Burst Index calculation as we discuss with GG
%             % per condition
%             % coordinated dRipples
%             x = coordinatedB(:,2);
%             count = 0;
%             y = ripplesD(:,2);
%             for s = 1 : length(x)
%                 seed = x(s);
%                 tmp = sum(and(y > seed-0.2 , y < seed+0.2));
%                 if tmp>1
%                     break
%                     count = count + 1;
%                 end
%             end
%             Burst_Index_cooridnatedB_dRipples = [Burst_Index_cooridnatedB_dRipples ; count/length(x)];
%             clear ccg time x y count
%             
%             x = coordinatedR(:,2);
%             count = 0;
%             y = ripplesD(:,2);
%             for s = 1 : length(x)
%                 seed = x(s);
%                 tmp = sum(and(y > seed-0.2 , y < seed+0.2));
%                 if tmp>1
%                     count = count + 1;
%                 end
%             end
%             Burst_Index_cooridnatedR_dRipples = [Burst_Index_cooridnatedR_dRipples ; count/length(x)];
%             clear ccg time x y count
%             
%             x = coordinatedA(:,2);
%             count = 0;
%             y = ripplesD(:,2);
%             for s = 1 : length(x)
%                 seed = x(s);
%                 tmp = sum(and(y > seed-0.2 , y < seed+0.2));
%                 if tmp>1
%                     count = count + 1;
%                 end
%             end
%             Burst_Index_cooridnatedA_dRipples = [Burst_Index_cooridnatedA_dRipples ; count/length(x)];
%             clear ccg time x y count
%             
%             
%             x = coordinatedB_V(:,2);
%             count = 0;
%             y = ripplesV(:,2);
%             for s = 1 : length(x)
%                 seed = x(s);
%                 tmp = sum(and(y > seed-0.2 , y < seed+0.2));
%                 if tmp>1
%                     count = count + 1;
%                 end
%             end
%             Burst_Index_cooridnatedB_vRipples = [Burst_Index_cooridnatedB_vRipples ; count/length(x)];
%             clear ccg time x y count
%             
%             x = coordinatedR_V(:,2);
%             count = 0;
%             y = ripplesV(:,2);
%             for s = 1 : length(x)
%                 seed = x(s);
%                 tmp = sum(and(y > seed-0.2 , y < seed+0.2));
%                 if tmp>1
%                     count = count + 1;
%                 end
%             end
%             Burst_Index_cooridnatedR_vRipples = [Burst_Index_cooridnatedR_vRipples ; count/length(x)];
%             clear ccg time x y count
%             
%             x = coordinatedA_V(:,2);
%             count = 0;
%             y = ripplesV(:,2);
%             for s = 1 : length(x)
%                 seed = x(s);
%                 tmp = sum(and(y > seed-0.2 , y < seed+0.2));
%                 if tmp>1
%                     count = count + 1;
%                 end
%             end
%             Burst_Index_cooridnatedA_vRipples = [Burst_Index_cooridnatedA_vRipples ; count/length(x)];
%             clear ccg time x y count
%             
%             % uncoordinated
%             x = uncoordinatedB(:,2);
%             count = 0;
%             y = ripplesD(:,2);
%             for s = 1 : length(x)
%                 seed = x(s);
%                 tmp = sum(and(y > seed-0.2 , y < seed+0.2));
%                 if tmp>1
%                     count = count + 1;
%                 end
%             end
%             Burst_Index_uncooridnatedB_dRipples = [Burst_Index_uncooridnatedB_dRipples ; count/length(x)];
%             clear ccg time x y count
%             
%             x = uncoordinatedR(:,2);
%             count = 0;
%             y = ripplesD(:,2);
%             for s = 1 : length(x)
%                 seed = x(s);
%                 tmp = sum(and(y > seed-0.2 , y < seed+0.2));
%                 if tmp>1
%                     count = count + 1;
%                 end
%             end
%             Burst_Index_uncooridnatedR_dRipples = [Burst_Index_uncooridnatedR_dRipples ; count/length(x)];
%             clear ccg time x y count
%             
%             x = uncoordinatedA(:,2);
%             count = 0;
%             y = ripplesD(:,2);
%             for s = 1 : length(x)
%                 seed = x(s);
%                 tmp = sum(and(y > seed-0.2 , y < seed+0.2));
%                 if tmp>1
%                     count = count + 1;
%                 end
%             end
%             Burst_Index_uncooridnatedA_dRipples = [Burst_Index_uncooridnatedA_dRipples ; count/length(x)];
%             clear ccg time x y count
%             
%             
%             x = uncoordinatedB_V(:,2);
%             count = 0;
%             y = ripplesV(:,2);
%             for s = 1 : length(x)
%                 seed = x(s);
%                 tmp = sum(and(y > seed-0.2 , y < seed+0.2));
%                 if tmp>1
%                     count = count + 1;
%                 end
%             end
%             Burst_Index_uncooridnatedB_vRipples = [Burst_Index_uncooridnatedB_vRipples ; count/length(x)];
%             clear ccg time x y count
%             
%             x = uncoordinatedR_V(:,2);
%             count = 0;
%             y = ripplesV(:,2);
%             for s = 1 : length(x)
%                 seed = x(s);
%                 tmp = sum(and(y > seed-0.2 , y < seed+0.2));
%                 if tmp>1
%                     count = count + 1;
%                 end
%             end
%             Burst_Index_uncooridnatedR_vRipples = [Burst_Index_uncooridnatedR_vRipples ; count/length(x)];
%             clear ccg time x y count
%             
%             x = uncoordinatedA_V(:,2);
%             count = 0;
%             y = ripplesV(:,2);
%             for s = 1 : length(x)
%                 seed = x(s);
%                 tmp = sum(and(y > seed-0.2 , y < seed+0.2));
%                 if tmp>1
%                     count = count + 1;
%                 end
%             end
%             Burst_Index_uncooridnatedA_vRipples = [Burst_Index_uncooridnatedA_vRipples ; count/length(x)];
%             clear ccg time x y count
%             
%             % All
%             x = Restrict(ripplesD(:,2),NREM_B);
%             count = 0;
%             y = ripplesD(:,2);
%             for s = 1 : length(x)
%                 seed = x(s);
%                 tmp = sum(and(y > seed-0.2 , y < seed+0.2));
%                 if tmp>1
%                     count = count + 1;
%                 end
%             end
%             Burst_Index_allB_dRipples = [Burst_Index_allB_dRipples ; count/length(x)];
%             clear ccg time x y count
%             
%             x = Restrict(ripplesD(:,2),NREM_R);
%             count = 0;
%             y = ripplesD(:,2);
%             for s = 1 : length(x)
%                 seed = x(s);
%                 tmp = sum(and(y > seed-0.2 , y < seed+0.2));
%                 if tmp>1
%                     count = count + 1;
%                 end
%             end
%             Burst_Index_allR_dRipples = [Burst_Index_allR_dRipples ; count/length(x)];
%             clear ccg time x y count
%             
%             x = Restrict(ripplesD(:,2),NREM_A);
%             count = 0;
%             y = ripplesD(:,2);
%             for s = 1 : length(x)
%                 seed = x(s);
%                 tmp = sum(and(y > seed-0.2 , y < seed+0.2));
%                 if tmp>1
%                     count = count + 1;
%                 end
%             end
%             Burst_Index_allA_dRipples = [Burst_Index_allA_dRipples ; count/length(x)];
%             clear ccg time x y count
%             
%             
%             x = Restrict(ripplesV(:,2),NREM_B);
%             count = 0;
%             y = ripplesV(:,2);
%             for s = 1 : length(x)
%                 seed = x(s);
%                 tmp = sum(and(y > seed-0.2 , y < seed+0.2));
%                 if tmp>1
%                     count = count + 1;
%                 end
%             end
%             Burst_Index_allB_vRipples = [Burst_Index_allB_vRipples ; count/length(x)];
%             clear ccg time x y count
%             
%             x = Restrict(ripplesV(:,2),NREM_R);
%             count = 0;
%             y = ripplesV(:,2);
%             for s = 1 : length(x)
%                 seed = x(s);
%                 tmp = sum(and(y > seed-0.2 , y < seed+0.2));
%                 if tmp>1
%                     count = count + 1;
%                 end
%             end
%             Burst_Index_allR_vRipples = [Burst_Index_allR_vRipples ; count/length(x)];
%             clear ccg time x y count
%             
%             x = Restrict(ripplesV(:,2),NREM_A);
%             count = 0;
%             y = ripplesV(:,2);
%             for s = 1 : length(x)
%                 seed = x(s);
%                 tmp = sum(and(y > seed-0.2 , y < seed+0.2));
%                 if tmp>1
%                     count = count + 1;
%                 end
%             end
%             Burst_Index_allA_vRipples = [Burst_Index_allA_vRipples ; count/length(x)];
%             clear ccg time x y count
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


% plot(tttt,mean(Auto.dRipples.all./max(Auto.dRipples.all),2),'k'),hold on
% plot(tttt,mean(Auto.vRipples.all./max(Auto.vRipples.all),2),'r'),xlim([-1 1])


%% Cross-Correlograms
% dorsal-ventral Ripples per condition
x = dRipples(:,2);
y = vRipples(:,2);
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',1,'mode','ccg');
ccg = ccg(:,1,2)%./sum(ccg(:,1,2));
figure,plot(tttt,ccg,'k','LineWidth',2),hold on
% ax = gca; % axes handle
% ax.YAxis.Exponent = 0;
xline(0,'--')

% dorsal-ventral Ripples per condition
x = dRipplesB(:,2);
y = vRipplesB(:,2);
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',2,'mode','ccg');
ccg = ccg(:,1,2)./sum(ccg(:,1,2));
figure,plot(tttt,ccg,'k','LineWidth',2),hold on
ax = gca; % axes handle
ax.YAxis.Exponent = 0;

% dorsal-ventral Ripples
x = dRipplesR(:,2);
y = vRipplesR(:,2);
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',2,'mode','ccg');
ccg = ccg(:,1,2)./sum(ccg(:,1,2));
plot(tttt,ccg,'b','LineWidth',2)
ax = gca; % axes handle
ax.YAxis.Exponent = 0;

% dorsal-ventral Ripples
x = dRipplesA(:,2);
y = vRipplesA(:,2);
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',2,'mode','ccg');
ccg = ccg(:,1,2)./sum(ccg(:,1,2));
plot(tttt,ccg,'r','LineWidth',2)
ax = gca; % axes handle
ax.YAxis.Exponent = 0;


%% Cross-Correlograms
% dorsal-ventral coordinated and uncoordinated Ripples
data = [dRipples_coordinatedB ; dRipples_coordinatedR ; dRipples_coordinatedA];
data1 = [dRipples_uncoordinatedB ; dRipples_uncoordinatedR ; dRipples_uncoordinatedA];
data2 = [vRipples_coordinatedB ; vRipples_coordinatedR ; vRipples_coordinatedA];
data3 = [vRipples_uncoordinatedB ; vRipples_uncoordinatedR ; vRipples_uncoordinatedA];
x = unique(data2(:,2));
y =  unique(data3(:,2));
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',2,'mode','ccg');
ccg = ccg(:,1,2)./sum(ccg(:,1,2));
figure,plot(tttt,ccg,'r','LineWidth',2),hold on
% ax = gca; % axes handle
% ax.YAxis.Exponent = 0;
xline(0,'--')

x = unique(data2(:,2));
y =  unique(data3(:,2));
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',2,'mode','ccg');
ccg = ccg(:,1,2)./sum(ccg(:,1,2));
plot(tttt,ccg,'r','LineWidth',2),hold on
% ax = gca; % axes handle
% ax.YAxis.Exponent = 0;
xline(0,'--')

% dorsal-ventral coordinated Ripples
x = unique(data(:,2));
y =  unique(data2(:,2));
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',2,'mode','ccg');
ccg = ccg(:,1,2)./sum(ccg(:,1,2));
figure,plot(tttt,ccg,'b','LineWidth',2),hold on
% ax = gca; % axes handle
% ax.YAxis.Exponent = 0;
xline(0,'--')


% dorsal-ventral coordinated Ripples per condition
data = [dRipples_uncoordinatedB];
data1 = [dRipples_uncoordinatedB];
x = unique(data(:,2));
y =  unique(data1(:,2));
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',1,'smooth',2,'mode','ccg');
ccg = ccg(:,1,1)./sum(ccg(:,1,1));
figure,plot(tttt,ccg,'k','LineWidth',2),hold on
% ax = gca; % axes handle
% ax.YAxis.Exponent = 0;
xline(0,'--')

% dorsal-ventral coordinated Ripples per condition
data = [dRipples_uncoordinatedR];
data1 = [dRipples_uncoordinatedR];
x = unique(data(:,2));
y =  unique(data1(:,2));
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',1,'smooth',2,'mode','ccg');
ccg = ccg(:,1,1)./sum(ccg(:,1,1));
plot(tttt,ccg,'b','LineWidth',2),hold on
% ax = gca; % axes handle
% ax.YAxis.Exponent = 0;
xline(0,'--')

% dorsal-ventral coordinated Ripples per condition
data = [dRipples_uncoordinatedA];
data1 = [dRipples_uncoordinatedA];
x = unique(data(:,2));
y =  unique(data1(:,2));
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',1,'smooth',2,'mode','ccg');
ccg = ccg(:,1,1)./sum(ccg(:,1,1));
plot(tttt,ccg,'r','LineWidth',2),hold on
% ax = gca; % axes handle
% ax.YAxis.Exponent = 0;
xline(0,'--')

%% IRI
tmp1 = dIRI;
tmp2 = vIRI;

[y,x]=histcounts(tmp1,100,'BinLimits',[0 1],'Normalization','probability')
% plot(x(2:end),Smooth(y,2),'r','LineWidth',2),xlim([0.015 0.1]),hold on
bar(x(2:end),Smooth(y,2),'FaceColor','r','FaceAlpha',0.5,'EdgeColor','none'),hold on
clear x y

[y,x]=histcounts(tmp2,100,'BinLimits',[0 1],'Normalization','probability')
% plot(x(2:end),Smooth(y,2),'b','LineWidth',2),xlim([0.015 0.1])
bar(x(2:end),Smooth(y,2),'FaceColor','b','FaceAlpha',0.5,'EdgeColor','none')%,xlim([0.015 0.1])
clear x y tmp1 tmp2


% IRI distributions per conditions
tmp1 = diff(vRipplesB(:,2));
% tmp1(tmp1>2) = [];
tmp2 = diff(vRipplesR(:,2));
% tmp2(tmp2>2) = [];
tmp3 = diff(vRipplesA(:,2));
% tmp3(tmp3>2) = [];

[y,x]=histcounts(tmp1,200,'BinLimits',[0 1],'Normalization','probability')
figure,plot(x(2:end),Smooth(y,1),'k','LineWidth',2),hold on
% bar(x(2:end),y,'FaceColor','k','FaceAlpha',0.5),hold on
clear x y

[y,x]=histcounts(tmp2,200,'BinLimits',[0 1],'Normalization','probability')
plot(x(2:end),Smooth(y,1),'b','LineWidth',2),hold on
% bar(x(2:end),y,'FaceColor','b','FaceAlpha',0.5),hold on
clear x y

[y,x]=histcounts(tmp3,200,'BinLimits',[0 1],'Normalization','probability')
plot(x(2:end),Smooth(y,1),'r','LineWidth',2)
% bar(x(2:end),y,'FaceColor','r','FaceAlpha',0.5),hold on
clear x y


% IRI distributions per conditions
tmp1 = diff(dRipplesB(:,2));
% tmp1(tmp1>2) = [];
tmp2 = diff(dRipplesR(:,2));
% tmp2(tmp2>2) = [];
tmp3 = diff(dRipplesA(:,2));
% tmp3(tmp3>2) = [];

[y,x]=histcounts(tmp1,200,'BinLimits',[0 1],'Normalization','probability')
figure,plot(x(2:end),Smooth(y,1),'k','LineWidth',2),hold on
% bar(x(2:end),y,'FaceColor','k','FaceAlpha',0.5),hold on
clear x y

[y,x]=histcounts(tmp2,200,'BinLimits',[0 1],'Normalization','probability')
plot(x(2:end),Smooth(y,1),'b','LineWidth',2),hold on
% bar(x(2:end),y,'FaceColor','b','FaceAlpha',0.5),hold on
clear x y

[y,x]=histcounts(tmp3,200,'BinLimits',[0 1],'Normalization','probability')
plot(x(2:end),Smooth(y,1),'r','LineWidth',2)
% bar(x(2:end),y,'FaceColor','r','FaceAlpha',0.5),hold on
clear x y

% Plot cummulative IRI graph
figure,boxplot([dIRIB dIRIR dIRIA vIRIB vIRIR vIRIA])
data = [dIRIB dIRIR dIRIA vIRIB vIRIR vIRIA];

%% Durations
tmp1 = dRipples(:,3)-dRipples(:,1);
tmp2 = vRipples(:,3)-vRipples(:,1);

[y,x]=histcounts(tmp1,100,'BinLimits',[0 0.1],'Normalization','probability')
% plot(x(2:end),Smooth(y,2),'r','LineWidth',2),xlim([0.015 0.1]),hold on
bar(x(2:end),Smooth(y,2),'FaceColor','r','FaceAlpha',0.5,'EdgeColor','none'),hold on
clear x y

[y,x]=histcounts(tmp2,100,'BinLimits',[0 0.1],'Normalization','probability')
% plot(x(2:end),Smooth(y,2),'b','LineWidth',2),xlim([0.015 0.1])
bar(x(2:end),Smooth(y,2),'FaceColor','b','FaceAlpha',0.5,'EdgeColor','none'),xlim([0.015 0.1])
clear x y tmp1 tmp2

%% Ripples cooridnated
data = [dRipples_coordinatedB dRipples_coordinatedR dRipples_coordinatedA vRipples_coordinatedB vRipples_coordinatedR vRipples_coordinatedA];
figure,boxplot(data)
[p,t,stats] = anova1(data);
[c,m,h,nms] = multcompare(stats);

%% burst
figure,boxplot([dBurstIndexB dBurstIndexR dBurstIndexA vBurstIndexB vBurstIndexR vBurstIndexA])

figure,boxplot([rateD rateV])


%% Burst Index of coordinated dorsal and ventral ripples (if is positive means that dRipples occure first)

data = [Burst_Index_uncooridnated_dRipples , Burst_Index_cooridnated_dRipples Burst_Index_uncooridnated_vRipples , Burst_Index_cooridnated_vRipples];
figure,boxplot(data)

% per condition
% just coordinated
data = [Burst_Index_cooridnatedB_dRipples , Burst_Index_cooridnatedR_dRipples , Burst_Index_cooridnatedA_dRipples Burst_Index_cooridnatedB_vRipples , Burst_Index_cooridnatedR_vRipples , Burst_Index_cooridnatedA_vRipples];
subplot(1,3,1),boxplot(data),ylim([0 0.8])
% anova1(data)

% just uncoordinated
data = [Burst_Index_uncooridnatedB_dRipples , Burst_Index_uncooridnatedR_dRipples , Burst_Index_uncooridnatedA_dRipples Burst_Index_uncooridnatedB_vRipples , Burst_Index_uncooridnatedR_vRipples , Burst_Index_uncooridnatedA_vRipples];
subplot(1,3,2),boxplot(data),ylim([0 0.8])
% anova1(data)


% just coordinated
data = [Burst_Index_allB_dRipples , Burst_Index_allR_dRipples , Burst_Index_allA_dRipples Burst_Index_allB_vRipples , Burst_Index_allR_vRipples , Burst_Index_allA_vRipples];
subplot(1,3,3),boxplot(data)
% anova1(data)

%% Correlograms session by session
% Autocorrelograms without subsampling
x = sum(dRipples_coordinated_single,2)./sum(sum(dRipples_coordinated_single));
y = sum(dRipples_uncoordinated_single,2)./sum(sum(dRipples_uncoordinated_single));
figure,subplot(4,2,1),plot([-1:0.01:1],x)
hold on
plot([-1:0.01:1],y),ylim([0 0.02])
clear x y

x = sum(vRipples_coordinated_single,2)./sum(sum(vRipples_coordinated_single));
y = sum(vRipples_uncoordinated_single,2)./sum(sum(vRipples_uncoordinated_single));
subplot(4,2,2),plot([-1:0.01:1],x)
hold on
plot([-1:0.01:1],y),ylim([0 0.02])

% Autocorrelograms session by session with subsampling
x = sum(dRipples_coordinated_single_ds,2)./sum(sum(dRipples_coordinated_single_ds));
y = sum(dRipples_uncoordinated_single_ds,2)./sum(sum(dRipples_uncoordinated_single_ds));
subplot(4,2,3),plot([-1:0.01:1],x)
hold on
plot([-1:0.01:1],y),ylim([0 0.02])
clear x y

x = sum(vRipples_coordinated_single_ds,2)./sum(sum(vRipples_coordinated_single_ds));
y = sum(vRipples_uncoordinated_single_ds,2)./sum(sum(vRipples_uncoordinated_single_ds));
subplot(4,2,4),plot([-1:0.01:1],x)
hold on
plot([-1:0.01:1],y),ylim([0 0.02])

% Correlograms session by session without subsampling
x = sum(dRipples_coordinated_single_cross_with_all,2)./sum(sum(dRipples_coordinated_single_cross_with_all));
y = sum(dRipples_uncoordinated_single_cross_with_all,2)./sum(sum(dRipples_uncoordinated_single_cross_with_all));
subplot(4,2,5),plot([-1:0.01:1],x)
hold on
plot([-1:0.01:1],y),ylim([0 0.02])
clear x y

x = sum(vRipples_coordinated_single_cross_with_all,2)./sum(sum(vRipples_coordinated_single_cross_with_all));
y = sum(vRipples_uncoordinated_single_cross_with_all,2)./sum(sum(vRipples_uncoordinated_single_cross_with_all));
subplot(4,2,6),plot([-1:0.01:1],x)
hold on
plot([-1:0.01:1],y),ylim([0 0.02])

%% Correlograms session by session with subsampling
x = sum(dRipples_coordinated_single_cross_with_uncoordinated,2)./sum(sum(dRipples_coordinated_single_cross_with_uncoordinated));
subplot(4,2,7),plot([-1:0.01:1],x), ylim([0 0.02])
clear x y

x = sum(vRipples_coordinated_single_cross_with_uncoordinated,2)./sum(sum(vRipples_coordinated_single_cross_with_uncoordinated));
subplot(4,2,8),plot([-1:0.01:1],x),ylim([0 0.02])
clear x y
