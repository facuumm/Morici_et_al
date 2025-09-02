clear
clc
close all

%% Parameters
path = {'E:\Rat103\usable';'E:\Rat126\Ephys\in_Pyr';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path
minimal_speed = 7; % minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods
Re = [];
Av = [];
meanSpeed = [];

curveA = []; curveR = [];
quietA = []; quietR = [];
speedA = []; speedR = [];
Rew = [];
N = 0;
%% Main loop, to iterate across sessions
for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    num_assembliesR = [];
    num_assembliesA = [];
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
        clear y A R
        
        % Defining what condition was first
        if aversiveTS_run(1) < rewardTS_run(1)
            config = 1;
        else
            config = 2;
        end
        
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
        
        Rew = [Rew ; size(Rewards_filt,1)];
        
        % periods of movment during eacj condition
        if rewardTS_run(1) < aversiveTS_run(1)
            load('laps1.mat','posx','posy');
            [camaraR,~] = find((camara(:,1)-rewardTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
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
        
        
        binsX = 400;
        binsY = 200;
        
        figure
        subplot(121),plot(behavior.pos.aversive(:,2),behavior.pos.aversive(:,3))
        subplot(122),plot(behavior.pos.reward(:,2),behavior.pos.reward(:,3))
        
        figure,
        X = ZeroToOne(behavior.pos.aversive(:,2));
        Y = ZeroToOne(behavior.pos.aversive(:,3));
        subplot(222)
        plot(X,Y)
        
        X = Bin(X , [0 1] , binsX);
        Y = Bin(Y , [0 1] , binsY);
        map.x = linspace(0,1,binsX);
        map.y = linspace(0,1,binsY);
%         map.count = Accumulate([x y],n,nBins);
        map.time = Accumulate([X Y],(1/30),[binsX binsY]);
        valid = map.time > (1/30);
        map.time = Smooth(Interpolate2(map.x,map.y,map.time,valid,'interpolate',5),0,'type','ll')';
        subplot(221)
        PlotColorMap(map.time),caxis([0 5])
        
        X = ZeroToOne(behavior.pos.reward(:,2));
        Y = ZeroToOne(behavior.pos.reward(:,3));
        subplot(224)
        plot(X,Y)
        
        X = Bin(X , [0 1] , binsX);
        Y = Bin(Y , [0 1] , binsY);
        map.x = linspace(0,1,binsX);
        map.y = linspace(0,1,binsY);
%         map.count = Accumulate([x y],n,nBins);
        map.time = Accumulate([X Y],(1/30),[binsX binsY]);
        valid = map.time > (1/30);
        map.time = Smooth(Interpolate2(map.x,map.y,map.time,valid,'interpolate',5),0,'type','ll')';
        subplot(223)
        PlotColorMap(map.time),caxis([0 5])
        
        
        
        
        
        % Restriction to quiet periods during AVERSIVE
        x = InIntervals(behavior.speed.aversive(:,1),behavior.quiet.aversive);  
        % Histogram of active periods 
        [y e] = histcounts(behavior.pos.aversive(not(x),2),60);
        y = y* (1/30);
        curveA = [curveA ; y]; clear y
        
        x = InIntervals(behavior.speed.aversive(:,1),behavior.quiet.aversive);  
        % Histogram of active periods 
        [y e] = histcounts(behavior.pos.aversive((x),2),60);
        y = y* (1/30);
        quietA = [quietA ; y]; clear y
        
        % average speed calculation
        [h,c] = histcounts(behavior.speed.aversive(not(x),2),[0:1:100]);
        m = mean(behavior.speed.aversive(not(x),2));
        m1 = mean(behavior.speed.aversive(:,2));
        m2 = (sum(not(x))/length(x))*100;
        Av = [Av ; h./sum(h)];
        % Average speed across positions
        [Y , E] = discretize(behavior.pos.aversive(:,2),[0:3:200]);
        count = [];
        for i = 1:length(E)
            count = [count ; mean(behavior.speed.aversive(Y==i,2),'omitnan')];
        end
        speedA = [speedA ; count'];
        clear x h Y E count
        
        % Restriction to quiet periods during REWARD
        x = InIntervals(behavior.speed.reward(:,1),behavior.quiet.reward); 
        % Histogram of active periods 
        [y e] = histcounts(behavior.pos.reward(not(x),2),60);
        y = y* (1/30);
        curveR = [curveR ; y]; clear y
        
        x = InIntervals(behavior.speed.reward(:,1),behavior.quiet.reward); 
        % Histogram of active periods 
        [y e] = histcounts(behavior.pos.reward((x),2),60);
        y = y* (1/30);
        quietR = [quietR ; y]; clear y
        
        % average speed calculation
        [h,c] = histcounts(behavior.speed.reward(not(x),2),[0:1:100]);
        mm = mean(behavior.speed.reward(not(x),2));
        mm1 = mean(behavior.speed.reward(:,2));
        mm2 = (sum(not(x))/length(x))*100;
        Re = [Re ; h./sum(h)];
        % Average speed across positions
        [Y , E] = discretize(behavior.pos.reward(:,2),[0:3:200]);
        count = [];
        for i = 1:length(E)
            count = [count ; mean(behavior.speed.reward(Y==i,2),'omitnan')];
        end
        speedR = [speedR ; count'];
        clear x h Y count
        
        meanSpeed = [meanSpeed ; m mm m1 mm1 m2 mm2];
        
        save([cd,'\behavioral_dataVF.mat'],'Rewards_filt','Shocks_filt','behavior')

        
        disp(' ')
        clear A aversiveTS aversiveTS_run baselineTS rewardTS rewardTS_run m mm m1 mm1 m2 mm2
        clear behavior bins Cell_type_classification cellulartype cond
        clear is K Kinfo group_dHPC group_vHPC matrixC matrixD matrixV
        clear NREM REM WAKE segmentation tmp cond
        clear spiketrains_dHPC spiketrains_vHPC opts MUA
        clear patterns Thresholded i  ii numberD numberV movement cross crossN
        clear Spikes bins Clusters Shocks_filt Rewards_filt config n_SU_D n_SU_V
        clear clusters coordinated coordinated_ripple_bursts coordinatedV
        clear cooridnated_event coordinatedV_refined coordinatedV_refined
    	N = N+1;
            
    end
    if tt < length(path)
        curveA = [curveA ; nan(3,size(curveA,2))];
        curveR = [curveR ; nan(3,size(curveR,2))];
        
        quietA = [quietA ; nan(3,size(quietA,2))];
        quietR = [quietR ; nan(3,size(quietR,2))];
        
        speedA = [speedA ; nan(3,size(speedA,2))];
        speedR = [speedR ; nan(3,size(speedR,2))];
    end

end

%% Plot Ocupancy
figure
subplot(121),imagesc([3:3:180],[1:1:size(curveA,1)],curveR)
subplot(122),imagesc([3:3:180],[1:1:size(curveA,1)],curveA)

figure
plot(e(2:end),mean(curveR,'omitnan')),ylim([0 24]),xlim([min(e(2:end)) max(e(2:end))])
hold on
plot(e(2:end),mean(curveA,'omitnan')),ylim([0 24]),xlim([min(e(2:end)) max(e(2:end))])

%% Plot average speed across space
figure
subplot(121),imagesc(E,[1:1:size(speedR,1)],speedR)
subplot(122),imagesc(E,[1:1:size(speedA,1)],speedA)

figure
m = mean(speedR,'omitnan');
plot(E,m,'b'),ylim([0 65]),xlim([min(E) max(E)]),hold on
s = std(speedR,'omitnan')./sqrt(N);
ciplot(m-s , m +s , E), alpha 0.5

m = mean(speedA,'omitnan');
plot(E,m,'r'),ylim([0 65]),xlim([min(E) max(E)]),hold on
s = std('omitnan')./sqrt(N);
ciplot(m-s , m +s , E,'r'), alpha 0.5

% Plot average speed
figure
m = mean(Re);
s = std(Re)./sqrt(size(Re,1));
plot(c(1:end-1),m,'b'),hold on
ciplot(m-s,m+s,c(1:end-1),'b'),alpha 0.1
xlim([0 100]),ylim([0 0.1])

m = mean(Av);
s = std(Av)./sqrt(size(Av,1));
plot(c(1:end-1),m,'r'),hold on
ciplot(m-s,m+s,c(1:end-1),'r'),alpha 0.1


%% Plot Ocupancy
figure
subplot(121),imagesc([3:3:180],[1:1:size(quietR,1)],quietR)
subplot(122),imagesc([3:3:180],[1:1:size(quietA,1)],quietA)

figure
plot(e(2:end),mean(quietR,'omitnan'),'b'),ylim([0 150]),xlim([min(e(2:end)) max(e(2:end))])
hold on
plot(e(2:end),mean(quietA,'omitnan'),'r'),ylim([0 150]),xlim([min(e(2:end)) max(e(2:end))])


%% Plot mean speed, percentage of Quiet and Active periods
figure
subplot(142)
scatter([ones(size(meanSpeed,1),1);ones(size(meanSpeed,1),1)*2] , [meanSpeed(:,1);meanSpeed(:,2)],'filled','MarkerEdgeColor','none','MarkerFaceColor','k'),xlim([0 3]),ylim([0 60])
hold on
scatter([1;2] , [mean(meanSpeed(:,1));mean(meanSpeed(:,2))],'filled','MarkerEdgeColor','none','MarkerFaceColor','r'),,xlim([0 3]),ylim([0 60])

subplot(141)
scatter([ones(size(meanSpeed,1),1);ones(size(meanSpeed,1),1)*2] , [meanSpeed(:,3);meanSpeed(:,4)],'filled','MarkerEdgeColor','none','MarkerFaceColor','k'),,xlim([0 3]),ylim([0 60])
hold on
scatter([1;2] , [mean(meanSpeed(:,3));mean(meanSpeed(:,4))],'filled','MarkerEdgeColor','none','MarkerFaceColor','r'),xlim([0 3]),ylim([0 60])

[h p] = ttest2(meanSpeed(:,1),meanSpeed(:,2));
[h p] = ttest2(meanSpeed(:,3),meanSpeed(:,4));

subplot(143)
scatter([ones(size(meanSpeed,1),1);ones(size(meanSpeed,1),1)*2] , [meanSpeed(:,5);meanSpeed(:,6)],'filled','MarkerEdgeColor','none','MarkerFaceColor','k'),,xlim([0 3]),ylim([0 100])
hold on
scatter([1;2] , [mean(meanSpeed(:,5));mean(meanSpeed(:,6))],'filled','MarkerEdgeColor','none','MarkerFaceColor','r'),,xlim([0 3]),ylim([0 100])

subplot(144)
scatter([ones(size(meanSpeed,1),1);ones(size(meanSpeed,1),1)*2] , [100-meanSpeed(:,5);100-meanSpeed(:,6)],'filled','MarkerEdgeColor','none','MarkerFaceColor','k'),,xlim([0 3]),ylim([0 100])
hold on
scatter([1;2] , [100-mean(meanSpeed(:,5));100-mean(meanSpeed(:,6))],'filled','MarkerEdgeColor','none','MarkerFaceColor','r'),,xlim([0 3]),ylim([0 100])


[h p] = ttest2(meanSpeed(:,5),meanSpeed(:,6));
