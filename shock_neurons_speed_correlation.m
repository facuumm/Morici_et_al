
%Variables
path = {'E:\Rat126\Ephys\in_Pyr';'E:\Rat103\usable';'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr\'};%List of folders from the path
criteria_fr = 0;
criteria_type = 0; % criteria for celltype (0:pyr, 1:int, 2:all)
win = 1;

%Storage
Speed = [];
dHPC.shock = [];     dHPC.no = [];
vHPC.shock = [];     vHPC.no = [];

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
        baselineTS = baselineTS./1000;                  aversiveTS = aversiveTS./1000;                    rewardTS = rewardTS./1000;
        aversiveTS_run = aversiveTS_run./1000;          rewardTS_run = rewardTS_run./1000;
        
        s = [];
        for i = 1 : size(Shocks_filt,1)
            tmp = Restrict(behavior.speed.aversive,[Shocks_filt(i)-40 Shocks_filt(i)+40]);
            if size(tmp,1)>=2400
                s = [s , tmp(1:2400,2)];
            end
        end
        %         figure,plot([-20 : 1/30:20-1/30],nanmean(s'))
        Speed = [Speed ; nanmean(s')];
        
        behavior.speed.aversive(InIntervals(behavior.speed.aversive(:,1),[Shocks_filt Shocks_filt+1]),2) = nan;
        
        T = [behavior.pos.aversive(1,1) : win : behavior.pos.aversive(end,1)+win];
        meanSpeed = [];
        for i = 1 : length(T)-1
            tmp = Restrict(behavior.speed.aversive,[T(i) T(i+1)]);
            meanSpeed = [meanSpeed ; nanmean(tmp(:,2))];
        end
        
        
        %% Spikes
        % Load Units
        disp('Uploading Spiking activity')
        cd 'Spikesorting'
        [clusters , numberD , numberV , spks , spks_dHPC , spks_vHPC , cellulartype] = load_SU_FM(cd,criteria_type,criteria_fr,aversiveTS_run,rewardTS_run);
        
        if exist('dHPC_responsivness_all.mat')
            load('dHPC_responsivness_all.mat')
%             clu = dHPC_resp.id(dHPC_resp.resp_ave==1);
            clu = dHPC_resp.id;
            dHPC_resp.Speedcorrelation.ShockSU = [];
            dHPC_resp.Speedcorrelation.noShockSU = [];
            if size(clu,1)>=1
                freq=1/win;
                spiketrains = [];
                speedCorrelation = [];
                for i=1:size(clu,1)
                    %                     time1 = sum(Shocks_filt+1 - Shocks_filt);
                    %                     x = Restrict(spks(spks(:,1)==clu(i),2),[Shocks_filt Shocks_filt+1]);
                    %                     x = size(x,1);
                    %                     
                    %                     y = InvertIntervals([Shocks_filt Shocks_filt-1] , behavior.pos.aversive(1,1) , behavior.pos.aversive(end,1));
                    %                     time2 = sum(y(:,2)-y(:,1));
                    %                     y = Restrict(spks(spks(:,1)==clu(i),2),y);
                    %                     y = size(y,1); y = y/time2;
                    %
                    %                     if y > 0
                    [spiketrain,bins]=binspikes(spks(spks(:,1)==clu(i),2),freq,[behavior.pos.aversive(1,1) : 1 : behavior.pos.aversive(end,1)]);
                    spiketrain = spiketrain(1:size(meanSpeed,1),:);
                    f = fitlm(meanSpeed,spiketrain);
                    
%                     if f.Rsquared.Ordinary>0.1
%                         figure,
%                         subplot(211),scatter((meanSpeed),(spiketrain),'filled','r'), hold on
%                         legend(['R2=',num2str(f.Rsquared.Ordinary) , ' / p=',num2str(coefTest(f))],'Location','Best')
%                         subplot(212),plot(f)
%                     end
                    
                    dHPC.shock = [dHPC.shock ];
                    speedCorrelation = [speedCorrelation ; f.Rsquared.Ordinary coefTest(f) clu(i) dHPC_resp.resp_ave(dHPC_resp.id==clu(i))];
                    %                     else
                    %                         dHPC.shock = [dHPC.shock ; nan nan];
                    %                     end
                end
                dHPC_resp.Speedcorrelation.all = speedCorrelation;
            end
            
%             save([cd,'\dHPC_responsivness_all.mat'],'dHPC_resp')
        end
        
        if exist('vHPC_responsivness_all.mat')
            load('vHPC_responsivness_all.mat')
%             clu = vHPC_resp.id(vHPC_resp.resp_ave==1);
            clu = vHPC_resp.id;
            vHPC_resp.Speedcorrelation.ShockSU = [];
            vHPC_resp.Speedcorrelation.noShockSU = [];
            if size(clu,1)>=1
                freq=1/win;
                spiketrains = [];
                speedCorrelation = [];
                for i=1:size(clu,1)
                    %                     time1 = sum(Shocks_filt+1 - Shocks_filt);
                    %                     x = Restrict(spks(spks(:,1)==clu(i),2),[Shocks_filt Shocks_filt+1]);
                    %                     x = size(x,1);
                    %                     
                    %                     y = InvertIntervals([Shocks_filt Shocks_filt-1] , behavior.pos.aversive(1,1) , behavior.pos.aversive(end,1));
                    %                     time2 = sum(y(:,2)-y(:,1));
                    %                     y = Restrict(spks(spks(:,1)==clu(i),2),y);
                    %                     y = size(y,1); y = y/time2;
                    %                     
                    %                     if y > 0
                    [spiketrain,bins]=binspikes(spks(spks(:,1)==clu(i),2),freq,[behavior.pos.aversive(1,1) : 1 : behavior.pos.aversive(end,1)]);
                    spiketrain = spiketrain(1:size(meanSpeed,1),:);
                    f = fitlm(meanSpeed,(spiketrain));
                    vHPC.shock = [vHPC.shock ; f.Rsquared.Ordinary coefTest(f)];
                    
                    if f.Rsquared.Ordinary>0.05
                        figure
                        subplot(211),scatter((meanSpeed),(spiketrain),'filled','b'), hold on
                        legend(['R2=',num2str(f.Rsquared.Ordinary) , ' / p=',num2str(coefTest(f))],'Location','Best')
                        subplot(212),plot(f)
                    end
                    speedCorrelation = [speedCorrelation ; f.Rsquared.Ordinary coefTest(f) clu(i) vHPC_resp.resp_ave(vHPC_resp.id==clu(i))];
                    %                     else
                    %                         vHPC.shock = [vHPC.shock ; nan nan];
                    %                     end
                end
                vHPC_resp.Speedcorrelation.all = speedCorrelation;
            end
            
%             save([cd,'\vHPC_responsivness_all.mat'],'vHPC_resp')
        end
        
    end
end

figure,
m = nanmean(Speed);
s = nansem(Speed);
T = [-40 : 1/30:40-1/30];
plot(T,m,'k'),hold on
ciplot(m-s , m+s , T , 'k'), alpha 0.5
xline(0,'--'),xline(1,'--')
xlim([-10 40])
ylim([2 16])


figure
subplot(121),bar([(sum(vHPC.shock(:,2)<0.05)/length(vHPC.shock))*100;100],'stacked')
subplot(122),bar([(sum(dHPC.shock(:,2)<0.05)/length(dHPC.shock))*100;100],'stacked')

figure
subplot(121),bar([(sum(vHPC.no(:,2)<0.05)/length(vHPC.no))*100;100],'stacked')
subplot(122),bar([(sum(dHPC.no(:,2)<0.05)/length(dHPC.no))*100;100],'stacked')


figure
histogram(vHPC(vHPC(:,2)<0.05 , 1),30,'Normalization','probability')

