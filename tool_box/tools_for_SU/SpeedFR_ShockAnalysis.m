function [dHPC , vHPC] = SpeedFR_ShockAnalysis(path)
% This function correlates speed and Firing Rate, save the R2 of that
% linear regression. Also save the Response during the Shock.
%
% Syntax: [dHPC , vHPC] = SpeedFR_ShockAnalysis(path)
%
% --- INPUTS ---
% path: structure, it contains the path where is located each session to
%       analyze.
%
% --- OUTPUT ----
% dHPC and vHPC: structur, they contain a matrix where I saved the R2 of
%                the linear regression between speed and firing Rate. Each
%                row is a neuron.
%                Rsquared / Pearson_Coeff / Shock_Response
%
% Morici Juan Facundo 10/2024

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
            tmp = Restrict(behavior.speed.aversive,[Shocks_filt(i)-2 Shocks_filt(i)+2]);
            if length(tmp) >= 120
                s = [s , tmp(1:120,2)];
            end
        end
        speed = nanmean(s'); clear s
        
        sub = [];
        count = 0;
        for i = 1 : 40
            if i == 1
            sub = speed(i);
            tmp = 3;
            else
            sub = [sub , nanmean(speed(tmp : tmp+3))];
            tmp = tmp+3;
            end
        end
        
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
            clu = dHPC_resp.id;%(dHPC_resp.resp_ave==1);
            if size(clu,1)>=1
                freq=1/win;
                spiketrains = [];
                speedCorrelation = [];
                for i=1:size(clu,1)
                    [~ , m] = min(abs([-2:0.1:2]-0));
                    [~ , mm] = min(abs([-2:0.1:2]-1));
                    
                    [resp , x] = max(dHPC_resp.curve_ave(i,m:mm));
%                     [resp] = mean(dHPC_resp.curve_ave(i,m:mm));
                    s = sub(x+m);
                    
                    [spiketrain,bins]=binspikes(spks(spks(:,1)==clu(i),2),freq,[behavior.pos.aversive(1,1) behavior.pos.aversive(end,1)]);
                    spiketrain = spiketrain(1:size(meanSpeed,1),:);
                    f = fitlm(meanSpeed,spiketrain);
                    p = coefTest(f);
                    z = dHPC_resp.resp_ave(i);
                    
                    dHPC.shock = [dHPC.shock ; f.Rsquared.Ordinary coefTest(f) resp s z p]; clear m mm resp s
                end
            end
        end
        
        if exist('vHPC_responsivness_all.mat')
            load('vHPC_responsivness_all.mat')
            clu = vHPC_resp.id;%(vHPC_resp.resp_ave==1);
            if size(clu,1)>=1
                freq=1/win;
                spiketrains = [];
                for i=1:size(clu,1)
                    
                    [~ , m] = min(abs([-2:0.1:2]-0));
                    [~ , mm] = min(abs([-2:0.1:2]-1));
                    
                    [resp , x] = max(vHPC_resp.curve_ave(i,m:mm));
%                     [resp] = mean(vHPC_resp.curve_ave(i,m:mm));
                    s = sub(x+m);
                    
                    [spiketrain,bins]=binspikes(spks(spks(:,1)==clu(i),2),freq,[behavior.pos.aversive(1,1) behavior.pos.aversive(end,1)]);
                    spiketrain = zscore(spiketrain(1:size(meanSpeed,1),:));
                    f = fitlm(meanSpeed,(spiketrain));
                    p = coefTest(f);
                    z = vHPC_resp.resp_ave(i);
                    vHPC.shock = [vHPC.shock ; f.Rsquared.Ordinary coefTest(f) resp s z p]; clear m mm resp s p

                end
            end
        end
    end
end
%  
% figure
% subplot(121)
% v = vHPC.shock;
% v(90,:) = [];
% fitlm(v(:,4) , v(:,3))
% plot(ans) , ylim([-0.5 4.5]) , xlim([0 70])
% subplot(122)
% scatter(v(:,4) , v(:,3),'filled') , ylim([-0.5 4.5]) , xlim([0 70])
% 
% 
% figure
% subplot(121)
% d = dHPC.shock;
% d(220,:) = [];
% fitlm(d(:,4) , d(:,3))
% plot(ans) , ylim([-0.5 4.5]) , xlim([0 70])
% subplot(122)
% scatter(d(:,4) , d(:,3),'filled') , ylim([-0.5 4.5]) , xlim([0 70])
% 
% y = [dHPC.shock(:,1) ; vHPC.shock(:,1)];
% x = [ones(length(dHPC.shock(:,1)),1) ; ones(length(vHPC.shock(:,1)),1)*2];
% figure
% boxplot(y,x)
% 
% figure,
% scatter(x,y,[],"filled",'jitter','on', 'jitterAmount',0.1),xlim([0 3]),hold on
% scatter([1 2] , [nanmean(dHPC.shock(:,1)) ; nanmean(vHPC.shock(:,1))],'filled')
% 
% [h p] = kstest2(dHPC.shock(:,1) , vHPC.shock(:,1))
% 
% figure
% cdfplot(dHPC.shock(:,1)),hold on
% cdfplot(vHPC.shock(:,1)),hold on
% % xlim([0 0.4])
% set(gca, 'XScale', 'log')
% 
% 
% 
% subsampling
R = [];
for j = 1 : 1000
i = vHPC.shock(:,6)<0.05;
ii = vHPC.shock(:,5)==1;
iii = and(i,ii);
x = vHPC.shock(iii==1,1);
y = vHPC.shock(iii,3); y = y-min(y); y = y/max(y);
% x(y==1) = [];
% y(y==1)=[];  y = y-min(y); y = y/max(y);
% x(y==1) = [];
% y(y==1)=[];  y = y-min(y); y = y/max(y);
x = vHPC.shock(iii,1);
y = vHPC.shock(iii,3); y = y-min(y); y = y/max(y);
x = x(not(isoutlier(y,"percentiles",[0 95]))); 
y = y(not(isoutlier(y,"percentiles",[0 95]))); y = y-min(y); y = y/max(y);
n = length(y);

i = dHPC.shock(:,6)<0.05;
ii = dHPC.shock(:,5)==1;
iii = and(i,ii);
x = dHPC.shock(iii,1);
y = dHPC.shock(iii,3); y = y-min(y); y = y/max(y);
x = x(not(isoutlier(y,"percentiles",[0 95]))); 
y = y(not(isoutlier(y,"percentiles",[0 95]))); y = y-min(y); y = y/max(y);
nn = randperm(length(x));
x = x(nn,:);    y = y(nn,:);
x = x(1:n);     y = y(1:n);

mdl3 = fitlm(x,y);
p = coefTest(mdl3);
R = [R ; mdl3.Rsquared.Ordinary p];
end

figure
[h hh] = kstest2(R(:,1),mdl1.Rsquared.Ordinary)
histogram(R(:,1),20,'Normalization','probability')
hold on,xline(mdl1.Rsquared.Ordinary)

% figure
% [h hh] = kstest2(R(:,2),coefTest(mdl1))
% histogram(R(:,2),50,'Normalization','probability')
% hold on,xline(coefTest(mdl1))
% cdfplot(R(:,2)),hold on,xline(0.05)



figure,
i = dHPC.shock(:,6)<0.05;
ii = dHPC.shock(:,5)==1;
iii = and(i,ii);
x = dHPC.shock(iii,1);
y = dHPC.shock(iii,3); y = y-min(y); y = y/max(y);
x = x(not(isoutlier(y,"percentiles",[0 95]))); 
y = y(not(isoutlier(y,"percentiles",[0 95]))); y = y-min(y); y = y/max(y);
% x(y==1) = [];
% y(y==1)=[];  y = y-min(y); y = y/max(y);
subplot(121),scatter(x,y,'filled'),xlim([0 0.3])
mdl1 = fitlm(x,y);


i = vHPC.shock(:,6)<0.05;
ii = vHPC.shock(:,5)==1;
iii = and(i,ii);
x = vHPC.shock(iii==1,1);
y = vHPC.shock(iii,3); y = y-min(y); y = y/max(y);
% x(y==1) = [];
% y(y==1)=[];  y = y-min(y); y = y/max(y);
% x(y==1) = [];
% y(y==1)=[];  y = y-min(y); y = y/max(y);
x = vHPC.shock(iii,1);
y = vHPC.shock(iii,3); y = y-min(y); y = y/max(y);
x = x(not(isoutlier(y,"percentiles",[0 95]))); 
y = y(not(isoutlier(y,"percentiles",[0 95]))); y = y-min(y); y = y/max(y); = fitlm(x,y);
subplot(122),scatter(x,y,'filled'),xlim([0 0.15])

figure
subplot(121),
h = plot(mdl1);
delete(h(1))
xlim([0 0.3]), ylim([0 1])

subplot(122),
h = plot(mdl2);
delete(h(1))
xlim([0 0.15]), ylim([0 1])

end