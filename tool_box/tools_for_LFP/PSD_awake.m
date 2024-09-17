function [power , f] = PSD_awake(path)
% Calculates PSD for dorsal and/or ventral lfp signal during the entier
% session, or during periods of movment.
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% OUTPUTS
% power: structure, it contains the spectrogram of each session.
%
% f: vector, it contains the frecuency dimension.
%
% Morici Juan Facundo, 09/2024


power.dHPC.aversive.all = [];           power.dHPC.reward.all = [];
power.vHPC.aversive.all = [];           power.vHPC.reward.all = [];

power.dHPC.aversive.movement = [];      power.dHPC.reward.movement = [];
power.vHPC.aversive.movement = [];      power.vHPC.reward.movement = [];

% Main loop, to iterate across sessions
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
        
        % Load LFP
        if exist('lfp1.mat')
            load([cd,'\lfp1.mat'])
        else exist('lfp.mat')
            load([cd,'\lfp.mat'])
        end
        
        %Loading TS of the sessions
        disp('Uploading session time stamps')
        load('session_organization.mat')
        load('behavioral_data.mat')
        
        baselineTS = baselineTS./1000;                  aversiveTS = aversiveTS./1000;                    rewardTS = rewardTS./1000;
        aversiveTS_run = aversiveTS_run./1000;          rewardTS_run = rewardTS_run./1000;
        
        if exist('dHPC')
            
%             time = SubtractIntervals([dHPC(1,1) dHPC(end,1)], [Shocks_filt-0.5 Shocks_filt+1.5]);
%             dHPC = (Restrict(dHPC , time));
            
            [spectrogram,f,s] = MTSpectrum(Restrict(dHPC,aversiveTS_run),'frequency',1250,'range',[0 40]);
            power.dHPC.aversive.all = [power.dHPC.aversive.all , spectrogram']; clear spectrogram f s
            
            [spectrogram,f,s] = MTSpectrum(Restrict(dHPC,rewardTS_run),'frequency',1250,'range',[0 40]);
            power.dHPC.reward.all = [power.dHPC.reward.all , spectrogram']; clear spectrogram f s
            
            [spectrogram,f,s] = MTSpectrum(Restrict(dHPC,movement.aversive),'frequency',1250,'range',[0 40]);
            power.dHPC.aversive.movement = [power.dHPC.aversive.movement , spectrogram']; clear spectrogram f s
            
            [spectrogram,f,s] = MTSpectrum(Restrict(dHPC,movement.reward),'frequency',1250,'range',[0 40]);
            power.dHPC.reward.movement = [power.dHPC.reward.movement , spectrogram']; clear spectrogram s
        end
        
        if exist('vHPC1')
%             time = SubtractIntervals([vHPC1(1,1) vHPC1(end,1)], [Shocks_filt-0.5 Shocks_filt+1.5]);
%             vHPC1 = (Restrict(vHPC1 , time));
            
            [spectrogram,f,s] = MTSpectrum(Restrict(vHPC1,aversiveTS_run),'frequency',1250,'range',[0 40]);
            power.vHPC.aversive.all = [power.vHPC.aversive.all , spectrogram']; clear spectrogram f s
            
            [spectrogram,f,s] = MTSpectrum(Restrict(vHPC1,rewardTS_run),'frequency',1250,'range',[0 40]);
            power.vHPC.reward.all = [power.vHPC.reward.all , spectrogram']; clear spectrogram f s
            
            [spectrogram,f,s] = MTSpectrum(Restrict(vHPC1,movement.aversive),'frequency',1250,'range',[0 40]);
            power.vHPC.aversive.movement = [power.vHPC.aversive.movement , spectrogram']; clear spectrogram f s
            
            [spectrogram,f,s] = MTSpectrum(Restrict(vHPC1,movement.reward),'frequency',1250,'range',[0 40]);
            power.vHPC.reward.movement = [power.vHPC.reward.movement , spectrogram']; clear spectrogram s
        end
        
        clear dHPC vHPC1 vHPC2 vHPC3 vHPC4 vHPC5
        clear aversiveTS aversiveTS_run rewardTS rewardTS_run baselineTS
        clear behavior movement config Rewards_filt Shocks_filt minimal_speed minimal_speed_time
        clear segments
        
    end
end


figure
subplot(221)
m = nanmean(power.dHPC.aversive.all');
s = nansem(power.dHPC.aversive.all');
plot(f,m,'r'),hold on
ciplot(m-s , m+s , f , 'r'), alpha 0.1

m = nanmean(power.dHPC.reward.all');
s = nansem(power.dHPC.reward.all');
subplot(221)
plot(f,m,'b'),hold on
ciplot(m-s , m+s , f , 'b'), alpha 0.1

subplot(222)
m = nanmean(power.vHPC.aversive.all');
s = nansem(power.vHPC.aversive.all');
plot(f,m,'r'),hold on
ciplot(m-s , m+s , f , 'r'), alpha 0.1

m = nanmean(power.vHPC.reward.all');
s = nansem(power.vHPC.reward.all');
plot(f,m,'b'),hold on
ciplot(m-s , m+s , f , 'b'), alpha 0.1

subplot(223)
m = nanmean(power.dHPC.aversive.movement');
s = nansem(power.dHPC.aversive.movement');
plot(f,m,'r'),hold on
ciplot(m-s , m+s , f , 'r'), alpha 0.1

m = nanmean(power.dHPC.reward.movement');
s = nansem(power.dHPC.reward.movement');
plot(f,m,'b'),hold on
ciplot(m-s , m+s , f , 'b'), alpha 0.1 

subplot(224)
m = nanmean(power.vHPC.aversive.movement');
s = nansem(power.vHPC.aversive.movement');
plot(f,m,'r'),hold on
ciplot(m-s , m+s , f , 'r'), alpha 0.1

m = nanmean(power.vHPC.reward.movement');
s = nansem(power.vHPC.reward.movement');
plot(f,m,'b'),hold on
ciplot(m-s , m+s , f , 'b'), alpha 0.1

end