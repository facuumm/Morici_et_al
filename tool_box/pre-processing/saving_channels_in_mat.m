clear
clc
close all

%Parameters
path = 'D:\ANR\controls';

%Channles for ripple detection
% 22
% 10
dorsal = [93 57 33 43 17 23];
dorsal = [0 1 24 5 21];
% ventral = [53];
% dorsal = 1
num = 35; %number of channels including accelorometers

fs = 1250; %Sampling frequency
dt = 1/fs;

%List of folders from the path
files = dir(path);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);

final = cell(length(subFolders)-2,3);
clear files dirFlags
D = 1;
V = 1;
%% Load LFP
for t = 1 : length(subFolders)-2
    session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
    cd(session)
    x = dir([cd,'\*.lfp']);
    
%     %dHPC
%     tmp = (LoadBinary([x.folder,'\',x.name],'channels',dorsal(1),'frequency',fs,'nChannels',num))*0.195;
%     Time = [1/fs:1/fs:length(tmp)*1/fs]'; %Construct time vector
%     dHPC = [Time tmp];
    
    %vHPC1
    tmp = (LoadBinary([x.folder,'\',x.name],'channels',dorsal(2),'frequency',fs,'nChannels',num))*0.195;
    Time = [1/fs:1/fs:length(tmp)*1/fs]'; %Construct time vector
    vHPC1 = [Time tmp];
    
    %vHPC2
    tmp = (LoadBinary([x.folder,'\',x.name],'channels',dorsal(3),'frequency',fs,'nChannels',num))*0.195;
    vHPC2 = [Time tmp];
    
    %vHPC3
    tmp = (LoadBinary([x.folder,'\',x.name],'channels',dorsal(4),'frequency',fs,'nChannels',num))*0.195;
    vHPC3 = [Time tmp];
    
    %vHPC4
    tmp = (LoadBinary([x.folder,'\',x.name],'channels',dorsal(5),'frequency',fs,'nChannels',num))*0.195;
    vHPC4 = [Time tmp];
    
%     %vHPC5
%     tmp = (LoadBinary([x.folder,'\',x.name],'channels',dorsal(6),'frequency',fs,'nChannels',num))*0.195;
%     vHPC5 = [Time tmp];
    
%     save([cd,'\lfp1.mat'],'vHPC1','vHPC2','vHPC3','vHPC4','vHPC5') 
    save([cd,'\lfp.mat'],'vHPC1','vHPC2','vHPC3','vHPC4') 

    t
    figure,
    subplot(6,1,1),plot(dHPC(:,1),dHPC(:,2))
    subplot(6,1,2),plot(vHPC1(:,1),vHPC1(:,2))
    subplot(6,1,3),plot(vHPC2(:,1),vHPC2(:,2))
    subplot(6,1,4),plot(vHPC3(:,1),vHPC3(:,2))
    subplot(6,1,5),plot(vHPC4(:,1),vHPC4(:,2))
    subplot(6,1,6),plot(vHPC5(:,1),vHPC5(:,2))
%     
    clear dHPC vHPC1 vHPC2 vHPC3 vHPC4 vHPC5 Time
end

clear aversiveTS baselineTS channels_to_use dHPC_REM dHPC_SWS dorsal dt final freq
clear fs I M num path session states subFolders sws Rem t TxGcomod TxGphases
clear rewardTS Spectrum 