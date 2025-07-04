%% Instructions
% This script takes a path where are all your folders of recordings
% sessions. It take all the subfolders inside each one and take the x and Y
% positions to calculate velocity.
% Using the position it give the begening and end of each lap.

clear
clc
close all

%% Parameters
% change currentdir deppending of the rat to be analyzed
currentdir = ('E:\Rat165\in_pyr');
listing = dir(currentdir);

% Parameters
%Heatmaps
X = 60; % Number of bins to construct heatmaps
Y = 20; % Number of bins to construct heatmaps
mySmooth = 9; %sice of gaussian smooth for heatmaps
cmin = 0;
cmax = 10;

%Time
dt = 1/30; %1sec/30frames
tt = 120; %how many time befor the begining of each lap I want to plot

%Parameters for Gaussian
sigma = 20; % pick sigma value for the gaussian
gaussFilter = gausswin(6*sigma + 1)';
gaussFilter = gaussFilter / sum(gaussFilter); % normalize

%Pixel/cm ratio (21cm = 100px)
pixcorrection = 21/100;

%Platform limits to count the number of laps
xi = 20;
xf = 140;

%% Generate cells where data is going to be stored
store = cell(length(listing)-2,1); % storage for neck position (x and y)
count = cell(length(listing)-2,1); % storage for alternation number
index = cell(length(listing)-2,1); % storage for alternation session position (start and end)
v = cell(length(listing)-2,1); % storage for velocity

for j = 1 : length(listing)-2
    % Get a logical vector that tells which is a directory.
    subFolders =  dir([listing(j+2).folder,'\',listing(j+2).name]);
    dirFlags = [subFolders.isdir];
    % Extract only those that are directories.
    subFolders = subFolders(dirFlags);
    clear dirFlags
    
    store{j,1} = cell(length(subFolders)-2,1); % storage for neck position (x and y)
    count{j,1} = cell(length(subFolders)-2,1); % storage for alternation number
    index{j,1} = cell(length(subFolders)-2,1); % storage for alternation session position (start and end)
    v{j,1} = cell(length(subFolders)-2,1); % storage for velocity
    
    % LAPS = [];
    % NUMLAPS = [];
    % POSX = [];
    % POSY = [];
    % VELOCITY = [];
    session = 1;
    
    for i=1:length(subFolders)-2
        if or(i == 2, i == 4)
%             if ~contains(subFolders(i+2).name,'Spikesorting')
                % open the data file for each video
                %     tmp0 = dir(fullfile([subFolders(i+2).folder,'\',subFolders(i+2).name],'*.csv')); %if labeling video is not in the path
                
                tmp0 = dir(fullfile([subFolders(i+2).folder,'\',subFolders(i+2).name],'*.mp4'));
                if ~isempty(tmp0)
                    nbvideos = length(tmp0)/2;
                    
                    %     nbvideos = 1;%if labeling video is not in the path
                    tmp = cell(1,nbvideos);
                    
                    tmp{1,1} = [tmp0(2).folder,'\',tmp0(2).name(1:end-12),'.csv'];
                    tmp{1,1} = readtable(tmp{1,1},'HeaderLines',3);
                    %     tmp{1,1} = readtable([tmp0(2).folder,'\',tmp0(2).name],'NumHeaderLines',3);%if labeling video is not in the path
                    
                    
                    store{j,1}{i,1} =  tmp{1}{:,5:7}; %store the position of the animal
                    
                    tmpo = store{j,1}{i,1}(:,1)*pixcorrection;
                    [id,t] = find(or(tmpo(:,1)<=xi, tmpo(:,1)>=xf));
                    z = diff(id);
                    b = [];
                    for k = 1:length(z)
                        dif = abs(tmpo(id(k))-tmpo(id(k+1)));
                        if (z(k)~=1) && dif > 100
                            b = [b;id(k-1),id(k+1)];
                        end
                    end
                    count{j,1}{i,1} = length(b);
                    index{j,1}{i,1} = b;
                    laps = b;
                    numlaps = length(b);
                    clear k dif id tmpo z
                    
                    
                    % Plot velocity and position
                    posx = store{j,1}{i,1}(:,1);
                    posy = store{j,1}{i,1}(:,2);
                    t = [0:dt:(length(posx)*dt)-dt];
                    
                    for f = 1:(length(posx)-1),
                        v{j,1}{i,1}(f,1) = real(sqrt((abs(posx(f+1)-posx(f)))*pixcorrection/(t(f+1)-t(f))^2 + (posy(f+1)-posy(f))*pixcorrection/(t(f+1)-t(f))^2));
                    end
                    v{j,1}{i,1} = conv(v{j,1}{i,1}, gaussFilter, 'same');
                    velocity = v{j,1}{i,1};
                    
                    figure(),
                    plot([dt : dt : length(v{j,1}{i,1})*dt],v{j,1}{i,1}(:,1))
                    hold on
                    c = (posx*pixcorrection)-min(posx*pixcorrection); %take the min as zero
                    plot([dt : dt : length(c)*dt],c), xlabel('Time (sec)'),ylabel('Distance (cm)'),title ('Velocity and Position in function of Time'),ylim([0 190]),yline(xi,'--'),yline(xf,'--')
                    
                    
                    posx = store{j,1}{i,1}(:,1)*pixcorrection;
                    posy = store{j,1}{i,1}(:,2)*pixcorrection;
                    posx_snout = tmp{1}{:,2}*pixcorrection;
                    posy_snout = tmp{1}{:,3}*pixcorrection;
                    %             save ([subFolders(i+2).folder,'\',subFolders(i+2).name,'\laps.mat'],'laps','numlaps','posx','posy','velocity');
                    %
                    %             save ([subFolders(i+2).folder,'\',subFolders(i+2).name,'\laps.mat'],'laps','numlaps','posx','posy','velocity');
                    if i == 2 || i == 4 % save only the 2nd and 4th sessions in root folder
                        save ([listing(j+2).folder,'\',listing(j+2).name,'\laps',num2str(session),'.mat'],'laps','numlaps','posx','posy','velocity','posx_snout','posy_snout');
                        session = session+1;
                    end
                    clear tmp0 tmpo posx posy t c laps numlaps velocity
                end
            end
%         end
    end
        clear session
        j
end
clear tt X xf xi Y sigma j i de mySmooth pixcorrection