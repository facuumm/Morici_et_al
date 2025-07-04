%% Instructions
% This script generates mat file with the timestamps of the start and end
% of the Camara, Shock, left and Rigth valves

clear
clc
close all

%% Parameters
% change currentdir deppending of the rat to be analyzed
currentdir = ('E:\Rat165\in_pyr\');
listing = dir(currentdir);


for j = 1 : length(listing)-2
    session = [listing(j+2).folder,'\',listing(j+2).name];
    cd(session)
    fileinfo = dir([cd,'\*digitalin.dat']);
%     fileinfo = dir([currentdir,'\',listing(j+2).name,'\*digitalin.dat']);
%     fileinfo = dir([currentdir,'\',listing(j+2).name,'\*.lfp']);
    file = [fileinfo.folder,'\',fileinfo.name];
    
    num_samples = fileinfo.bytes/2; % int32 = 4 bytes
    fid = fopen(file, 'r');
    digitalin = fread(fid, num_samples, 'uint16');
    fclose(fid);
    % ch = (bitand(digitalin, 2.^[0:1:15]) > 0); % ch has a value of 0-15 here
    ch = (bitand(digitalin, 2.^[0:1:3]) > 0); % I just upload the first 4 channels
    clear digitalin fid
    
    FS = 20000;
    dt = 1/FS;
    
    time = [0:dt:length(ch)*dt-dt];
    
    camara = time(ToIntervals(ch(:,1)));
    shock = time(ToIntervals(ch(:,2)));
    leftvalve = time(ToIntervals(ch(:,3)));
    rightvalve = time(ToIntervals(ch(:,4)));
    
    
    
%     figure (),
%     subplot(411),plot(time,ch(:,1)),title('video')
%     subplot(412),plot(time,ch(:,2)),title('shock')
%     subplot(413),plot(time,ch(:,3)),title('left-valve')
%     subplot(414),plot(time,ch(:,4)),title('right-valve')
%     
    %% save position of shocks in the intan
    save([cd,'\digitalin.mat'],'camara','shock','leftvalve','rightvalve')
    clear camara shock leftvalve rightvalve ch FS 
end
