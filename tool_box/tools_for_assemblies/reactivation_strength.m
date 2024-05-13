function [R] = reactivation_strength(patterns , cond , SpikeTrain , Is , th , type , config , normalization, templates)
% This function calculate Reactivation Strength (van de Ven et al (2016)).
% Find the peaks of the zscored assemblies strength and calculate the mean
% during Pre sleep and Post sleep. Then, it calculates Post-Pre.
%
% [R] = reactivation_strength(patterns , cond , SpikeTrain , is , th , type)
%
% --- Inputs ---
% patterns: matrix with the weigths for each cell in each assembly.
%           Structure: Single-Units x Assemblies (rows x column)
%
% cond: logical, to select which pattern will be used (1, include / 0, not include)
%       Row vector with the same elements as number of assemblies.
%       Example,
%       patterns:     A1    A2    A3
%                 SU1 0.6   0.1   0.1       ---> cond: 0   1   0
%                 SU2 0.3   0.6   0.1       In this case, only the second
%                 SU3 0.1   0.3   0.3       assembly (A2) will be selected.
%                 SU4 0.1   0.1   0.6
%
% SpikeTrain: matrix, First column, time bins, Rest columns Spike Counts
%             Example:   t    SU1   SU2   SU3   SU4
%                        0     3     5     1     3    This Spike Train was
%                       0.1    1     2     1     0    constructed using a
%                       0.2    3     5     1     3    0.1 sec time bin.
%                       ...   ...   ...   ...   ...
%
% Is: structure, it contains the periods of time (in sec) of interest.
%     Struture:
%                 Is.baseline     Should be a vector with the same length
%                 Is.aversive     of the SpikeTrain with logical values to
%                 Is.reward       include or not the time bin.
%                 Is.runaversive
%                 Is.runreward
%
%     Example:
%             is.aversive:
%                         b1  b2  b3  b4  b5  b6  b7  b8  b9
%                         0   0   1   1   1   1   0   0   0
%
% th: float/int, threshold value for peaks detection. Note that this
%     function zscored the assemblies activity.
%
% type: String, 'A' if is Aversive assemblies 'R' if they are Reward
%
% config: int, if is 1 Aversive occurs first, if is 2, Reward was first
%
% normalization: int, if is 1 then (mean(Post) - mean(Pre))/mean(Post+Pre)
%                     if is 0 then (mean(Post) - mean(Pre))
%
% templates: 2-columns matrix, first column for dHPC and second for vHPC.
%            for more details see assembly_activity_only_joint
%
% --- OUTPUT ---
% R: column vector storing all the Reactivation Strength values for sleep
%    and awake periods.
%
%    S1 ---- *C1*  ---- S2 ----  ---- C2 ---- S3
%
%    * Assemblies detected during this session. It could do the same if C2
%    is the condition of interest. Just be sure you are choosing correctly
%    the type and config inputs.
%
%    Structure:  mean(S1)-mean(S2)   mean(C1)   mean(C2)
%                         R1           A11        A21
%                         R2           A12        A22
%                         R3           A13        A23
%                         ...          ...        ...
%
% requirments:
%       assembly_activity.m from Lopes-dos-Santos et al 2013 (*, see below)
%
% Morici Juan Facundo 02/2024

bins = SpikeTrain(:,1);
dt = bins(2)-bins(1); % delta time
spks = SpikeTrain(:,2:end);
iterations = 1;
prc = 75;
if not(isempty(templates))
    a = assembly_activity_only_joint(patterns(:,cond) , spks',templates(:,1),templates(:,2));
else
    a = assembly_activity(patterns(:,cond) , spks');
end

A = a;
a = zscore(a,1,2);
R = [];

for i = 1:size(a,1)

    MB = A(i,Is.baseline);% MB = MB(MB>=0);
    MR = A(i,Is.reward);% MR = MR(MR>=0);
    MA = A(i,Is.aversive);% MA = MA(MA>=0);
    
    % mean activation
    MB = nanmean(MB);
    MR = nanmean(MR);
    MA = nanmean(MA);

    % using the time vector for plotting   
    [pks.baseline,loc.baseline] = findpeaks(a(i,Is.baseline),bins(Is.baseline),'MinPeakHeight',th);
    [pks.reward,loc.reward] = findpeaks(a(i,Is.reward),bins(Is.reward),'MinPeakHeight',th);
    [pks.aversive,loc.aversive] = findpeaks(a(i,Is.aversive),bins(Is.aversive),'MinPeakHeight',th);
    
    [pks.runreward,loc.reward] = findpeaks(a(i,Is.runreward),bins(Is.runreward),'MinPeakHeight',th);
    [pks.runaversive,loc.aversive] = findpeaks(a(i,Is.runaversive),bins(Is.runaversive),'MinPeakHeight',th);
    
    [pks.all,loc.all] = findpeaks(a(i,:),bins,'MinPeakHeight',th);
    
    
    % shuffle peaks to create surrogated
    surrogated = [];
    for ii = 1 : iterations
        tmp = Shuffle2D([loc.all , pks.all']);
        
        pks.b = Restrict(tmp , ToIntervals(bins,Is.baseline));
        pks.b = pks.b(:,2);
        
        pks.a = Restrict(tmp , ToIntervals(bins,Is.aversive));
        pks.a = pks.a(:,2);
        
        pks.r = Restrict(tmp , ToIntervals(bins,Is.reward));
        pks.r = pks.r(:,2);
        
        
        if type == 'A' % check if is aversive assembly
            if config == 1
                if normalization
                    surrogated = [surrogated ; (nanmean(pks.a) - nanmean(pks.b))/(nanmean([pks.a ; pks.b]))];
                else
                    surrogated = [surrogated ; (nanmean(pks.a) - nanmean(pks.b))];
                end
            else
                if normalization
                    surrogated = [surrogated ; (nanmean(pks.a) - nanmean(pks.r))/(nanmean([pks.a ; pks.r]))];
                else
                    surrogated = [surrogated ; (nanmean(pks.a) - nanmean(pks.r))];
                end
            end
            
        elseif type == 'R' % check if is reward assembly
            if config == 2
                if normalization
                    surrogated = [surrogated ; (nanmean(pks.r) - nanmean(pks.b))/(nanmean([pks.r ; pks.b]))];
                else
                    surrogated = [surrogated ; (nanmean(pks.r) - nanmean(pks.b))];
                end
            else
                if normalization
                    surrogated = [surrogated ; (nanmean(pks.r) - nanmean(pks.a))/(nanmean([pks.r ; pks.a]))];
                else
                    surrogated = [surrogated ; (nanmean(pks.r) - nanmean(pks.a))];
                end
            end
        end
        clear tmp
    end
    
    s = nanstd(surrogated);
%     m = nanmean(surrogated);
    m = prctile(surrogated,prc);
    m = 0;
    clear s
    
    % calculation of time across conditions to normalize pks
    timeB = (sum(Is.baseline)*dt)/60;
    timeR = (sum(Is.reward)*dt)/60;
    timeA = (sum(Is.aversive)*dt)/60;
    
    %% Calculation of Reactivation
    if type == 'A' % check if is aversive assembly
        if config == 1
            if normalization
                strength = (nanmean(pks.aversive) - nanmean(pks.baseline))/(nanmean([pks.aversive , pks.baseline]));
                
%                 tmp1 = Restrict(Is.timestamps.sleep.aversive,[Is.timestamps.aversiveSleep(1) Is.timestamps.aversiveSleep(1)+1800]);
%                 tmp2 = Restrict(Is.timestamps.sleep.baseline,[Is.timestamps.baselineSleep(2)-1800 Is.timestamps.baselineSleep(2)]);
                tmp1 = Is.timestamps.sleep.aversive;
                tmp2 = Is.timestamps.sleep.baseline;

 
                p1 = Restrict(loc.all,tmp1);
                p2 = Restrict(loc.all,tmp2);
                
                strength1 = ((size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) - (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60))) / ((size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) + (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60)));
%                 strength1 = (length(pks.aversive)/timeA - length(pks.baseline)/timeB);

                FR.pre = (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));
                FR.post = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60));
                
                R = [R ; strength FR.pre FR.post nanmean(pks.runaversive) nanmean(pks.runreward) strength > m strength1 MA MB];
                clear strength FR strength1
            else
                strength = (nanmean(pks.aversive) - nanmean(pks.baseline));
%                 strength1 = (length(pks.aversive)/timeA - length(pks.baseline)/timeB)/(sum([length(pks.aversive) , length(pks.baseline)])/(timeA+timeB));
%                 strength1 = (length(pks.aversive)/timeA - length(pks.baseline)/timeB);

%                 tmp1 = Restrict(Is.timestamps.sleep.aversive,[Is.timestamps.aversiveSleep(1) Is.timestamps.aversiveSleep(1)+1800]);
%                 tmp2 = Restrict(Is.timestamps.sleep.baseline,[Is.timestamps.baselineSleep(2)-1800 Is.timestamps.baselineSleep(2)]);
                tmp1 = Is.timestamps.sleep.aversive;
                tmp2 = Is.timestamps.sleep.baseline;

                p1 = Restrict(loc.all,tmp1);
                p2 = Restrict(loc.all,tmp2);
                
                strength1 = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) - (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));

                FR.pre = (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));
                FR.post = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60));
                
                R = [R ; strength FR.pre FR.post nanmean(pks.runaversive) nanmean(pks.runreward) strength > m strength1 MA MB];
                clear strength FR strength1
            end
        else
            if normalization
                strength = (nanmean(pks.aversive) - nanmean(pks.reward))/(nanmean([pks.aversive , pks.reward]));
%                 strength1 = (length(pks.aversive)/timeA - length(pks.reward)/timeR)/(sum([length(pks.aversive) , length(pks.reward)])/(timeA+timeR));
%                 strength1 = (length(pks.aversive)/timeA - length(pks.reward)/timeR);
               
%                 tmp1 = Restrict(Is.timestamps.sleep.aversive,[Is.timestamps.aversiveSleep(1) Is.timestamps.aversiveSleep(1)+1800]);
%                 tmp2 = Restrict(Is.timestamps.sleep.reward,[Is.timestamps.rewardSleep(2)-1800 Is.timestamps.rewardSleep(2)]);
                tmp1 = Is.timestamps.sleep.aversive;
                tmp2 = Is.timestamps.sleep.reward;

                p1 = Restrict(loc.all,tmp1);
                p2 = Restrict(loc.all,tmp2);
                
                strength1 = ((size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) - (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60))) / ((size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) + (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60)));

                FR.pre = (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));
                FR.post = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60));
                
                R = [R ; strength FR.pre FR.post nanmean(pks.runaversive) nanmean(pks.runreward) strength > m strength1 MA MR];
                clear strength FR
            else
                strength = (nanmean(pks.aversive) - nanmean(pks.reward));
%                 strength1 = (length(pks.aversive)/timeA - length(pks.reward)/timeR)/(sum([length(pks.aversive) , length(pks.reward)])/(timeA+timeR));
%                 strength1 = (length(pks.aversive)/timeA - length(pks.reward)/timeR);

%                 tmp1 = Restrict(Is.timestamps.sleep.aversive,[Is.timestamps.aversiveSleep(1) Is.timestamps.aversiveSleep(1)+1800]);
%                 tmp2 = Restrict(Is.timestamps.sleep.reward,[Is.timestamps.rewardSleep(2)-1800 Is.timestamps.rewardSleep(2)]);
                tmp1 = Is.timestamps.sleep.aversive;
                tmp2 = Is.timestamps.sleep.reward;
                
                p1 = Restrict(loc.all,tmp1);
                p2 = Restrict(loc.all,tmp2);
                
                strength1 = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) - (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));


                FR.pre = (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));
                FR.post = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60));
                
                R = [R ; strength FR.pre FR.post nanmean(pks.runaversive) nanmean(pks.runreward) strength > m strength1 MA MR];

                clear strength FR strength1
            end
        end
        
    elseif type == 'R' % check if is reward assembly
        if config == 2
            if normalization
                strength = (nanmean(pks.reward) - nanmean(pks.baseline))/(nanmean([pks.reward , pks.baseline]));
%                 strength1 = (length(pks.reward)/timeR - length(pks.baseline)/timeB)/(sum([length(pks.reward) , length(pks.baseline)])/(timeR+timeB));
%                 strength1 = (length(pks.reward)/timeR - length(pks.baseline)/timeB);

%                 tmp1 = Restrict(Is.timestamps.sleep.reward,[Is.timestamps.rewardSleep(1) Is.timestamps.rewardSleep(1)+1800]);
%                 tmp2 = Restrict(Is.timestamps.sleep.baseline,[Is.timestamps.baselineSleep(2)-1800 Is.timestamps.baselineSleep(2)]);
                tmp1 = Is.timestamps.sleep.reward;
                tmp2 = Is.timestamps.sleep.baseline;
                
                p1 = Restrict(loc.all,tmp1);
                p2 = Restrict(loc.all,tmp2);
                
                strength1 = ((size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) - (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60))) / ((size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) + (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60)));

                FR.pre = (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));
                FR.post = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60));
                
                R = [R ; strength FR.pre FR.post nanmean(pks.runreward) nanmean(pks.runaversive) strength > m strength1 MR MB];
                clear strength FR strength1
            else
                strength = (nanmean(pks.reward) - nanmean(pks.baseline));
%                 strength1 = (length(pks.reward)/timeR - length(pks.baseline)/timeB)/(sum([length(pks.reward) , length(pks.baseline)])/(timeR+timeB));
%                 strength1 = (length(pks.reward)/timeR - length(pks.baseline)/timeB);

%                 tmp1 = Restrict(Is.timestamps.sleep.reward,[Is.timestamps.rewardSleep(1) Is.timestamps.rewardSleep(1)+1800]);
%                 tmp2 = Restrict(Is.timestamps.sleep.baseline,[Is.timestamps.baselineSleep(2)-1800 Is.timestamps.baselineSleep(2)]);
                tmp1 = Is.timestamps.sleep.reward;
                tmp2 = Is.timestamps.sleep.baseline;
                
                p1 = Restrict(loc.all,tmp1);
                p2 = Restrict(loc.all,tmp2);
                
                strength1 = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) - (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));

                FR.pre = (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));
                FR.post = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60));                
                
                R = [R ; strength FR.pre FR.post nanmean(pks.runreward) nanmean(pks.runaversive) strength > m strength1 MR MB];
                clear strength FR strength1
            end
        else
            if normalization
                strength = (nanmean(pks.reward) - nanmean(pks.aversive))/(nanmean([pks.reward , pks.aversive]));
%                 strength1 = (length(pks.reward)/timeR - length(pks.aversive)/timeA)/(sum([length(pks.reward) , length(pks.aversive)])/(timeR+timeA));
%                 strength1 = (length(pks.reward)/timeR - length(pks.aversive)/timeA);

%                 tmp1 = Restrict(Is.timestamps.sleep.reward,[Is.timestamps.rewardSleep(1) Is.timestamps.rewardSleep(1)+1800]);
%                 tmp2 = Restrict(Is.timestamps.sleep.aversive,[Is.timestamps.aversiveSleep(2)-1800 Is.timestamps.aversiveSleep(2)]);
                tmp1 = Is.timestamps.sleep.reward;
                tmp2 = Is.timestamps.sleep.aversive;

                p1 = Restrict(loc.all,tmp1);
                p2 = Restrict(loc.all,tmp2);
                
                strength1 = ((size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) - (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60))) / ((size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) + (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60)));

                FR.pre = (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));
                FR.post = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60));
                
                R = [R ; strength FR.pre FR.post  nanmean(pks.runreward) nanmean(pks.runaversive) strength > m strength1 MR MA];
                clear strength FR strength1
            else
                strength = (nanmean(pks.reward) - nanmean(pks.aversive));
%                 strength1 = (length(pks.reward)/timeR - length(pks.aversive)/timeA)/(sum([length(pks.reward) , length(pks.aversive)])/(timeR+timeA));
%                 strength1 = (length(pks.reward)/timeR - length(pks.aversive)/timeA);

%                 tmp1 = Restrict(Is.timestamps.sleep.reward,[Is.timestamps.rewardSleep(1) Is.timestamps.rewardSleep(1)+1800]);
%                 tmp2 = Restrict(Is.timestamps.sleep.aversive,[Is.timestamps.aversiveSleep(2)-1800 Is.timestamps.aversiveSleep(2)]);
                tmp1 = Is.timestamps.sleep.reward;
                tmp2 = Is.timestamps.sleep.aversive;
                
                p1 = Restrict(loc.all,tmp1);
                p2 = Restrict(loc.all,tmp2);
                
                strength1 = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) - (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));

                FR.pre = (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));
                FR.post = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60));
                
                R = [R ; strength FR.pre FR.post  nanmean(pks.runreward) nanmean(pks.runaversive) strength > m strength1 MR MA];
                clear strength FR strength1
            end
        end
    end
    clear m pks loc MB MR MA
end

end