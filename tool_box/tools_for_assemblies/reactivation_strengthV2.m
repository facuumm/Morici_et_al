function [R P] = reactivation_strengthV2(patterns , cond , SpikeTrain , Is , th , type , config , normalization, templates)
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

tempo = 3600; % to restrictin time the amount of sleep that I want to include
w = 10;
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
% A = a;
R = [];
P = [];
for i = 1:size(a,1)
    % using the time vector for plotting
    [pks.runreward,loc.reward] = findpeaks(a(i,Is.runreward),bins(Is.runreward),'MinPeakHeight',th);
    [pks.runaversive,loc.aversive] = findpeaks(a(i,Is.runaversive),bins(Is.runaversive),'MinPeakHeight',th);
    
    [pks.all,loc.all] = findpeaks(a(i,:),bins,'MinPeakHeight',th);
    
    %% Calculation of Reactivation
    if type == 'A' % check if is aversive assembly
        if config == 1
            tmp1 = Restrict(Is.timestamps.sleep.aversive,[Is.timestamps.aversiveRun(1,2) (Is.timestamps.aversiveRun(1,2))+tempo]);
            tmp2 = Restrict(Is.timestamps.sleep.baseline,[(Is.timestamps.aversiveRun(1,1))-tempo Is.timestamps.aversiveRun(1,1)]);
            
            % concatenation of peaks
            temporal = [];
            iterator = 0;
            for ii = 1 : size(tmp1,1)
                t = Restrict(loc.all,tmp1(ii,:));
                t = (t - tmp1(ii,1));
                temporal = [temporal ; t + iterator];
                iterator = (tmp1(ii,2)-tmp1(ii,1))+iterator; % concatenation of time segments
                clear t
            end
%             h = histcounts(temporal,'BinEdges',[0 : w : tempo]);
            h = histcounts(temporal,10);
            P = [P , h']; clear h
            
            if normalization
                M1 = Restrict([loc.all pks.all'],tmp1);
                M2 = Restrict([loc.all pks.all'],tmp2);
                strength = (nanmean(M1(:,2)) - nanmean(M2(:,2)))/(nanmean(M1(:,2)) + nanmean(M2(:,2)));
                
                p1 = Restrict(loc.all,tmp1);
                p2 = Restrict(loc.all,tmp2);
                
                strength1 = ((size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) - (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60))) / ((size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) + (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60)));
                
                FR.pre = (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));
                FR.post = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60));
                
                t = Restrict([bins,A(i,:)'],Is.timestamps.sleep.aversive);
                t(:,1) = [dt:dt:size(t,1)*dt];
                edges = (0:w:tempo);
                [~,~,l]=histcounts(t(:,1),edges);
                t(l==0,:) = [];
                l(l==0) = [];
                meany = accumarray(l(:),t(:,2))./accumarray(l(:),1);
                
%                 c = (tempo/w) - size(meany,1);
%                 if c == 0
%                      R = [R , meany];
%                 else
%                     R = [R , [meany;nan(c,size(meany,2))]];
%                 end
%                 clear t edges l meany
                
                R = [R ; strength FR.pre FR.post nanmean(pks.runaversive) nanmean(pks.runreward) strength1];
                clear strength FR strength1 p1 p2 tmp1 tmp2
            else
                %                 strength = (nanmean(pks.aversive) - nanmean(pks.baseline));
                %                 tmp1 = Is.timestamps.sleep.aversive;
                %                 tmp2 = Is.timestamps.sleep.baseline;
                
                M1 = Restrict([loc.all pks.all'],tmp1);
                M2 = Restrict([loc.all pks.all'],tmp2);
                strength = (nanmean(M1(:,2)) - nanmean(M2(:,2)));
                
                p1 = Restrict(loc.all,tmp1);
                p2 = Restrict(loc.all,tmp2);
                
                strength1 = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) - (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));
                
                FR.pre = (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));
                FR.post = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60));
                
                t = Restrict([bins,A(i,:)'],Is.timestamps.sleep.aversive);
                t(:,1) = [dt:dt:size(t,1)*dt];
                edges = (0:w:tempo);
                [~,~,l]=histcounts(t(:,1),edges);
                t(l==0,:) = [];
                l(l==0) = [];
                meany = accumarray(l(:),t(:,2))./accumarray(l(:),1);
                
%                 c = (tempo/w) - size(meany,1);
%                 if c == 0
%                      R = [R , meany];
%                 else
%                     R = [R , [meany;nan(c,size(meany,2))]];
%                 end
%                 clear t edges l meany                
                R = [R ; strength FR.pre FR.post nanmean(pks.runaversive) nanmean(pks.runreward) strength1];
                clear strength FR strength1 p1 p2 tmp1 tmp2
            end
        else
            tmp1 = Restrict(Is.timestamps.sleep.aversive,[Is.timestamps.aversiveRun(1,2) (Is.timestamps.aversiveRun(1,2))+tempo]);
            tmp2 = Restrict(Is.timestamps.sleep.reward,[(Is.timestamps.aversiveRun(1,1))-tempo Is.timestamps.aversiveRun(1,1)]);
            
            % concatenation of peaks
            temporal = [];
            iterator = 0;
            for ii = 1 : size(tmp1,1)
                t = Restrict(loc.all,tmp1(ii,:));
                t = (t - tmp1(ii,1));
                temporal = [temporal ; t + iterator];
                iterator = (tmp1(ii,2)-tmp1(ii,1))+iterator; % concatenation of time segments
                clear t
            end
%             h = histcounts(temporal,'BinEdges',[0 : w : tempo]);
            h = histcounts(temporal,10);
            P = [P , h']; clear h
            
            if normalization
                M1 = Restrict([loc.all pks.all'],tmp1);
                M2 = Restrict([loc.all pks.all'],tmp2);
                strength = (nanmean(M1(:,2)) - nanmean(M2(:,2)))/(nanmean(M1(:,2)) + nanmean(M2(:,2)));
                
                p1 = Restrict(loc.all,tmp1);
                p2 = Restrict(loc.all,tmp2);
                
                strength1 = ((size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) - (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60))) / ((size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) + (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60)));
                
                FR.pre = (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));
                FR.post = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60));
                
                t = Restrict([bins,A(i,:)'],Is.timestamps.sleep.aversive);
                t(:,1) = [dt:dt:size(t,1)*dt];
                edges = (0:w:tempo);
                [~,~,l]=histcounts(t(:,1),edges);
                t(l==0,:) = [];
                l(l==0) = [];
                meany = accumarray(l(:),t(:,2))./accumarray(l(:),1);
                
%                 c = (tempo/w) - size(meany,1);
%                 if c == 0
%                      R = [R , meany];
%                 else
%                     R = [R , [meany;nan(c,size(meany,2))]];
%                 end
%                 clear t edges l meany                
                R = [R ; strength FR.pre FR.post nanmean(pks.runaversive) nanmean(pks.runreward) strength1];
                clear strength FR strength1 p1 p2 tmp1 tmp2
            else
                M1 = Restrict([loc.all pks.all'],tmp1);
                M2 = Restrict([loc.all pks.all'],tmp2);
                strength = (nanmean(M1(:,2)) - nanmean(M2(:,2)));
                
                p1 = Restrict(loc.all,tmp1);
                p2 = Restrict(loc.all,tmp2);
                
                strength1 = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) - (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));
                
                
                FR.pre = (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));
                FR.post = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60));
                
                t = Restrict([bins,A(i,:)'],Is.timestamps.sleep.aversive);
                t(:,1) = [dt:dt:size(t,1)*dt];
                edges = (0:w:tempo);
                [~,~,l]=histcounts(t(:,1),edges);
                t(l==0,:) = [];
                l(l==0) = [];
                meany = accumarray(l(:),t(:,2))./accumarray(l(:),1);
                
%                 c = (tempo/w) - size(meany,1);
%                 if c == 0
%                      R = [R , meany];
%                 else
%                     R = [R , [meany;nan(c,size(meany,2))]];
%                 end
%                 clear t edges l meany
                R = [R ; strength FR.pre FR.post nanmean(pks.runaversive) nanmean(pks.runreward) strength1];
                clear strength FR strength1 p1 p2 tmp1 tmp2
            end
        end
        
    elseif type == 'R' % check if is reward assembly
        if config == 2
            tmp1 = Restrict(Is.timestamps.sleep.reward,[Is.timestamps.rewardRun(1,2) (Is.timestamps.rewardRun(1,2))+tempo]);
            tmp2 = Restrict(Is.timestamps.sleep.baseline,[(Is.timestamps.rewardRun(1,1))-tempo Is.timestamps.rewardRun(1,1)]);
            % concatenation of peaks
            temporal = [];
            iterator = 0;
            for ii = 1 : size(tmp1,1)
                t = Restrict(loc.all,tmp1(ii,:));
                t = (t - tmp1(ii,1));
                temporal = [temporal ; t + iterator];
                iterator = (tmp1(ii,2)-tmp1(ii,1))+iterator; % concatenation of time segments
                clear t
            end
%             h = histcounts(temporal,'BinEdges',[0 : w : tempo]);
            h = histcounts(temporal,10);
            P = [P , h']; clear h
            
            if normalization
                M1 = Restrict([loc.all pks.all'],tmp1);
                M2 = Restrict([loc.all pks.all'],tmp2);
                strength = (nanmean(M1(:,2)) - nanmean(M2(:,2)))/(nanmean(M1(:,2)) + nanmean(M2(:,2)));
                
                p1 = Restrict(loc.all,tmp1);
                p2 = Restrict(loc.all,tmp2);
                
                strength1 = ((size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) - (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60))) / ((size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) + (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60)));
                
                FR.pre = (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));
                FR.post = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60));
                
                t = Restrict([bins,A(i,:)'],Is.timestamps.sleep.reward);
                t(:,1) = [dt:dt:size(t,1)*dt];
                edges = (0:w:tempo);
                [~,~,l]=histcounts(t(:,1),edges);
                t(l==0,:) = [];
                l(l==0) = [];
                meany = accumarray(l(:),t(:,2))./accumarray(l(:),1);
                
%                 c = (tempo/w) - size(meany,1);
%                 if c == 0
%                      R = [R , meany];
%                 else
%                     R = [R , [meany;nan(c,size(meany,2))]];
%                 end
%                 clear t edges l meany                
                R = [R ; strength FR.pre FR.post nanmean(pks.runreward) nanmean(pks.runaversive) strength1];
                clear strength FR strength1 p1 p2 tmp1 tmp2
            else
                M1 = Restrict([loc.all pks.all'],tmp1);
                M2 = Restrict([loc.all pks.all'],tmp2);
                strength = (nanmean(M1(:,2)) - nanmean(M2(:,2)));
                
                p1 = Restrict(loc.all,tmp1);
                p2 = Restrict(loc.all,tmp2);
                
                strength1 = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) - (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));
                
                FR.pre = (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));
                FR.post = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60));
                
                t = Restrict([bins,A(i,:)'],Is.timestamps.sleep.reward);
                t(:,1) = [dt:dt:size(t,1)*dt];
                edges = (0:w:tempo);
                [~,~,l]=histcounts(t(:,1),edges);
                t(l==0,:) = [];
                l(l==0) = [];
                meany = accumarray(l(:),t(:,2))./accumarray(l(:),1);
                
%                 c = (tempo/w) - size(meany,1);
%                 if c == 0
%                      R = [R , meany];
%                 else
%                     R = [R , [meany;nan(c,size(meany,2))]];
%                 end
%                 clear t edges l meany
                R = [R ; strength FR.pre FR.post nanmean(pks.runreward) nanmean(pks.runaversive) strength1];
                clear strength FR strength1 p1 p2 tmp1 tmp2
            end
        else
            tmp1 = Restrict(Is.timestamps.sleep.reward,[Is.timestamps.rewardRun(1,2) (Is.timestamps.rewardRun(1,2))+tempo]);
            tmp2 = Restrict(Is.timestamps.sleep.aversive,[(Is.timestamps.rewardRun(1,1))-tempo Is.timestamps.rewardRun(1,1)]);
            % concatenation of peaks
            temporal = [];
            iterator = 0;
            for ii = 1 : size(tmp1,1)
                t = Restrict(loc.all,tmp1(ii,:));
                t = (t - tmp1(ii,1));
                temporal = [temporal ; t + iterator];
                iterator = (tmp1(ii,2)-tmp1(ii,1))+iterator; % concatenation of time segments
                clear t
            end
%             h = histcounts(temporal,'BinEdges',[0 : w : tempo]);
            h = histcounts(temporal,10);
            P = [P , h']; clear h
            
            if normalization
                M1 = Restrict([loc.all pks.all'],tmp1);
                M2 = Restrict([loc.all pks.all'],tmp2);
                strength = (nanmean(M1(:,2)) - nanmean(M2(:,2)))/(nanmean(M1(:,2)) + nanmean(M2(:,2)));
                
                p1 = Restrict(loc.all,tmp1);
                p2 = Restrict(loc.all,tmp2);
                
                strength1 = ((size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) - (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60))) / ((size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) + (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60)));
                
                FR.pre = (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));
                FR.post = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60));
                
                t = Restrict([bins,A(i,:)'],Is.timestamps.sleep.reward);
                t(:,1) = [dt:dt:size(t,1)*dt];
                edges = (0:w:tempo);
                [~,~,l]=histcounts(t(:,1),edges);
                t(l==0,:) = [];
                l(l==0) = [];
                meany = accumarray(l(:),t(:,2))./accumarray(l(:),1);
                
%                 c = (tempo/w) - size(meany,1);
%                 if c == 0
%                      R = [R , meany];
%                 else
%                     R = [R , [meany;nan(c,size(meany,2))]];
%                 end
%                 clear t edges l meany
                R = [R ; strength FR.pre FR.post  nanmean(pks.runreward) nanmean(pks.runaversive) strength1];
                clear strength FR strength1 p1 p2 tmp1 tmp2
            else
                M1 = Restrict([loc.all pks.all'],tmp1);
                M2 = Restrict([loc.all pks.all'],tmp2);
                strength = (nanmean(M1(:,2)) - nanmean(M2(:,2)));
                
                p1 = Restrict(loc.all,tmp1);
                p2 = Restrict(loc.all,tmp2);
                
                strength1 = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60)) - (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));
                
                FR.pre = (size(p2,1)/(sum(tmp2(:,2)-tmp2(:,1))/60));
                FR.post = (size(p1,1)/(sum(tmp1(:,2)-tmp1(:,1))/60));
                
                t = Restrict([bins,A(i,:)'],Is.timestamps.sleep.reward);
                t(:,1) = [dt:dt:size(t,1)*dt];
                edges = (0:w:tempo);
                [~,~,l]=histcounts(t(:,1),edges);
                t(l==0,:) = [];
                l(l==0) = [];
                meany = accumarray(l(:),t(:,2))./accumarray(l(:),1);
                
%                 c = (tempo/w) - size(meany,1);
%                 if c == 0
%                      R = [R , meany];
%                 else
%                     R = [R , [meany;nan(c,size(meany,2))]];
%                 end
%                 clear t edges l meany
                
                R = [R ; strength FR.pre FR.post  nanmean(pks.runreward) nanmean(pks.runaversive) strength1];
                clear strength FR strength1 p1 p2 tmp1 tmp2
            end
        end
    end
    clear m pks loc MB MR MA
end

end