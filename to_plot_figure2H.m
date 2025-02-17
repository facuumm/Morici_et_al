% coordinated
CoorAversive = reactivation.aversive.vHPC;
CoorReward = reactivation.reward.vHPC;

% uncoordinated dorsal
UnCoorDAversive = reactivation.aversive.vHPC;
UnCoorDReward = reactivation.reward.vHPC;

% uncoordinated ventral
UnCoorVAversive = reactivation.aversive.vHPC;
UnCoorVReward = reactivation.reward.vHPC;


%% Figure separated
figure,
subplot(131),
grps = [ones(size(CoorReward(:,1),1),1) ; ones(size(CoorAversive(:,1),1),1)*2];
Y = [CoorReward(:,1);CoorAversive(:,1)];
scatter(grps,Y,"filled",'jitter','on', 'jitterAmount',0.1),hold on
scatter([1 2] , [nanmean(CoorReward(:,1)) nanmean(CoorAversive(:,1))],'filled'),xlim([0 3]),ylim([-1 1])
[h p] = ttest(CoorReward(:,1),0,'tail','right')
[h p] = ttest(CoorAversive(:,1),0,'tail','right')
[h p] = ttest2(CoorAversive(:,1),CoorReward(:,1))

subplot(132),
grps = [ones(size(UnCoorDReward(:,1),1),1) ; ones(size(UnCoorDAversive(:,1),1),1)*2];
Y = [UnCoorDReward(:,1);UnCoorDAversive(:,1)];
scatter(grps,Y,"filled",'jitter','on', 'jitterAmount',0.1),hold on
scatter([1 2] , [nanmean(UnCoorDReward(:,1)) nanmean(UnCoorDAversive(:,1))],'filled'),xlim([0 3]),ylim([-1 1])
[h p] = ttest(UnCoorDReward(:,1),0,'tail','right')
[h p] = ttest(UnCoorDAversive(:,1),0,'tail','right')

subplot(133),
grps = [ones(size(UnCoorVReward(:,1),1),1) ; ones(size(UnCoorVAversive(:,1),1),1)*2];
Y = [UnCoorVReward(:,1);UnCoorVAversive(:,1)];
scatter(grps,Y,"filled",'jitter','on', 'jitterAmount',0.1),hold on
scatter([1 2] , [nanmean(UnCoorVReward(:,1)) nanmean(UnCoorVAversive(:,1))],'filled'),xlim([0 3]),ylim([-1 1])
[h p] = ttest(UnCoorVReward(:,1),0,'tail','right')
[h p] = ttest(UnCoorVAversive(:,1),0,'tail','right')


u = [UnCoorDAversive(:,1) , UnCoorVAversive(:,1)]; u(sum(isnan(u),2)>0,:) = [];
UnCoorAversive = nanmean(u');

u = [UnCoorDReward(:,1) , UnCoorVReward(:,1)]; u(sum(isnan(u),2)>0,:) = [];
UnCoorReward = nanmean(u');

%% Figure all
x = CoorReward(:,1);
x(isnan(x)) = [];
y = CoorAversive(:,1);
y(isnan(y)) = [];

x1 = [UnCoorDReward(:,1) ; UnCoorVReward(:,1)];
x1(isnan(x1)) = [];
y1 = [UnCoorDAversive(:,1) ; UnCoorVAversive(:,1)];
y1(isnan(y1)) = [];

figure,
subplot(121),
grps = [ones(size(CoorReward(:,1),1),1) ; ones(size(CoorAversive(:,1),1),1)*2];
Y = [CoorReward(:,1);CoorAversive(:,1)];
scatter(grps,Y,"filled",'jitter','on', 'jitterAmount',0.1),hold on
scatter([1 2] , [nanmean(CoorReward(:,1)) nanmean(CoorAversive(:,1))],'filled'),xlim([0 3]),ylim([-1 1])

subplot(122),
grps = [ones(size(x1,1),1) ; ones(size(y1,1),1)*2];
Y = [x1(:,1);y1(:,1)];
scatter(grps,Y,"filled",'jitter','on', 'jitterAmount',0.1),hold on
scatter([1 2] , [nanmean(x1(:,1)) nanmean(y1(:,1))],'filled'),xlim([0 3]),ylim([-1 1])

data = [x ; y ; x1(:,1) ; y1(:,1)];
ref = [ones(length(x),1)*2 ; ones(length(y),1)*2 ; ones(length(x1),1) ; ones(length(y1),1)];
ref = [ref , [ones(length(x),1) ; ones(length(y),1)*2 ; ones(length(x1),1) ; ones(length(y1),1)*2]];


perform2WayANOVANonPaired(data, ref(:,1), ref(:,2))

save([cd,'\Extended_data_ventral_Assemblies.mat'])