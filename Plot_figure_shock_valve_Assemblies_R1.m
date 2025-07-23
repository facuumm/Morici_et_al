% %% Decorrelation of speed and response
% time = [-5:0.1:5];
% [h hh] = min(abs(time-(1)));
% Decorrelated.joint.aversive = [];
% for i = 1 : size(Response.joint.aversive,2)
%     id = modulated.joint.aversive(i,:);
%     template = Response.speed.id;
%     template = and(id(2)==template(:,1) , id(2)==template(:,2));
%     [deco] = decorr(Response.joint.aversive(:,i),Response.speed.data(:,template));
% %     n = nanmean(deco(1:hh));
%     Decorrelated.joint.aversive = [Decorrelated.joint.aversive , deco]; clear deco n
% end
% 
% Decorrelated.joint.reward = [];
% for i = 1 : size(Response.joint.reward,2)
%     id = modulated.dHPC.reward(i,:);
%     template = Response.speed.id;
%     template = and(id(2)==template(:,1) , id(2)==template(:,2));
%     [deco] = decorr(Response.joint.reward(:,i),Response.speed.data(:,template));
% %     n = nanmean(deco(1:hh));
%     Decorrelated.joint.reward = [Decorrelated.joint.reward , deco]; clear deco n
% end
% 
% Decorrelated.dHPC.aversive = [];
% for i = 1 : size(Response.dHPC.aversive,2)
%     id = modulated.dHPC.aversive(i,:);
%     template = Response.speed.id;
%     template = and(id(2)==template(:,1) , id(2)==template(:,2));
%     [deco] = decorr(Response.dHPC.aversive(:,i),Response.speed.data(:,template));
% %     n = nanmean(deco(1:hh));
%     Decorrelated.dHPC.aversive = [Decorrelated.dHPC.aversive , deco]; clear deco n
% end
% 
% Decorrelated.dHPC.reward = [];
% for i = 1 : size(Response.dHPC.reward,2)
%     id = modulated.dHPC.reward(i,:);
%     template = Response.speed.id;
%     template = and(id(2)==template(:,1) , id(2)==template(:,2));
%     [deco] = decorr(Response.dHPC.reward(:,i),Response.speed.data(:,template));
% %     n = nanmean(deco(1:hh));
%     Decorrelated.dHPC.reward = [Decorrelated.dHPC.reward , deco]; clear deco n
% end
% 
% Decorrelated.vHPC.aversive = [];
% for i = 1 : size(Response.vHPC.aversive,2)
%     id = modulated.vHPC.aversive(i,:);
%     template = Response.speed.id;
%     template = and(id(2)==template(:,1) , id(2)==template(:,2));
%     [deco] = decorr(Response.vHPC.aversive(:,i),Response.speed.data(:,template));
% %     n = nanmean(deco(1:hh));
%     Decorrelated.vHPC.aversive = [Decorrelated.vHPC.aversive , deco]; clear deco n
% end
% 
% Decorrelated.vHPC.reward = [];
% for i = 1 : size(Response.vHPC.reward,2)
%     id = modulated.vHPC.reward(i,:);
%     template = Response.speed.id;
%     template = and(id(2)==template(:,1) , id(2)==template(:,2));
%     [deco] = decorr(Response.dHPC.reward(:,i),Response.speed.data(:,template));
% %     n = nanmean(deco(1:hh));
%     Decorrelated.vHPC.reward = [Decorrelated.vHPC.reward , deco]; clear deco n
% end
% 
%% Plot Aversive Assemblies
figure,
subplot(131)
i = modulated.joint.aversive(:,1) == 1;
m = nanmean(Response.joint.aversive(:,i)');
s = nansem(Response.joint.aversive(:,i)');
plot([-5:0.1:5],m,'k'), hold on
ciplot(m-s , m+s , [-5:0.1:5] , 'k'),alpha 0.10
ylim([-1 3])
xline(0), hold on, xline(1)

subplot(132)
i = modulated.joint.aversive(:,1) == -1;
m = nanmean(Response.joint.aversive(:,i)');
s = nansem(Response.joint.aversive(:,i)');
plot([-5:0.1:5],m,'k'), hold on
ciplot(m-s , m+s , [-5:0.1:5] , 'k'),alpha 0.10
ylim([-1 3])
xline(0), hold on, xline(1)

subplot(133)
i = modulated.joint.aversive(:,1) == 0;
m = nanmean(Response.joint.aversive(:,i)');
s = nansem(Response.joint.aversive(:,i)');
plot([-5:0.1:5],m,'k'), hold on
ciplot(m-s , m+s , [-5:0.1:5] , 'k'),alpha 0.10
ylim([-1 3])
xline(0), hold on, xline(1)

figure,
p1 = (sum(modulated.joint.aversive(:,1) == 1)/length(modulated.joint.aversive(:,1)))*100;
p2 = (sum(modulated.joint.aversive(:,1) == -1)/length(modulated.joint.aversive(:,1)))*100;
p3 = 100 - p1 - p2;
pie([p1 p2 p3] , {'Increased' , 'Decreased' , 'None'})

%% dHPC
figure,
subplot(131)
i = modulated.dHPC.aversive(:,1) == 1;
m = nanmean(Response.dHPC.aversive(:,i)');
s = nansem(Response.dHPC.aversive(:,i)');
plot([-5:0.1:5],m,'k'), hold on
ciplot(m-s , m+s , [-5:0.1:5] , 'k'),alpha 0.10
ylim([-1 3])
xline(0), hold on, xline(1)

subplot(132)
i = modulated.dHPC.aversive(:,1) == -1;
m = nanmean(Response.dHPC.aversive(:,i)');
s = nansem(Response.dHPC.aversive(:,i)');
plot([-5:0.1:5],m,'k'), hold on
ciplot(m-s , m+s , [-5:0.1:5] , 'k'),alpha 0.10
ylim([-1 3])
xline(0), hold on, xline(1)

subplot(133)
i = modulated.dHPC.aversive(:,1) == 0;
m = nanmean(Response.dHPC.aversive(:,i)');
s = nansem(Response.dHPC.aversive(:,i)');
plot([-5:0.1:5],m,'k'), hold on
ciplot(m-s , m+s , [-5:0.1:5] , 'k'),alpha 0.10
ylim([-1 3])
xline(0), hold on, xline(1)

figure,
p1 = (sum(modulated.dHPC.aversive(:,1) == 1)/length(modulated.dHPC.aversive(:,1)))*100;
p2 = (sum(modulated.dHPC.aversive(:,1) == -1)/length(modulated.dHPC.aversive(:,1)))*100;
p3 = 100 - p1 - p2;
pie([p1 p2 p3] , {'Increased' , 'Decreased' , 'None'})

%% vHPC
figure,
subplot(131)
i = modulated.vHPC.aversive(:,1) == 1;
m = nanmean(Response.vHPC.aversive(:,i)');
s = nansem(Response.vHPC.aversive(:,i)');
plot([-5:0.1:5],m,'k'), hold on
ciplot(m-s , m+s , [-5:0.1:5] , 'k'),alpha 0.10
ylim([-1 3])
xline(0), hold on, xline(1)

subplot(132)
i = modulated.vHPC.aversive(:,1) == -1;
m = nanmean(Response.vHPC.aversive(:,i)');
s = nansem(Response.vHPC.aversive(:,i)');
plot([-5:0.1:5],m,'k'), hold on
ciplot(m-s , m+s , [-5:0.1:5] , 'k'),alpha 0.10
ylim([-1 3])
xline(0), hold on, xline(1)

subplot(133)
i = modulated.vHPC.aversive(:,1) == 0;
m = nanmean(Response.vHPC.aversive(:,i)');
s = nansem(Response.vHPC.aversive(:,i)');
plot([-5:0.1:5],m,'k'), hold on
ciplot(m-s , m+s , [-5:0.1:5] , 'k'),alpha 0.10
ylim([-1 3])
xline(0), hold on, xline(1)

figure,
p1 = (sum(modulated.vHPC.aversive(:,1) == 1)/length(modulated.vHPC.aversive(:,1)))*100;
p2 = (sum(modulated.vHPC.aversive(:,1) == -1)/length(modulated.vHPC.aversive(:,1)))*100;
p3 = 100 - p1 - p2;
pie([p1 p2 p3] , {'Increased' , 'Decreased' , 'None'})


%% Plot Reward Assemblies
figure,
subplot(131)
i = modulated.joint.reward(:,1) == 1;
m = nanmean(Response.joint.reward(:,i)');
s = nansem(Response.joint.reward(:,i)');
plot([-5:0.1:5],m,'k'), hold on
ciplot(m-s , m+s , [-5:0.1:5] , 'k'),alpha 0.10
ylim([-1.5 3])
xline(0), hold on, xline(1)

subplot(132)
i = modulated.joint.reward(:,1) == -1;
m = nanmean(Response.joint.reward(:,i)');
s = nansem(Response.joint.reward(:,i)');
plot([-5:0.1:5],m,'k'), hold on
ciplot(m-s , m+s , [-5:0.1:5] , 'k'),alpha 0.10
ylim([-1.5 3])
xline(0), hold on, xline(1)

subplot(133)
i = modulated.joint.reward(:,1) == 0;
m = nanmean(Response.joint.reward(:,i)');
s = nansem(Response.joint.reward(:,i)');
plot([-5:0.1:5],m,'k'), hold on
ciplot(m-s , m+s , [-5:0.1:5] , 'k'),alpha 0.10
ylim([-1.5 3])
xline(0), hold on, xline(1)

figure,
p1 = (sum(modulated.joint.reward(:,1) == 1)/length(modulated.joint.reward(:,1)))*100;
p2 = (sum(modulated.joint.reward(:,1) == -1)/length(modulated.joint.reward(:,1)))*100;
p3 = 100 - p1 - p2;
pie([p1 p2 p3] , {'Increased' , 'Decreased' , 'None'})

%% dHPC
figure,
subplot(131)
i = modulated.dHPC.reward(:,1) == 1;
m = nanmean(Response.dHPC.reward(:,i)');
s = nansem(Response.dHPC.reward(:,i)');
plot([-5:0.1:5],m,'k'), hold on
ciplot(m-s , m+s , [-5:0.1:5] , 'k'),alpha 0.10
ylim([-1.5 3])
xline(0), hold on, xline(1)

subplot(132)
i = modulated.dHPC.reward(:,1) == -1;
m = nanmean(Response.dHPC.reward(:,i)');
s = nansem(Response.dHPC.reward(:,i)');
plot([-5:0.1:5],m,'k'), hold on
ciplot(m-s , m+s , [-5:0.1:5] , 'k'),alpha 0.10
ylim([-1.5 3])
xline(0), hold on, xline(1)

subplot(133)
i = modulated.dHPC.reward(:,1) == 0;
m = nanmean(Response.dHPC.reward(:,i)');
s = nansem(Response.dHPC.reward(:,i)');
plot([-5:0.1:5],m,'k'), hold on
ciplot(m-s , m+s , [-5:0.1:5] , 'k'),alpha 0.10
ylim([-1.5 3])
xline(0), hold on, xline(1)

figure,
p1 = (sum(modulated.dHPC.reward(:,1) == 1)/length(modulated.dHPC.reward(:,1)))*100;
p2 = (sum(modulated.dHPC.reward(:,1) == -1)/length(modulated.dHPC.reward(:,1)))*100;
p3 = 100 - p1 - p2;
pie([p1 p2 p3] , {'Increased' , 'Decreased' , 'None'})

%% vHPC
figure,
subplot(131)
i = modulated.vHPC.reward(:,1) == 1;
m = nanmean(Response.vHPC.reward(:,i)');
s = nansem(Response.vHPC.reward(:,i)');
plot([-5:0.1:5],m,'k'), hold on
ciplot(m-s , m+s , [-5:0.1:5] , 'k'),alpha 0.10
ylim([-1.5 3])
xline(0), hold on, xline(1)

subplot(132)
i = modulated.vHPC.reward(:,1) == -1;
m = nanmean(Response.vHPC.reward(:,i)');
s = nansem(Response.vHPC.reward(:,i)');
plot([-5:0.1:5],m,'k'), hold on
ciplot(m-s , m+s , [-5:0.1:5] , 'k'),alpha 0.10
ylim([-1.5 3])
xline(0), hold on, xline(1)

subplot(133)
i = modulated.vHPC.reward(:,1) == 0;
m = nanmean(Response.vHPC.reward(:,i)');
s = nansem(Response.vHPC.reward(:,i)');
plot([-5:0.1:5],m,'k'), hold on
ciplot(m-s , m+s , [-5:0.1:5] , 'k'),alpha 0.10
ylim([-1.5 3])
xline(0), hold on, xline(1)

figure,
p1 = (sum(modulated.vHPC.reward(:,1) == 1)/length(modulated.vHPC.reward(:,1)))*100;
p2 = (sum(modulated.vHPC.reward(:,1) == -1)/length(modulated.vHPC.reward(:,1)))*100;
p3 = 100 - p1 - p2;
pie([p1 p2 p3] , {'Increased' , 'Decreased' , 'None'})

%% Plot and Analyse
figure,
grps = [ones(size(values.dHPC.shock,1),1) ; ones(size(values.dHPC.decorrelated,1),1)*2 ; ones(size(values.vHPC.shock,1),1)*3 ; ones(size(values.vHPC.decorrelated,1),1)*4];
y = [values.dHPC.shock ; values.dHPC.decorrelated ; values.vHPC.shock ; values.vHPC.decorrelated]; 
scatter(grps,y,'filled','jitter','on', 'jitterAmount',0.1),ylim([-0.5 6]) , xlim([0 5]),hold on
scatter([1 2 3 4] , [nanmean(values.dHPC.shock) ; nanmean(values.dHPC.decorrelated) ; nanmean(values.vHPC.shock) ; nanmean(values.vHPC.decorrelated)] , 'filled')

% Run analysis
y = [values.dHPC.shock ; values.dHPC.decorrelated ; values.vHPC.shock ; values.vHPC.decorrelated]; 

structure = [ones(size(values.dHPC.shock,1),1) ; ones(size(values.dHPC.decorrelated,1),1) ; ...
             2*ones(size(values.vHPC.shock,1),1) ; 2*ones(size(values.vHPC.decorrelated,1),1)];

condition = [ones(size(values.dHPC.shock,1),1) ; 2*ones(size(values.dHPC.decorrelated,1),1) ; ...
             ones(size(values.vHPC.shock,1),1) ; 2*ones(size(values.vHPC.decorrelated,1),1)];

perform2WayUnpairedANOVA(y, structure, condition)
