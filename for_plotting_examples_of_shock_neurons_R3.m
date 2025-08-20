% For plotting and checking if shuffle works
% dHPC
T = [-5 : 0.1 : 5];
for i = 1 : size(I.dHPC,1)
    if curve.id.dHPC(i,2) == 1
        figure, title('Up-modulated')
        plot(T,curve.dHPC(i,:),'LineWidth',2,'Color','r'),hold on
        plot(T,shuffle.dHPC(:,1,i),'LineWidth',2,'Color','k')
        ciplot(shuffle.dHPC(:,1,i)-shuffle.dHPC(:,2,i) , shuffle.dHPC(:,1,i)+shuffle.dHPC(:,2,i) , T , 'k'), alpha 0.2
        xline(0) , xline(1)
        ylim([-0.6 3])
        yticks([-0.6 0 0.6 1.2 1.8 2.4 3])
        ylabel('Response (z-score)')
        xlabel('Time from the Shock (sec)')        
    elseif curve.id.dHPC(i,2) == -1
        figure, title('Down-modulated')
        plot(T,curve.dHPC(i,:),'LineWidth',2,'Color','g'),hold on
        plot(T,shuffle.dHPC(:,1,i),'LineWidth',2,'Color','k')
        ciplot(shuffle.dHPC(:,1,i)-shuffle.dHPC(:,2,i) , shuffle.dHPC(:,1,i)+shuffle.dHPC(:,2,i) , T , 'k'), alpha 0.2
        xline(0) , xline(1)
        ylim([-0.6 3])
        yticks([-0.6 0 0.6 1.2 1.8 2.4 3])
        ylabel('Response (z-score)')
        xlabel('Time from the Shock (sec)')
    else
        figure, title('No-modulated')
        plot(T,curve.dHPC(i,:),'LineWidth',2,'Color','y'),hold on
        plot(T,shuffle.dHPC(:,1,i),'LineWidth',2,'Color','k')
        ciplot(shuffle.dHPC(:,1,i)-shuffle.dHPC(:,2,i) , shuffle.dHPC(:,1,i)+shuffle.dHPC(:,2,i) , T , 'k'), alpha 0.2
        xline(0) , xline(1)
        ylim([-0.6 3])
        yticks([-0.6 0 0.6 1.2 1.8 2.4 3])
        ylabel('Response (z-score)')
        xlabel('Time from the Shock (sec)')
    end
end


%vHPC
% For plotting and checking if shuffle works
T = [-5 : 0.1 : 5];
for i = 1 : size(I.vHPC,1)
    if curve.id.vHPC(i,2) == 1
        figure, title('Up-modulated')
        plot(T,curve.vHPC(i,:),'LineWidth',2,'Color','r'),hold on
        plot(T,shuffle.vHPC(:,1,i),'LineWidth',2,'Color','k')
        ciplot(shuffle.vHPC(:,1,i)-shuffle.vHPC(:,2,i) , shuffle.vHPC(:,1,i)+shuffle.vHPC(:,2,i) , T , 'k'), alpha 0.2
        xline(0) , xline(1)
        ylim([-0.6 3])
        yticks([-0.6 0 0.6 1.2 1.8 2.4 3])
        ylabel('Response (z-score)')
        xlabel('Time from the Shock (sec)')        
    elseif curve.id.vHPC(i,2) == -1
        figure, title('Down-modulated')
        plot(T,curve.vHPC(i,:),'LineWidth',2,'Color','g'),hold on
        plot(T,shuffle.vHPC(:,1,i),'LineWidth',2,'Color','k')
        ciplot(shuffle.vHPC(:,1,i)-shuffle.vHPC(:,2,i) , shuffle.vHPC(:,1,i)+shuffle.vHPC(:,2,i) , T , 'k'), alpha 0.2
        xline(0) , xline(1)
        ylim([-0.6 3])
        yticks([-0.6 0 0.6 1.2 1.8 2.4 3])
        ylabel('Response (z-score)')
        xlabel('Time from the Shock (sec)')
    else
        figure, title('No-modulated')
        plot(T,curve.vHPC(i,:),'LineWidth',2,'Color','y'),hold on
        plot(T,shuffle.vHPC(:,1,i),'LineWidth',2,'Color','k')
        ciplot(shuffle.vHPC(:,1,i)-shuffle.vHPC(:,2,i) , shuffle.vHPC(:,1,i)+shuffle.vHPC(:,2,i) , T , 'k'), alpha 0.2
        xline(0) , xline(1)
        ylim([-0.6 3])
        yticks([-0.6 0 0.6 1.2 1.8 2.4 3])
        ylabel('Response (z-score)')
        xlabel('Time from the Shock (sec)')
    end
end