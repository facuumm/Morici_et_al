function [cooridnated_event , coordinatedX , coordinatedY] = Find_Coordinated_events(x,y,win)
% This function detect coordinated events coming from two matrices of
% timestamps. It use the peaks timestamps to determine if two events are
% cooridnated.
%
% syntax:
% [cooridnated_event , coordinatedX , coordinatedY] = Find_Coordinated_events(x,y,win)
%
% --- INPUTS ---
% x and y: matrices of timestamps respecting following shape [events x timeStamps]
%           Begining    Peak     End
%            Time1      Time1   Time1
%            Time2      Time2   Time2
%             ...        ...     ...
%            TimeN      TimeN   TimeN
%
% win: float, define the time window before and after the reference to
%      determine if and event in the other matrix is cooridnated.
%
% --- OUTPUTS ---
% coordinated_event: matrices contianing timestamps of merged events
%                    coordinated coming from X and Y.
%
% cooridnatedX and coordinatedY: matrices containing the same data but from
%                                the events of x and y matrices that are
%                                coordinated.
%
% Morici Juan Facundo 10/2024


coordinatedX = [];
coordinatedY = [];
cooridnated_event = [];
for i = 1:length(x)
    r = x(i,:);
    tmp = sum(and(y(:,2)>= r(1,2)-win, y(:,2)<= r(1,2)+win));
    if tmp>0
        z = y(and(y(:,2)>= r(1,2)-win, y(:,2)<= r(1,2)+win),:);
        coordinatedY = [coordinatedY ; z];
        [p,indice] = min(abs(r(2)-z(:,2)));
        coordinatedX = [coordinatedX ; r];
        
        peak = min(min(z(indice,2),r(2)))+abs(z(indice,2)-r(2))/2;
        low = min([r(1) , z(indice,1)]);
        up = max([r(3) , z(indice,3)]);
        cooridnated_event = [cooridnated_event ; low , peak , up];
        
        clear tmp2 tmp1 p indice z peak low up
    end
    clear r
end
clear x tmp i

[C,IA,IC] = unique(coordinatedY(:,1));
coordinatedY  = coordinatedY(IA,:); clear C IA IC
end