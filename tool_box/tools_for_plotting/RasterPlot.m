function RasterPlot(Times1,Times2,w)
% This function construct a rasterplot using Times1 to lock the occurrence 
% of Times2
%
% INPUTS
% Times1 / Times2: Column vector. Times1 = reference / Times2 = referenced.
% w = float, duration fo the time window sourrounding reference.
%
% OUTPUT
% Raster Plot, in y-axis Times1 events, y x-axis time.
%
% Morci Juan Facundo 11/2023

figure
for i = 1 : length(Times1)
    period = [Times1(i)-w Times1(i)+w];
    spks = Restrict(Times2,period);
    spks = spks-Times1(i);
    spks = spks';
    xspikes = repmat(spks,3,1);
    yspikes = nan(size(xspikes));
    
    if not(isempty(yspikes))
        yspikes(1,:) = i-1;
        yspikes(2,:) = i;
    end
    
    plot(xspikes,yspikes,'Color','k','LineWidth',2),hold on
    clear period spks xspikes yspikes
end
xlim([-w w])

end

