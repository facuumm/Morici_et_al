function [events] = FindReplay(Signal,Threshold,Durations,Points,restrict)

%FindRipples - Find Replays events in MUA vector.
%
%  USAGE
%
%    [events] = FindReplay(signal,Threshold,Durations,restrict)
%
%    Putative replay events are detected by thresholding the vector,
%    detecting the start, peak, and stop of each event.
%
%    signal         MUA, counts of all the SUs from your recording
%    Threshold      [start/stop Threshold , peak Threshold] in SD
%    Durations      [minInter , minDuration , maxDuration] in seconds
%    Points         Int value, define the gaussian window size for smooth
%    restrict       To restrict your detection within a perdiod of time
%                   I recomend using it, if not is too time consuming
%
%  OUTPUT
%
%    events        for each replay, [start_t peak_t end_t]
%
% Morici Juan Facundo 06/2023


% Default values
lowThresholdFactor = Threshold(1);
highThresholdFactor = Threshold(2);
minInter = Durations (1);
minDuration = Durations (2);
maxDuration = Durations (3);

% Normalize signal
time = Signal(:,1);
signal = Signal(:,2);
w = gausswin(Points, 1);
signal = filter(w , 1 ,signal);
signal = zscore(signal);

% Detect replays periods by thresholding normalized squared signal
thresholded = signal > mean(signal) + (std(signal) * lowThresholdFactor);
start = find(diff(thresholded)>0);
stop = find(diff(thresholded)<0);

% Exclude last replay if it is incomplete
if length(stop) == length(start)-1,
	start = start(1:end-1);
end

% Exclude first replay if it is incomplete
if length(stop)-1 == length(start),
    stop = stop(2:end);
end

% Correct special case when both first and last replays are incomplete
if start(1) > stop(1),
	stop(1) = [];
	start(end) = [];
end
firstPass = [start,stop];
if isempty(firstPass),
% 	disp('Detection by thresholding failed');
	return
else
% 	disp(['After detection by thresholding: ' num2str(length(firstPass)) ' events.']);
end

% Merge replays if inter-replay period is too short (unless this would yield too long a ripple)
secondPass = firstPass;
iri = time(secondPass(2:end,1)) - time(secondPass(1:end-1,2));
duration = time(secondPass(2:end,2)) - time(secondPass(1:end-1,1));
toMerge = iri<minInter & duration<maxDuration;
while any(toMerge),
    % Get indices of first ripples in pairs to be merged
    rippleStart = strfind([0 toMerge'],[0 1])';
    % Incorporate second ripple into first in all pairs
    rippleEnd = rippleStart+1;
    secondPass(rippleStart,2) = secondPass(rippleEnd,2);
    % Remove second ripples and loop
    secondPass(rippleEnd,:) = [];
    iri = time(secondPass(2:end,1)) - time(secondPass(1:end-1,2));
    duration = time(secondPass(2:end,2)) - time(secondPass(1:end-1,1));
    toMerge = iri<minInter/1000 & duration<maxDuration/1000;
end


if isempty(secondPass),
% 	disp('Replay merge failed');
	return
else
% 	disp(['After replay merge: ' num2str(length(secondPass)) ' events.']);
end

% Discard replays with a peak power < highThresholdFactor
thirdPass = [];
peakPosition = [];
formax = signal > mean(signal) + (std(signal) * highThresholdFactor);
for i = 1:size(secondPass,1)
    [~,maxValue] = max(signal([secondPass(i,1):secondPass(i,2)]));
    if formax(secondPass(i,1)+maxValue)
        peakPosition = [peakPosition ; secondPass(i,1)+maxValue];
        thirdPass = [thirdPass ; secondPass(i,:)];
    end
end

if isempty(thirdPass),
	disp('Peak thresholding failed.');
	return
else
% 	disp(['After peak thresholding: ' num2str(length(thirdPass)) ' events.']);
end


% Discard ripples that are way too short
events = [time(thirdPass(:,1)) time(peakPosition) time(thirdPass(:,2))];
duration = events(:,3)-events(:,1);
events(duration<minDuration,:) = [];
% disp(['After min duration test: ' num2str(size(events,1)) ' events.']);

% Discard ripples that are way too long
duration = events(:,3)-events(:,1);
events(duration>maxDuration,:) = [];
% disp(['After max duration test: ' num2str(size(events,1)) ' events.']);

if exist('restrict','var')
events = Restrict(events,restrict);
% disp(['After Restriction to selected periods: ' num2str(size(events,1)) ' events.']);
end

end