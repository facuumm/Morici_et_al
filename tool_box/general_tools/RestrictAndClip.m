function samples = RestrictAndClip(samples,intervals,varargin)
%RestrictAndClip - Keep only samples that fall in a list of time intervals,
% trimming the intervals to fit within the sample range.
%
%    samples = RestrictAndClip(samples,intervals,<options>)
%
%    samples        samples to restrict (Nx1 or NxM matrix, first column is time)
%    intervals      list of (start,stop) pairs
%    <options>      optional list of property-value pairs
%
%    -------------------------------------------------------------------------
%     'shift'       shift remaining epochs together in time (default = 'off')
%    -------------------------------------------------------------------------
%
%  This function trims each interval to fall within the bounds of the sample data.
%  Intervals fully outside the sample range are discarded.
% Adaptation from Restrict function from FMA toolbox, to modify the
% begining of the interval if is outside of the intervals.
% Facundo Morici 07/2025

% Default values
shift = 'off';

% Check number of parameters
if nargin < 2 || mod(length(varargin),2) ~= 0
    error('Incorrect number of parameters (type ''help RestrictAndClip'' for details).');
end

% Check parameters
if ~isdmatrix(intervals) || size(intervals,2) ~= 2
    error('Incorrect intervals (type ''help RestrictAndClip'' for details).');
end

if size(samples,1) == 1
    samples = samples(:);
end

% Parse parameter list
for i = 1:2:length(varargin)
    if ~ischar(varargin{i})
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help RestrictAndClip'' for details).']);
    end
    switch(lower(varargin{i}))
        case 'shift'
            shift = varargin{i+1};
            if ~isastring(shift,'on','off')
                error('Incorrect value for property ''shift'' (type ''help RestrictAndClip'' for details).');
            end
        otherwise
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help RestrictAndClip'' for details).']);
    end
end

% -------------------------------
% Trim intervals to sample range
% -------------------------------
sample_start = min(samples(:,1));
sample_end = max(samples(:,1));

new_intervals = [];

for i = 1:size(intervals,1)
    int_start = intervals(i,1);
    int_end = intervals(i,2);

    if int_end <= sample_start || int_start >= sample_end
        continue;
    end

    clipped_start = max(int_start, sample_start);
    clipped_end = min(int_end, sample_end);

    new_intervals = [new_intervals; clipped_start, clipped_end];
end

intervals = new_intervals;

% -------------------------------
% Restrict
% -------------------------------
if isempty(intervals)
    samples = samples([],:);  % return empty with same shape
    return;
end

[status,interval,index] = InIntervals(samples,intervals);
samples = samples(status,:);

% -------------------------------
% Shift (optional)
% -------------------------------
if strcmp(shift,'on')
    interval = interval(status);
    start = intervals(2:end,1);
    stop = intervals(1:end-1,2);
    cumulativeShift = [0; cumsum(start - stop)];
    shifts = cumulativeShift(interval);
    samples(:,1) = samples(:,1) - shifts - intervals(1,1);
end