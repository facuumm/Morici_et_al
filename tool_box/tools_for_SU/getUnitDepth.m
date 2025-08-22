function [depth_um, depth_norm, x_um] = getUnitLocation(template, chanPos)
% getUnitLocation computes the anatomical location of a neuron based on its spike templates
%
% INPUTS:
%   template : [nTemplates x nTime x nChannels] matrix of spike templates
%              or [1 x nTime x nChannels] for a single template
%   chanPos  : [nChannels x 3] matrix:
%                 column 1 = channel ID
%                 column 2 = y-position (depth, um)
%                 column 3 = x-position (um)
%
% OUTPUTS:
%   depth_um   : estimated depth of the neuron in micrometers
%   depth_norm : depth normalized between 0 (superficial) and 1 (deep)
%   x_um       : estimated lateral (x) position in micrometers
%
% Notes:
%   - Works for positive or negative spikes.
%   - If multiple templates are given for the same neuron, averages them first.

    % Check input dimensions
    if ndims(template) ~= 3
        error('Template must be 3D: [nTemplates x nTime x nChannels]');
    end
    
    % Average across templates if more than one
    if size(template,1) > 1
        templateAvg = squeeze(mean(template,1));  % [nTime x nChannels]
    else
        templateAvg = squeeze(template);           % [nTime x nChannels]
    end
    
    % Compute peak-to-peak amplitude for each channel
    w = max(templateAvg,[],1) - min(templateAvg,[],1);  % 1 x nChannels
    w = w(:);  % make column vector
    
    % Extract channel positions
    xPos = double(chanPos(:,2));  % depth
    yPos = double(chanPos(:,3));  % lateral
    
    if length(yPos) ~= length(w)
        error('Number of channels in template and chanPos must match.');
    end
    
    % Weighted depth and lateral position
    depth_um = sum(w .* yPos) / sum(w);
    x_um     = sum(w .* xPos) / sum(w);
    
    % Normalized depth between 0 (superficial) and 1 (deep)
    depth_norm = (depth_um - min(yPos)) / (max(yPos) - min(yPos));
end