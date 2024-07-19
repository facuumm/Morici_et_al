function [time_projection] = assembly_activity_only_joint(AssemblyTemplates,SpikeCount,template, template1,components)
% I modified the function assembly_activity from Lopes-dos-Santos V et al (2013)
% to only calculate the assemblies activity dirven by cross-structural
% activity.
%
% INPUTS
%
% template and template1= column vector containing 1 or 0 depending on the 
%                         structure they contain.
%                         Example of a template for vHPC
%                                   0    dHPC1
%                                   0    dHPC2
%                                   0    dHPC3
%                                   0    dHPC4
%                                   0    dHPC5
%                                   1    vHPC1
%                                   1    vHPC2
%                                   1    vHPC3
%
% components: string, to select what component of the projector matrix to
%             use. 
%                   'all' --> dorsal+ventral activations
%                   'dHPC' --> dorsal activations
%                   'vHPC' --> ventral activations
%                   'both' --> just joint activations
%
% Please send bug reports to vitor@neuro.ufrn.br (Vítor)
% Modified to exclude parts of the projector matrix (Morici Juan Facundo 02-2024)

% Modification introduced by Morici Juan Facundo
if strcmp(components,'all')
    template = ones(size(template));
    template = (template*template');
elseif strcmp(components,'dHPC')
    template = (template*template');
elseif strcmp(components,'vHPC')
    template = (template1*template1');
elseif strcmp(components,'both')
    template = not((template*template') + (template1*template1'));
end


try
SpikeCount = zscore(SpikeCount')';
catch
    for neuron_i = 1:size(SpikeCount,1);
        SpikeCount(neuron_i,:) = zscore(SpikeCount(neuron_i,:));
    end
end

time_projection=zeros(size(AssemblyTemplates,2),length(SpikeCount));
for assembly_idx = 1:size(AssemblyTemplates,2)
    
    % computing projector
    ASSEMBLYPROJECTOR=AssemblyTemplates(:,assembly_idx)*AssemblyTemplates(:,assembly_idx)';
    ASSEMBLYPROJECTOR=squeeze(ASSEMBLYPROJECTOR)-diag(diag(squeeze(ASSEMBLYPROJECTOR)));
    ASSEMBLYPROJECTOR = ASSEMBLYPROJECTOR.*template;

    % computing activity time course
    for ntime=1:length(SpikeCount)
        
        time_projection(assembly_idx,ntime)=(SpikeCount(:,ntime)'*ASSEMBLYPROJECTOR*SpikeCount(:,ntime));
        
    end
    
end