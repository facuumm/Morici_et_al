function [cond] = classification_of_asselblies(thresholded,clusters);
% This function was created to reduce the extension of my pipelines.
% They might not easily implemented to other data sets.
% This function return the classification of assemblies (dHPC, vHPC, joint)
%
% Morici Facundo 21/08/2024


if not(isempty(thresholded))
    if size(clusters,1)>0
        cond1 =  sum(thresholded(1:size(clusters,1),:),1)>0; %checking of dHPC SU
        cond2 =  sum(thresholded(size(clusters,1)+1:end,:),1)>0; %checking of vHPC SU
        cond.dHPC = and(cond1 , not(cond2));
        cond.vHPC = and(cond2 , not(cond1));
        cond.both = and(cond1 , cond2); clear cond1 cond2
    else
        cond1 =  logical(zeros(1,size(thresholded,2))); %checking of dHPC SU
        cond2 =  logical(ones(1,size(thresholded,2))); %checking of vHPC SU
        cond.dHPC = and(cond1 , not(cond2));
        cond.vHPC = and(cond2 , not(cond1));
        cond.both = and(cond1 , cond2); clear cond1 cond2
    end
else
    cond1 =  false; %checking of dHPC SU
    cond2 =  logical(0); %checking of vHPC SU
    cond.dHPC = and(cond1 , not(cond2));
    cond.vHPC = and(cond2 , not(cond1));
    cond.both = and(cond1 , cond2); clear cond1 cond2
end

end