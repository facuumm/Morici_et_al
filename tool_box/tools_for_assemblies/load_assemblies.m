function [patterns , cond , Thresholded] = load_assemblies(p , name , clusters, numberD , tag)

load([p,'\',name])

thresholded = Th;
p = pat;
clear cond Th pat

p = p .* thresholded;

% Detection of members
if not(isempty(p))
    if numberD>0
        cond1 =  sum(thresholded(1:size(clusters.dHPC,1),:),1)>0; %checking of dHPC SU
        cond2 =  sum(thresholded(size(clusters.dHPC,1)+1:end,:),1)>0; %checking of vHPC SU
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

Thresholded.(tag) = thresholded;
patterns.(tag) = p;
end