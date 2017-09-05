function combScore=getCombinedScore(infoGain,distance,weightInfo,weightDist)
% calculates the combined score ofa region incorporating distance and
% information gain

% standardise infoGain 
infoGain=bsxfun(@minus,infoGain,min(infoGain));
infoGain=infoGain./max(infoGain);
infoGain(isnan(infoGain))=0; % to account for the first snap with all elements zero


distance=bsxfun(@minus,distance,min(distance));
distance=distance./max(distance);
distance(isnan(distance))=0;

combScore=infoGain.*weightInfo+distance.*weightDist;
end