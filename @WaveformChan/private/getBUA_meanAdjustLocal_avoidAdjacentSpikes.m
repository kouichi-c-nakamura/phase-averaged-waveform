function Wadj = getBUA_meanAdjustLocal_avoidAdjacentSpikes(W,indRep)
% Wadj = getBUA_meanAdjustLocal_avoidAdjacentSpikes(W,indRep)
%
% an alternative of getBUA_meanAdjustGlobal1
%
% See also
% getBUA_meanAdjustLocal
% getBUA_meanAdjustGlobal1,getBUA_meanAdjustGlobal2, getBUA

n = W.Length;

jump = find(diff(indRep) > 1);
indRepEpochsStart = [1;jump+1];
indRepEpochsEnd = [jump; length(indRep)];

realRepStart = indRep(indRepEpochsStart); % joined-up considered
realRepEnd   = indRep(indRepEpochsEnd);

Crep = arrayfun(@(x,y) x:y,realRepStart,realRepEnd,'UniformOutput',false);
realRep  = cellfun(@(x) trimExcess(x,n),Crep,'UniformOutput',false);
repMeans = cellfun(@(x) nanmean(W.Data(x)),realRep); 

%TODO logical vector for real replacements
% before and after is locally detemined by skipping realRep data points
% with for loop
% Compare the results with getBUA_meanAdjustLocal

tforiginal = true(W.Length,1);
tforiginal([realRep{:}]) = false;


srdWidth = round(0.003/W.SInterval);

befEnd   = realRepEnd- 1;

aftStart = realRepEnd + 1;

% to avoid "Index exceeds matrix dimensions." error, discard illegal window
ind1 = befEnd < 1;
ind2 = aftStart > W.Length;
befEnd(ind1) = [];
befEnd(ind2) = [];
aftStart(ind1) = [];
aftStart(ind2) = [];
realRep(ind1) = [];
realRep(ind2) = [];

befAftMeans = zeros(length(realRep),1);
meanDiff = zeros(length(realRep),1);
newdata = W.Data;

for i = 1:length(realRep)
    
    befInd = find(tforiginal(1:befEnd(i)),srdWidth,'last'); %NOTE find is slow
    if isempty(befInd)
        befInd = 1;
    end
    befMean = mean(W.Data(befInd));

    aftInd = find(tforiginal(1:aftStart(i)),srdWidth,'first');% SLOW
    if isempty(aftInd)
        aftInd = n;
    end
    aftMean = mean(W.Data(aftInd));

    befAftMeans(i) = mean([befMean,aftMean]);
    meanDiff(i) = befAftMeans(i) - repMeans(i);
    
    ind = realRep{i};
    newdata(ind) = W.Data(ind) + meanDiff(i);
   
end

Wadj = W;
Wadj.Data = newdata;

Wadj.ChanTitle = [W.ChanTitle, ' meanadjusted']; % or Header?

end