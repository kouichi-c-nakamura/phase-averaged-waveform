function Wadj = getBUA_meanAdjustLocal(W,indRep)
% Wadj = local_meanAdjustLocal(W,indRep)
%
% an alternative of getBUA_meanAdjustGlobal1
%
% See also
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

srdWidth = round(0.003/W.SInterval);

befStart = realRepStart- srdWidth;
befEnd   = realRepEnd- 1;
Cbef = arrayfun(@(x,y) x:y,befStart,befEnd,'UniformOutput',false);
befMeans = cellfun(@(x) nanmean(W.Data(trimExcess(x,n))), Cbef);

aftStart = realRepEnd + 1;
aftEnd   = realRepEnd + srdWidth;
Caft = arrayfun(@(x,y) x:y,aftStart,aftEnd,'UniformOutput',false);
aftMeans = cellfun(@(x) nanmean(W.Data(trimExcess(x,n))), Caft);

befAftMeans = mean([befMeans,aftMeans],2);
meanDiff   = befAftMeans - repMeans;


newdata = W.Data;
for i = 1:length(realRep)
    
   ind = realRep{i};
   newdata(ind) = W.Data(ind) + meanDiff(i);
   
end

Wadj = W;
Wadj.Data = newdata;

Wadj.ChanTitle = [W.ChanTitle, ' meanadjusted']; % or Header?

end