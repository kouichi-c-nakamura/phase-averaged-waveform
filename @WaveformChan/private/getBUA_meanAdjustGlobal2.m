function Wadj = getBUA_meanAdjustGlobal2(W,spikeEvent,befp,aftp)
% an alternative of getBUA_meanAdjustGlobal1
% Performance is almost identical
%
% See also
% getBUA_meanAdjustGlobal1

n = W.Length;
spkInd = find(spikeEvent);

repStart = spkInd - befp;
repEnd   = spkInd + aftp;
indRep = startsends2ind(repStart,repEnd,n);
meanRep   = nanmean(W.Data(indRep));


srdWidth = round(0.003/W.SInterval);

befStart = repStart- srdWidth;
befEnd   = repStart- 1;
indBef = startsends2ind(befStart,befEnd,n);
meanBef   = nanmean(W.Data(indBef)) ;

aftStart = repEnd + 1;
aftEnd   = repEnd + srdWidth;
indAft = startsends2ind(aftStart,aftEnd,n);
meanAft   = nanmean(W.Data(indAft));


meanBefAft = mean([meanBef,meanAft]);

meanDiff   = meanBefAft - meanRep;

Wadj = W;
Wadj.Data(indRep) = W.Data(indRep) + meanDiff; % correction of mean of replaced

Wadj.ChanTitle = [W.ChanTitle, ' meanadjusted']; % or Header?

end