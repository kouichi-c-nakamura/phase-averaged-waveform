function Wadj = getBUA_meanAdjustLocalConditional(bua,indRep)
% Wadj = getBUA_meanAdjustLocalConditional(bua,indRep)
%
% an alternative of getBUA_meanAdjustGlobal1
%
% bua     WaveformChan object, high-pass filtered, spike-removed signal
%         before full wave rectification or mean subtraction
%
% See also
% getBUA_meanAdjustGlobal1,getBUA_meanAdjustGlobal2, getBUA, getBUA_meanAdjustLocal

n = bua.Length;

rectified = bua;
rectified.Data = abs(bua.Data);

jump = find(diff(indRep) > 1);
indRepEpochsStart = [1;jump+1];
indRepEpochsEnd = [jump; length(indRep)];

realRepStart = indRep(indRepEpochsStart); % joined-up considered
realRepEnd   = indRep(indRepEpochsEnd);

Crep = arrayfun(@(x,y) x:y,realRepStart,realRepEnd,'UniformOutput',false);
realRep  = cellfun(@(x) trimExcess(x,n),Crep,'UniformOutput',false);
repMeans = cellfun(@(x) nanmean(rectified.Data(x)),realRep); 

srdWidth = round(0.003/rectified.SInterval);

befStart = realRepStart- srdWidth;
befEnd   = realRepEnd- 1;
Cbef = arrayfun(@(x,y) x:y,befStart,befEnd,'UniformOutput',false);
befMeans = cellfun(@(x) nanmean(rectified.Data(trimExcess(x,n))), Cbef);

aftStart = realRepEnd + 1;
aftEnd   = realRepEnd + srdWidth;
Caft = arrayfun(@(x,y) x:y,aftStart,aftEnd,'UniformOutput',false);
aftMeans = cellfun(@(x) nanmean(rectified.Data(trimExcess(x,n))), Caft);

befAftMeans = mean([befMeans,aftMeans],2);
meanDiff   = befAftMeans - repMeans;
% if meanDiff >= 0; % addition; end
% if meadDiff < 0; % subtracion; end

newdata = rectified.Data;
for i = 1:length(realRep)
    
    ind = realRep{i};
    
    % conditinal
    if min(rectified.Data(ind)) + meanDiff(i) >= 0
        newdata(ind) = rectified.Data(ind) + meanDiff(i);
    else
        newdata(ind) = rectified.Data(ind) - min(rectified.Data(ind));
    end
   %t = rectified.time;figure;plot(t(ind),rectified.Data(ind),'k',t(ind),newdata(ind),'r');
   
end

Wadj = rectified;
Wadj.Data = newdata - mean(newdata); 
%NOTE mean subtraction must be done after mean adjustment
% because the above processing assumes Y >= 0, but signals after mean
% subtraction can take negative values.

Wadj.ChanTitle = [rectified.ChanTitle, ' meanadjusted']; % or Header?

end