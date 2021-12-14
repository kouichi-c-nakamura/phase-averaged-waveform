function Wadj = getBUA_meanAdjustLocalConditional2(W,indRep,befp,aftp)
% Wadj = getBUA_meanAdjustLocalConditional2(W,indRep)
%
% an alternative of getBUA_meanAdjustGlobal1
%
% W     WaveformChan object, high-pass filtered, spike-removed signal
%         before full wave rectification or mean subtraction
%
% See also
% getBUA_meanAdjustGlobal1,getBUA_meanAdjustGlobal2, getBUA, getBUA_meanAdjustLocal

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
% if meanDiff >= 0; % move up; end
% if meadDiff < 0; % move down; end


%% Core algorithm

newdata = W.Data;
for i = 1:length(realRep)
    
    ind = realRep{i};
    
    if min(W.Data(ind)) + meanDiff(i) < min(W.Data) %TODO
        newdata(ind) = W.Data(ind) - (min(W.Data) - min(W.Data(ind)));
    elseif max(W.Data(ind)) + meanDiff(i) > max(W.Data)
        newdata(ind) = W.Data(ind) + max(W.Data) - max(W.Data(ind));
    else
        newdata(ind) = W.Data(ind) + meanDiff(i);
    end
    
    %NOTE only for debug
%     try
%         
%         t = W.time;
%         thisrange = ind(1)-befp:ind(end)+aftp;
%         plot(t(thisrange),W.Data(thisrange),'k',t(thisrange),newdata(thisrange),'r');
%         pause(2);
%     
%     catch
%         
%     end
   
end

Wadj = W;
Wadj.Data = newdata;
Wadj.ChanTitle = [W.ChanTitle, ' meanadjusted']; % or Header?

end