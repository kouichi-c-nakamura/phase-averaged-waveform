function plotReplacementsVsSurrounds(W,indRep)


% for i=1:length(realRep)
%     ind = realRep{i};
%     t = W.time;
%     thisrange = ind(1)-befp:ind(end)+aftp;
% end


np = W.Length;

jump = find(diff(indRep) > 1);
indRepEpochsStart = [1;jump+1];
indRepEpochsEnd = [jump; length(indRep)];

realRepStart = indRep(indRepEpochsStart); % joined-up considered
realRepEnd   = indRep(indRepEpochsEnd);

Crep = arrayfun(@(x,y) x:y,realRepStart,realRepEnd,'UniformOutput',false);
realRep  = cellfun(@(x) trimExcess(x,np),Crep,'UniformOutput',false);
repMeans = cellfun(@(x) nanmean(W.Data(x)),realRep); 

srdWidth = round(0.003/W.SInterval);

befStart = realRepStart- srdWidth;
befEnd   = realRepEnd- 1;
Cbef = arrayfun(@(x,y) x:y,befStart,befEnd,'UniformOutput',false);
befMeans = cellfun(@(x) nanmean(W.Data(trimExcess(x,np))), Cbef);

aftStart = realRepEnd + 1;
aftEnd   = realRepEnd + srdWidth;
Caft = arrayfun(@(x,y) x:y,aftStart,aftEnd,'UniformOutput',false);
aftMeans = cellfun(@(x) nanmean(W.Data(trimExcess(x,np))), Caft);

befAftMeans = mean([befMeans,aftMeans],2);
meanDiff   = befAftMeans - repMeans;


nrep = length(realRep);

%%
f1 = figure;
g = [repmat({'befAftMeans'},nrep,1);...
    repmat({'repMeans'},nrep,1);...
    repmat({'meanDiff'},nrep,1)];
boxplot([befAftMeans;repMeans;meanDiff],g);

ph = linealpha(repmat((1:3)',1,nrep),[befAftMeans,repMeans,meanDiff]',...
    [0.5 0.5 0.5],0.1);

meanmarks = line(1:3,[mean(befAftMeans),mean(repMeans),mean(meanDiff)],...
    'Marker','o','Color',defaultPlotColors(3));

hbar = line([0.5 3.5], [0 0], 'Color',[0.5 0.5 0.5],'LineStyle','--',...
    'Color','g');
xlim([0.5 3.5])
a = gca;
a.Children = a.Children ([end,1:end-1]);
set(a,'Box','off','TickDir','out')

ylabel('(mV)')
title(sprintf('%s%s | %s',W.Header.animal,W.Header.record,W.Header.title));


%% TODO max(W.Data), min(W.Data), max(W.Data(ind))), min(W.Data(ind))?



