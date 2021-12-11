function [results, handles] = K_PhaseWave(lfpwaveform,eegwaveform,sourceRate,newRate,b,a,varargin)
% K_PhaseWave is similar to K_PhaseHist, but works for a pair of ECoG and
% LFP waveform signals. It returns values for plotting "phase-triggered
% average waveforms".
%
% [results, handles] = K_PhaseWave(lfpwaveform,eegwaveform,sourceRate,newRate,b,a,varargin)
%
% INPUT ARGUMENTS
% lfpwaveform  vector of source data from LFP channel
%
% eegwaveform  vector of source data from EEG channel (lfpwaveform &
%              eegwaveform must have same length)
%
% sourceRate   sampling rate [Hz] of the input data event and eeg
%
% newRate      the new sampling rate [Hz] after resample
%              In many cases, 1024 is good.
%
% b, a         b and a as coefficients of a filter transfer function
%              b, a must be calculated for newRate rather than souceRate
%              b for numeraotr, a for denominator. You can get b and a by:
%
%              [b, a] = butter(n, Wn)
%
%              where n is the order of the Butterworth filter (lowpass), or
%              half the order(bandpass). Wn is normalized cuttoff
%              frequency, i.e. [cycles per sec] devided by Niquist
%              frequency [newRate/2].
%
%              Wn = frequencyHerz/(samplingrateHerz/2)
%
%              The following command can check the stablity of the filter
%              fvtool(b,a,'FrequencyScale','log', 'Analysis','info');
%
%
% OPTIONAL PARAMETER/VALUE PAIRS
%
% 'PlotLinear'   true | false (default)
%
% 'PlotCirc'     true | false (default)
%
% 'Histbin'      number of histogram bins (default = 72)
%
% 'HistType'     'line' (default) | 'bar'
%
% 'Randomization' 
%                'bootstrap' (default) | 'circshift' | 'none'
%
%
% OUTPUT ARGUMENTS
% results       Structure conttaining following fields
%
%     binmean       Mean per bin
%     binstd        SD per bin
%     binsem        SEM per bin
%     axrad         Phase axis (-pi to pi)
%
%     meanvec       Non-scalar structure with the following fields
%
%         vec       Vector represented as a complex number
%
%                     vec = mean(xy)*2, 
%
%                   where, xy = rad2cmp(x).*y, x is instantaneous phase (in
%                   radian), and y is the value of the waveform data
%                   (typically in mV), Y = rad2cmp(X) is a function defined
%                   by:
%
%                     y = exp(1i * x);
%
%                   where 1i is equivalent to 1*i, and i is the imaginary
%                   unit.
%
%
%         length   abs(vec)
%         radian   angle(vec)
%         degree   rad2deg(angle(vec))
%
%    bootstrap, circshift    
%                 These are results of shuffling by different methods
%         p_lessthan 
%                  Estimated range of P value
%         length   
%         radian
%         degree
%
% handles        Structure of graphic handles
%
%
%
%NOTE
% What about bias correction with ecdf? How does it apply to LFP ddata?
%
% Because here we do not really create histograms of events for Rayleigh
% test, but rather create average waveforms in relation to instantaneous
% phase in essense, even if the phases are non-uniform, average will be
% computed for each bin without bias.
%
%
% See also
% K_PhaseWaveWave_test, K_plotLinearPhaseWave, K_plotCircPhaseWave_one,
% K_plotCircPhaseHist_group,
% K_PhaseHist, ECoGLFPphase_draft, K_plotLinearPhaseHist, K_plotCircPhaseHist_one
% scr2016_06_23_121300_ECoGLFPphase_draft.mlx
% scr2016_07_04_223502_K_PhaseWaveWave_draft.mlx


narginchk(4, inf);

[lfpwaveform, eegwaveform, sourceRate, newRate, b, a, plotLinear, ...
    plotCirc, nbin, ColorSpec, plottype,randomiz] = ...
    local_parse(lfpwaveform, eegwaveform, sourceRate, newRate, b, a, varargin{:});


%% Resampling

if newRate ~= sourceRate
    lfp = resample(lfpwaveform, newRate, sourceRate); % resampled lfp waveform
    clear lfpwaveform
    
    eeg = resample(eegwaveform, newRate, sourceRate); % resampled eeg waveform
    clear eegwaveform
else
    lfp = lfpwaveform;
    eeg = eegwaveform;
end

%% Mean subtraction 
% without this, output potentials can have arbitorary zeros

lfp = lfp - mean(lfp);
eeg = eeg - mean(eeg);


%% core computation for Hilbert transform
eegf = filtfilt(b, a, eeg); % filtering without phase shift
eegrad = angle(hilbert(eegf)); % instantaneous phase value in radians
eegenv = abs(hilbert(eegf)); % amplitude envelope

lfpf = filtfilt(b,a,lfp);
lfprad = angle(hilbert(lfpf));
lfpenv = abs(hilbert(lfpf));


%% Bin by phase 

[results] = local_phaseBinning(eegrad,lfp,nbin);


%% Bootstrap

switch randomiz
    case 'bootstrap'
        %     tic
        %
        %     % for the sake of memory and speed use single here
        %     [rad_btstat,~] = bootstrp(1000, @(x) x, single(eegrad));
        %     [len_btstat,~] = bootstrp(1000, @(x) x, single(lfp));
        %
        %     telap = toc;
        %     fprintf('Bootstrapping took %.1f min\n',telap/60)
        %
        %
        %     btvec = rad2cmp(rad_btstat).*len_btstat;
        %     meanvecSh = mean(btvec,2);
        %
        
        meanvecSh = local_bootstrap(lfp,eegrad);
        [results.bootstrap.p_lessthan] = local_p_lessthan(meanvecSh,results);
        results.bootstrap.length = abs(meanvecSh);
        results.bootstrap.radian = angle(meanvecSh);
        results.bootstrap.degree = rad2deg(angle(meanvecSh));
        
        results.circshift.p_lessthan = [];
        results.circshift.length = [];
        results.circshift.radian = [];
        results.circshift.degree = [];
        
    case 'circshift'
        meanvecSh = local_shuffleByCircshift(lfp,eegrad);
        [results.circshift.p_lessthan] = local_p_lessthan(meanvecSh,results);
        results.circshift.length = abs(meanvecSh);
        results.circshift.radian = angle(meanvecSh);
        results.circshift.degree = rad2deg(angle(meanvecSh));
        
        results.bootstrap.p_lessthan = [];
        results.bootstrap.length = [];
        results.bootstrap.radian = [];
        results.bootstrap.degree = [];

    case 'none'
        meanvecSh = [];
        
        results.bootstrap.p_lessthan = [];
        results.bootstrap.length = [];
        results.bootstrap.radian = [];
        results.bootstrap.degree = [];
        
        results.circshift.p_lessthan = [];
        results.circshift.length = [];
        results.circshift.radian = [];
        results.circshift.degree = [];
        
end

if plotLinear
    
    hlin = K_plotLinearPhaseWave(results,'Plottype',plottype);
    
else
    hlin = [];
end

if plotCirc
    hcirc = K_plotCircPhaseWave_one(results);
    
    switch randomiz
        case 'bootstrap'
            [f4,linh3]= local_plotECDF(meanvecSh,results.meanvec.length);
        case 'circshift'
            %TODO 
            [f4,linh3]= local_plotECDF(meanvecSh,results.meanvec.length);
    end
    
    
else
    hcirc = [];
end

handles.hlin = hlin;
handles.hcirc =hcirc;

end

%--------------------------------------------------------------------------

function [targetwave, refwave, sourceRate, newRate, b, a, plotLinear, ...
    plotCirc, histbin, ColorSpec, histtype,randomiz] = ...
    local_parse(targetwave, refwave, sourceRate, newRate, b, a, varargin)
%
% See also
% K_PhaseHist/local_parse


p = inputParser;

p.addRequired('targetwave', @(x) isnumeric(x) && iscolumn(x));

p.addRequired('refwave', @(x) isnumeric(x) && iscolumn(x) && ...
    length(x) == length(targetwave));

vfscnumpos = @(x) isnumeric(x) && isscalar(x) && x > 0;

p.addRequired('sourceRate', vfscnumpos);
p.addRequired('newRate', vfscnumpos);

vfnumrow = @(x) isnumeric(x) && isrow(x);

p.addRequired('b', vfnumrow);
p.addRequired('a', vfnumrow);

p.addParameter('plotLinear', false, @(x) ~isempty(x) && isscalar(x) && ...
    x == 0 || x == 1);
p.addParameter('plotCirc', false, @(x) ~isempty(x) && isscalar(x) && ...
    x == 0 || x == 1);
p.addParameter('histBin', 72, @(x) isscalar(x) && isnumeric(x) && ...
    x > 0 && fix(x) == x);
p.addParameter('Color', 'b', @(x) iscolorspec(x));
p.addParameter('Threshold', [0 0], @(x) ~isempty(x) && ...
    isnumeric(x) && isrow(x) && length(x) ==2);
p.addParameter('histType', 'line', @(x) ~isempty(x) && ischar(x) && isrow(x) &&...
    ismember(lower(x),{'line','bar'}));
p.addParameter('Randomization','bootstrap', @(x) ...
    ismember(x,{'bootstrap','circshift','none'}));

p.parse(targetwave, refwave, sourceRate, newRate, b, a,varargin{:});

if ~isstable(b, a)
    %if ~K_isstable(b, a) % using fvtool
    warning(eid('filter:notstable'), ...
        'Filter is not stable. Reconsider the parameters.');
end

plotLinear = logical(p.Results.plotLinear);
plotCirc   = logical(p.Results.plotCirc);
histbin    = p.Results.histBin;
ColorSpec  = p.Results.Color;


histtype   = lower(p.Results.histType);
randomiz    = p.Results.Randomization;

end

%--------------------------------------------------------------------------

function y = rad2cmp(x)

y = exp(1i * x);

end

%--------------------------------------------------------------------------

function S = local_phaseBinning(x,y,nbin)

binrad = (2*pi)/nbin;
edges = (-pi:binrad:pi)';

axrad = edges(1:end-1) + binrad/2; % + binrad/2 to point the center of each bin

binind = discretize(x,edges);
n = length(edges) - 1;
binmean = zeros(n,1);
binstd  = zeros(n,1);
binsem  = zeros(n,1);

for i = 1:n
    ind = find(binind == i);
    binmean(i) = nanmean(y(ind));
    binstd(i) = nanstd(y(ind));
    binsem(i) = binstd(i)/sqrt(nnz(ind));
end

S.binmean = binmean;
S.binstd  = binstd;
S.binsem  = binsem;
S.axrad   = axrad;

%%

instPhaseCmp = rad2cmp(x);

xy = instPhaseCmp .* y;

meanvec = mean(xy)*2; %NOTE The length of mean(xy) *2 is roughly equivalent to the peak amplitude
% See also 
% local_bootstrap local_circshift

S.meanvec.vec = meanvec;
S.meanvec.length = abs(meanvec);
S.meanvec.radian = angle(meanvec);
S.meanvec.degree = rad2deg(angle(meanvec));

end

%--------------------------------------------------------------------------

function [f4,linh3]= local_plotECDF(meanvecSh,L)


%% p values and other outputs

[f,x]= ecdf(abs(meanvecSh));
f4 = figure;
axh4 = gca;
ax4.LabelFontSizeMultiplier = 1.1;
plot(x,f)
hold on
plot([L,L],[0 1])
text(L,0.53,'Real Data')

plot([max(abs(meanvecSh)),max(xlim)],[1,1],'Color',defaultPlotColors(1));

plot(xlim,[0.95,0.95],'--','Color',defaultPlotColors(2));
plot(xlim,[0.99,0.99],'--','Color',defaultPlotColors(3));
plot(xlim,[0.999,0.999],'--','Color',defaultPlotColors(4));

set(axh4,'Box','off','TickDir','out')
xlabel('Vector Length (mV)')
ylabel('Cumulative Distribution Function')


the95 = prctile(abs(meanvecSh),95);
the99 = prctile(abs(meanvecSh),99);
the999 = prctile(abs(meanvecSh),99.9);

linh3(1) = plot([the95,the95],[0 1],'--','Color',defaultPlotColors(2),...
    'DisplayName','p = 0.05');
linh3(2) = plot([the99,the99],[0 1],'--','Color',defaultPlotColors(3),...
    'DisplayName','p = 0.01');
linh3(3) = plot([the999,the999],[0 1],'--','Color',defaultPlotColors(4),...
    'DisplayName','p = 0.001');

leg = legend(linh3);
leg.Location = 'southeast';
leg.Box = 'off';

end

%--------------------------------------------------------------------------

function [meanvecSh] = local_bootstrap(lfp,eegrad)

tic

% for the sake of memory and speed use single here
[rad_btstat,~] = bootstrp(1000, @(x) x, single(eegrad));
[len_btstat,~] = bootstrp(1000, @(x) x, single(lfp));

elap = toc;
% fprintf('Bootstrapping took %.1f min\n',elap/60)

btvec = rad2cmp(rad_btstat).*len_btstat;
meanvecSh = mean(btvec,2)*2; %NOTE *2

end

%--------------------------------------------------------------------------

function meanvecSh = local_shuffleByCircshift(lfp,eegrad)

cycleStarts = find(diff(eegrad) < -6)+1;
Cy = length(cycleStarts);

instPhaseC = repmat({eegrad},1000,1);%TODO consider using single
meanvecSh = zeros(1000,1);

tic

parfor j = 1:1000
    
    cycleStarts_ = cycleStarts; % broadcast to temporal variable
    instPhase_ = eegrad; % broadcast to temporal variable
    instPhaseCj = instPhaseC{j}; % broadcast to temporal variable
    
    if round(j/100) > round((j-1)/100)
        fprintf('*')
    end
    for i = 1:Cy-1
        thisrange = cycleStarts_(i):cycleStarts_(i+1)-1;
        jump = randi(length(thisrange)) - 1;
        
        newInd = circshift(thisrange,jump);
        instPhaseCj(newInd) = instPhase_(thisrange);
        
    end
    
    instPhaseC{j} = instPhaseCj;
    
    instPhaseCmp0 = rad2cmp(instPhaseCj);
    
    xy = instPhaseCmp0 .* lfp;
    
    meanvecSh(j) = mean(xy)*2; %NOTE *2
end
elap = toc;

% fprintf('randi and circshift took %0.1f min\n',elap/60)


end

%--------------------------------------------------------------------------

function [p_lessthan] = local_p_lessthan(meanvecSh,results)
% [p_lessthan] = local_p_lessthan(meanvecSh,results)
%
%
% p_lessthan     Not exactly p value itself, but estimation of the range of 
%                p values based on random resampling.

[f,x]= ecdf(abs(meanvecSh));

meanvec = results.meanvec.vec;
L = abs(meanvec);
ind = find(L > x,1,'last');
if ind ~= length(f)
    p_lessthan = 1-f(ind);
else
    p_lessthan = 1 - f(end-1);
end

end
