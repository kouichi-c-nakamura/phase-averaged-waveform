function [fwrbuaLD,fwrbuaL,fwrbua,bua,fwrbuanan,spiketimes,figh,width_ms,...
    spikewindow,indRep,indRepD, filtered] = getBUA(wave,from,varargin)
% getBUA returns continuous BUA signal out of wideband LFP signal,
% according to Moran A, Bergman H, Israel Z, Bar-Gad I (2008) Subthalamic
% nucleus functional organization revealed by parkinsonian neuronal
% oscillations and synchrony. Brain 131:3395-3409.
%
%
%   _____ = getBUA(wave,'highpass',spiketimes)
%
%   _____ = getBUA(wave,'highpass thresholding',order,spikeThreshold)
%
%   _____ = getBUA(wave,'wideband',order,spiketimes)
%
%   _____ = getBUA(wave,'wideband thresholding',order,spikeThreshold)
%
%   _____ = getBUA(_____, spikewindow)
%   _____ = getBUA(_____, 'Name',Value)
% 
%   [fwrbuaLD,fwrbuaL,fwrbua,bua,fwrbuanan,spiketimes,figh,width_ms,...
%       spikewindow,indRep,indRepD] = ...
%       getBUA(_____)
%
% Preparation of FWR BUA is done in the following order:
%
% 1. High-pass filtering >300 Hz with butterworth filter (skipped if |from|
% is 'highpass')
%
% 2. Spike removal with createBUA() with spikewindow option that specifies
% the window size
%
% 3. Rectification (abs(x)) and mean subtraction (x - mean(x))
%
%
% *INPUT ARGUEMENTS*
% wave        A WaveformChan object
%
% from       'highpass' | 'highpass thresholding' | 'wideband'| 'wideband thresholding'
%             If 'highpass' is chosen, |wave| is considered already high
%             pass-filtered and the filtering step is skipped. Spike
%             removal is done with spiketimes. 
%
%             If 'wideband' is chosen, getBUA first high-pass filters
%             |wave| and then creates BUA. Spike removal is done with
%             spiketimes.
%
%             If 'wideband thresholding' is chosen, after filtering, getBUA
%             carries out spike detection according to |spikeThreshold|
%             using |findpeaks()|.
%
% spiketimes  A column of spike times or a MetaEventChan object with
%             matching Start, SRate and Length
%
% order       A positive integer for the order of butterworth filter.
%             For example, when sampling rate is 17 KHz, order = 9 is OK
%             (isstable(b,a) == 1), but can be smaller. Larger order
%             results in more effective filtering, but with drawbacks of
%             creating artefacts. Use fvtool(b,a,'Fs',samplingrate) for
%             validation
%
% spikeThredhold
%             Real value for detecting peak locations with findpeaks'
%             'MinPeakHeight' option. If positive or 0, then findpeaks will
%             be called to detect peaks. If negative, then findpeaks will
%             be called to detect troughs by using negative values of
%             signals.
%
%
% *OPTION*
% spikewindow Two element row vector of zero or positive numbers.
%             [before_ms, after_ms]
%             where before_ms and after_ms are durations relative to a
%             spike time in millisec during which data points are replaced
%             with randomly chosen background activity. See ceateBUA() for
%             more details.
%             Default [0.5, 2.5]
%
%
% *OPTIONAL PARAMETER/VALUE PAIRS*
%
% 'fvtool'    'on' | 'off' (default)
%             'on' to use fvtool to visually validate the highpass filter
%
% 'plotBW'    'on' | 'off' (default)
%             Show where spike removal occurred.
% 
% 'plotFindpeaks'
%             'on' | 'off' (default)
%             Show where detected (and to be removed) peaks are.
%
% 'plotTriggered'
%             'on' | 'off' (default)
%             Show waveform overdraw for the original waveform data and FWR
%             BUA signal relative to spiketimes.
%
% 'plotWith'  'overdraw' | 'std' (default)  | 'sem' | 'averageonly'
%             Overdraw plots always include average waveform, but you can
%             choose whant to show together with the average.
%
% 'newFs'     default 1024 (in Hz)
%             New sampling rate for down sampling to prepare fwrbuaLD
%
% 'overdrawWidth_sec'
%             Width in second for wide overdraw view of fwrbuaLD and fwrbuaL.
%
% 'meanadjustment'   
%             'global' | 'local' | 'localavoidneighbour' | 'off'  
%             | 'on' (default)
%             The mean value of insert fragments will be adusted to level
%             with the mean value of the surrouds before and after detected
%             spiketimes (3 msec window each)
%
%             'on' is equivalent to 'global'
%             'localconditional' gives generally the best results. For the
%             sake of comparisons and validations, othere options are kept
%             available.
%
%             See also
%             getBUA_meanAdjustGlobal1
%             getBUA_meanAdjustGlobal2
%             getBUA_meanAdjustLocal
%             getBUA_meanAdjustLocalConditional
%
%
% *OUTPUT ARGUMENTS*
% fwabuaLD    A WaveformChan object. Low-pass filtered and down sampled (to
%             newFs) fwrbua.
%
% fwrbuaL     A WaveformChan object. Low-pass filtered (<300 Hz) fwrbua.
%
% fwrbua      A WaveformChan object holding FWR BUA (full width rectified
%             background unit activity) of the original data. In order to
%             retrieve the envelope of background spiking activity at a low
%             frequency range, you ought to further down sample and/or low-pass
%             filter the FWR BUA data. ChanTitle is '*_FWRBUA', where * is
%             the original ChanTitle.
%
% bua         A WaveformChan object holding BUA, which is free of large 
%             amplitude spikes but has not yet processed for rectification or
%             mean subtraction. In order to retrieve the envelope of
%             background spiking activity at a low frequency range, you
%             ought to down sample and/or low-pass filter the FWR BUA data.
%             ChanTitle is '*_BUA' where * is the original ChanTitle.
%
% fwrbuanan   A WaveformChan object holding FWR BUA of the original data.
%             Instead of replacing fragments around spike event with
%             background, they were padded with NaN. fwrbuanan can be
%             useful in order to assess the effect of spike removal.
%
% figh        Structure of figure handles with fields:
%                 fvtool
%                 peaks
%                 overdraw0
%                 overdraw1
%                 overdraw2
%                 overdraw3
%                 overdraw4
%                 overdraw5
%                 BW
%
% width_ms    Actual width of spike removal windows in milliseconds
%             considering exention by overlaps
%
% spikewindow Ideal spike removal window in the format
%             [before_ms, after_ms]
%
% indRep,indRepD
%             Indices for replaced data points and those for down sampled
%             data, respectively
%
% filtered    Highpass filtered data.
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% kouichi.c.nakamura@gmail.com
% 23-Jun-2016 22:23:59
%
% See also
% createBUA2 (written by Izhar Bar-Gad, modified by Kouichi Nakamura), 
% eventTriggeredAverage, normalizedfreq, timestamps2binned, getBUA_meanAdjustGlobal1
% getBUA_meanAdjustGlobal2, getBUA_meanAdjustLocal,
% getBUA_meanAdjustLocalConditional


[wave,from,n,spkthre,spiketimes,spikewindow,fvtoolmode,plotBW,plotFindpeaks,...
    plotTriggered,showOverdraw,showErrorbar,newFs,overdrawWidth_sec,...
    meanadjustment] = parse(wave,from,varargin{:});


%% High-pass filtering (>300 Hz) of the broadband (LFP) signal,

figh.fvtool = [];
figh.overdraw0 = [];
figh.overdraw1 = [];
figh.overdraw2 = [];
figh.overdraw3 = [];
figh.BW = [];

switch from
    case {'wideband','wideband thresholding','highpass thresholding'}
        
        switch from
            case {'wideband','wideband thresholding'}
                Wn = normalizedfreq(300,wave.SRate);
                
                [b,a] = butter(n, Wn,'high'); % the dimension n can be empirically explored
                if strcmpi(fvtoolmode,'on')
                    fvtool(b,a,'Fs',wave.SRate);
                    xlim([0, 0.5]);
                end
                
                assert(isstable(b,a));
                filtered = filtfilt(b, a, wave.Data);
            case 'highpass thresholding'
                filtered = wave.Data;
        end
        
        %%  Thresholding with findpeaks
        if strcmp(from,'wideband thresholding') || strcmp(from,'highpass thresholding')
            
            wh = which('findpeaks'); %NOTE solve the conflict with Chronux's findpeaks
            if ~isempty(wh) && ~contains(wh, ...
                    fullfile('toolbox','signal','signal','findpeaks.m'))
                rmpath(fileparts(wh));
            end
            
            warning off 'signal:findpeaks:largeMinPeakHeight'
            if spkthre >= 0
                [~,locs] = findpeaks(filtered,'MinPeakHeight',spkthre);
                mk = 'v';
            else
                [~,locs] = findpeaks(-filtered,'MinPeakHeight',-spkthre);
                mk = '^';
            end
            t = wave.time;
            warning on 'signal:findpeaks:largeMinPeakHeight'
            
            spiketimes = t(locs);
            
            if strcmp(plotFindpeaks,'on')
                if ~isempty(spiketimes)
                    
                    wave.plot;
                    set(findobj(gca,'Type','line'),...
                        'DisplayName',wave.ChanTitle);
                    
                    title(['Peaks: ',wave.ChanTitle]);

                    tvec = wave.time;
                    line(tvec(locs),filtered(locs),'Marker',mk,...
                        'LineStyle','none',...
                        'MarkerEdgeColor',defaultPlotColors(1),...
                        'MarkerFaceColor',defaultPlotColors(1),...
                        'DisplayName','Peaks/Troughs');
                    
                    line([wave.Start,wave.MaxTime],[spkthre,spkthre],...
                        'LineStyle','--',...
                        'Color',[0.5 0.5 0.5],...
                        'DisplayName','threshold');
                    
                    figh.peaks = gcf;
                    figh.peaks.Tag = 'findpeaks';
                else
                    figh.peaks = gobjects(1);
                end
            else
                figh.peaks = gobjects(1);
            end
          
            figh = local_plotTriggered(wave,wave,spiketimes,showOverdraw,...
                showErrorbar,spikewindow,figh,plotTriggered,'overdraw0','Original',...
                0.03,0.015);

            filtd = WaveformChan(filtered,wave.Start,wave.SRate,wave.ChanTitle);
            figh = local_plotTriggered(filtd,wave,spiketimes,showOverdraw,...
                showErrorbar,spikewindow,figh,plotTriggered,...
                'overdraw1','High-pass filtered',0.03,0.015);


            addpath(fileparts(wh));
        end
        
    case 'highpass'
        filtered = wave;
end

%% Replace the high amplitude units with data from elsewhere in the signal
% Moran A, Bergman H, Israel Z, Bar-Gad I (2008) Subthalamic
% nucleus functional organization revealed by parkinsonian neuronal
% oscillations and synchrony. Brain 131:3395-3409.

if isa(spiketimes,'MetaEventChan')
    
    spikeEvent = spiketimes.Data;
    spiketimes = spiketimes.TimeStamps;
    
else
    
    [spikeEvent, spiketimes] = timestamps2binned(spiketimes,wave.Start,wave.MaxTime,wave.SRate);
    
end

befp = round(wave.SRate/1000*spikewindow(1));
aftp = round(wave.SRate/1000*spikewindow(2));

if isempty(spiketimes)
    BUA = filtered;
    if strcmp(plotTriggered, 'on') || strcmpi(plotBW,'on') || strcmp(plotFindpeaks,'on')
        fprintf('No spike over %f has been detected\n',spkthre)
    end
    befp = [];
    aftp = [];
else
    [BUA,~,~,~,befp2,aftp2] = createBUA2(filtered, spiketimes, wave.SRate, 0, spikewindow);
    BUA = BUA';
    %NOTE spikeTrain (createBUA2's 3rd output) and find(spikeEvent) does not
    % match!
    % The latter is smaller by befp - 1, i.e. appears to taking the befp
    % rather than spike position
    
    %Note spikeEvent equals to original createBUA's spikeTrain -1 (different way of binning time)
    %
    %
    
    assert(befp == befp2)
    assert(aftp == aftp2)
    clear befp2 aftp2
end

bua = wave;
bua.Data = BUA;
bua.ChanTitle = [bua.ChanTitle, '_BUA'];

indRep = local_getIndicesOfReplaced(wave.Length,spikeEvent,befp,aftp);


%% Rectification and mean subtraction

FWRBUA = abs(BUA) - mean(abs(BUA)); %NOTE mean of abs(BUA) (Moran et al 2008)

fwrbua = wave;
fwrbua.Data = FWRBUA;
fwrbua.ChanTitle = [fwrbua.ChanTitle, '_FWRBUA'];


[figh,xdata,ydata] = local_plotTriggered(fwrbua,wave,spiketimes,showOverdraw,showErrorbar,...
    spikewindow,figh,plotTriggered,'overdraw2','FWR BUA',0.030,0.015);

%% meanadjustment

switch lower(meanadjustment)
    case {'global','local','on','localavoidneighbour'}
        
        switch lower(meanadjustment)
            case {'global','on',}
                fwrbuaadj = getBUA_meanAdjustGlobal1(fwrbua,indRep,befp,aftp,xdata,ydata); % global means
                
            case 'local'
                fwrbuaadj = getBUA_meanAdjustLocal(fwrbua,indRep); % local means considering joined-up spike removal windows
                
            case {'localavoidneighbour'}
                fwrbuaadj = getBUA_meanAdjustLocal_avoidAdjacentSpikes(fwrbua,indRep); % VERY SLOW
         
        end
        
        fwrbua = fwrbuaadj;
        
    case {'off'}
        % nothing
    otherwise
        error('wrong value for meanadjustment')
end

fwrbuanan = local_replacewithnan(wave,filtered,indRep,bua.ChanTitle);

width_ms = local_getwidthOfSpikeRemoval(fwrbuanan);

figh = local_BW(fwrbuanan,spikeEvent,wave,spiketimes,spikewindow,figh,plotBW);

figh = local_plotTriggered(fwrbuanan,wave,spiketimes,showOverdraw,...
    showErrorbar,spikewindow,figh,plotTriggered,'overdraw3','FWR BUA NaN',...
    0.03,0.015);

%% Low pass filtering

cutoffHz = 300;
Fs = fwrbua.SRate;
Wn = cutoffHz/(Fs/2);
[b,a] = butter(n, Wn,'low');

if strcmpi(fvtoolmode,'on')
    fvtool(b,a,'Fs',wave.SRate);
    xlim([0, 0.5]);
end

assert(isstable(b, a))

fwrbuaL = fwrbua;
fwrbuaL.ChanTitle = [fwrbua.ChanTitle,'_Low'];
fwrbuaL.Data = filtfilt(b, a, fwrbua.Data);
if strcmpi(fvtoolmode,'on')
    fvtool(b,a,'Fs',fwrbuaL.SRate);
    xlim([0, 0.5]);
end

figh = local_plotTriggered(fwrbuaL,wave,spiketimes,showOverdraw,...
    showErrorbar,spikewindow,figh,plotTriggered,'overdraw4',...
    'FWR BUA Low-pass filtered',overdrawWidth_sec,overdrawWidth_sec/2);


%% Down Sampling

fwrbuaLD = fwrbuaL.resample(newFs);

spiketimesL = spiketimes; %NOTE rounding problems occur here
if max(spiketimesL) > fwrbuaLD.MaxTime
    spiketimesL(spiketimesL > fwrbuaLD.MaxTime) = [];
    spiketimesL = [spiketimesL;fwrbuaLD.MaxTime];
end
    
figh = local_plotTriggered(fwrbuaLD,wave,spiketimesL,showOverdraw,...
    showErrorbar,spikewindow,figh,plotTriggered,'overdraw5',...
    'FWR BUA Low-pass filtered, Down sampled',overdrawWidth_sec,overdrawWidth_sec/2);


time = wave.time;
repTFdown = timestamps2binned(time(indRep),fwrbuaLD.Start,fwrbuaLD.MaxTime,...
    fwrbuaLD.SRate,'ignore');
indRepD = find(repTFdown); %TODO


end




%% LOCAL FUNCTIONS --------------------------------------------------------

function [wave,from,n,spkthre,spiketimes,spikewindow,fvtoolmode,plotBW,...
    plotFindpeaks,plotTriggered,showOverdraw,showErrorbar,newFs,...
    overdrawWidth_sec,meanadjustment] = parse(wave,from,varargin)

p = inputParser;
assert(isa(wave, 'WaveformChan'));

from = lower(from);
assert(ismember(from,{'wideband','wideband thresholding','highpass',...
    'highpass thresholding'}));

vfeventchan = @(x) isscalar(x) && isa(x,'MetaEventChan') && wave.SRate == x.SRate ...
    && wave.Start == x.Start && wave.Length == x.Length;

switch from
    case 'highpass'
        % [fwrbua,bua,fwrbuanan] = getBUA(wave,'highpass',spiketimes)

        p.addRequired('spiketimes', @(x) (iscolumn(x) && isreal(x) ...
            && min(x) >= wave.StartTime && max(x) <= wave.MaxTime )...
            || vfeventchan(x));
        
    case 'wideband'
        % [fwrbua,bua,fwrbuanan] = getBUA(wave,'wideband',order,spiketimes)
        p.addRequired('order',@(x) (ischar(x) && isrow(x) && strcmpi(x,'')) ...
            || isscalar(x) && fix(x) == x && x > 0);
        
        p.addRequired('spiketimes', @(x) (iscolumn(x) && isreal(x) ...
            && min(x) >= wave.StartTime && max(x) <= wave.MaxTime )...
            || vfeventchan(x));

    case 'wideband thresholding'
        % [fwrbua,bua,fwrbuanan] = getBUA(wave,'wideband thresholding',order,spikeThreshold)
        p.addRequired('order',@(x) (ischar(x) && isrow(x) && strcmpi(x,'')) ...
            || isscalar(x) && fix(x) == x && x > 0);
        
        p.addRequired('spikeThreshold',@(x) isscalar(x) && isreal(x))

     case 'highpass thresholding'
        % [fwrbua,bua,fwrbuanan] = getBUA(wave,'wideband thresholding',order,spikeThreshold)
        p.addRequired('order',@(x) (ischar(x) && isrow(x) && strcmpi(x,'')) ...
            || isscalar(x) && fix(x) == x && x > 0);
        
        p.addRequired('spikeThreshold',@(x) isscalar(x) && isreal(x))       
    otherwise
        error('Value for "from" must be either ''%s'', ''%s'', or ''%s''',...
            'wideband','wideband thresholding','highpass','highpass thresholding')
end

p.addOptional('spikewindow',[0.5,2.5], @(x) isrow(x) && numel(x) == 2 && isreal(x) ...
    && all(x >= 0));

p.addParameter('fvtool','off',@(x) ismember(lower(x),{'on','off'}));

p.addParameter('plotBW','off',@(x) ismember(lower(x),{'on','off'}));

p.addParameter('plotFindpeaks','off',@(x) ischar(x) && isrow(x) &&...
    ismember(lower(x),{'on','off'}));

p.addParameter('plotTriggered','off',@(x) ischar(x) && isrow(x) &&...
    ismember(lower(x),{'on','off'}));

p.addParameter('plotWith','std',@(x) ischar(x) && isrow(x) &&...
    ismember(lower(x),{'overdraw','std','sem','averageonly'}));

p.addParameter('newFs',1024,@(x) isscalar(x) && x > 0);

p.addParameter('overdrawWidth_sec',0.5,@(x) isscalar(x) && x > 0);

p.addParameter('meanadjustment','on',@(x) ismember(x,...
    {'global','local','localavoidneighbour','off','on'}));

p.parse(varargin{:});

switch from
    case 'highpass'
        % [fwrbua,bua,fwrbuanan] = getBUA(wave,'highpass',spiketimes)

        spiketimes = p.Results.spiketimes;
        n = [];
        spkthre = [];
        
    case 'wideband'
        % [fwrbua,bua,fwrbuanan] = getBUA(wave,'wideband',order,spiketimes)
        n = p.Results.order;
        spiketimes = p.Results.spiketimes;
        spkthre = [];

    case 'wideband thresholding'
        % [fwrbua,bua,fwrbuanan] = getBUA(wave,'wideband thresholding',order,spikeThreshold)
        n = p.Results.order;
        spkthre = p.Results.spikeThreshold;
        spiketimes = [];
        
    case 'highpass thresholding'
        % [fwrbua,bua,fwrbuanan] = getBUA(wave,'wideband thresholding',order,spikeThreshold)
        n = p.Results.order;
        spkthre = p.Results.spikeThreshold;
        spiketimes = [];
end

spikewindow       = p.Results.spikewindow;
fvtoolmode        = p.Results.fvtool;
plotBW            = p.Results.plotBW;
plotFindpeaks     = p.Results.plotFindpeaks;
plotTriggered      = p.Results.plotTriggered;
plotWith          = p.Results.plotWith;
newFs             = p.Results.newFs;
overdrawWidth_sec = p.Results.overdrawWidth_sec;
meanadjustment    = p.Results.meanadjustment;

switch plotWith
    case 'overdraw'
        showOverdraw = 'on';
        showErrorbar = 'off';
    case 'std'
        showOverdraw = 'off';
        showErrorbar = 'std';
    case 'sem'
        showOverdraw = 'off';
        showErrorbar = 'sem';
    case 'averageonly'
        showOverdraw = 'off';
        showErrorbar = 'off';
end

end
%--------------------------------------------------------------------------

function fwrbuanan = local_replacewithnan(wave,filtered,indRep,namestem)
% fwrbuanan = local_replacewithnan(wave,filtered,indRep,namestem)
%

BUAnan = filtered; % data only
BUAnan(indRep) = NaN;
FWRBUANaN = abs(BUAnan) - nanmean(BUAnan); % full wave rectification and mean substraction
fwrbuanan = wave;
fwrbuanan.Data = FWRBUANaN;
fwrbuanan.ChanTitle = [namestem, '_FWRBUANaN'];

end

%--------------------------------------------------------------------------

function local_plotBWforNaN(nandata,spikeEvent,Fs,chantitle)

[trigT,~,~,~,seg] = eventTriggeredAverage(single(nandata),spikeEvent,Fs,0.1,0.05);

segnan  = ~isnan(squeeze(seg));
seewps = 1:size(segnan,2);
figure
imagesc(trigT,seewps,segnan')
xlim([-0.010, 0.010])
colormap([1 1 1; 0 0 0])
cbh = colorbar;
cbh.Ticks = [0.25,0.75];
cbh.TickDirection = 'out';
cbh.TickLabels = {'NaN','Data'};
set(gca,'TickDir','out','Box','off');
xlabel('Time relative to spk (msec)');
set(gca,'XTick',-0.010:0.002:0.010,'XTickLabel',-10:2:10);
ylabel('Sweeps');
title(sprintf('%s: event-triggered average\n',chantitle),...
    'Interpreter','none');

end

%--------------------------------------------------------------------------

function width_ms = local_getwidthOfSpikeRemoval(fwrbuanan)
%
% local_getwidthOfSpikeRemoval finds out beginning and end of NaNs and calculate
% the duration of continuous NaNs as actual width of spike removal window,
% which can be joined up together.

loginan = isnan(fwrbuanan.Data);
diffnan = diff(loginan); % note length is -1
ind_nan_start = find(diffnan == 1);
ind_nan_end = find(diffnan == -1);

if ~isempty(ind_nan_start) || ~isempty(ind_nan_end)
    if ind_nan_start(1) > ind_nan_end(1)
        % ignore the first of ind_nan_end
        ind_nan_end(1) = [];
    end
    
    if ind_nan_start(end) > ind_nan_end(end)
        % ignore the last of ind_nan_start
        ind_nan_start(end) = [];
    end
    
    assert(length(ind_nan_start) == length(ind_nan_end),'getBUA:indnan:doesntmatch',...
        'They are supposed to be the same length')
    
    interval_s = fwrbuanan.SInterval;
    width_ms = (interval_s * ind_nan_end - interval_s * ind_nan_start)*1000; 
else
    width_ms = [];
end


end

%--------------------------------------------------------------------------

function local_addVerticalLines(figh,spikewindow)

axh = findobj(figh,'Type','axes');

line(axh,repmat(-spikewindow(1),2,1),ylim(axh),'color',[0.5 0.5 0.5],...
    'tag','spike removal window');

line(axh,repmat(spikewindow(2),2,1),ylim(axh),'color',[0.5 0.5 0.5],...
    'tag','spike removal window');

end

%--------------------------------------------------------------------------

function [figh,xdata,ydata] = local_plotTriggered(waveform,original,spiketimes,showOverdraw,...
    showErrorbar,spikewindow,figh,plotTriggered,field,tag,width,offset)

xdada = [];
ydata = [];

if strcmp(plotTriggered, 'on')
    if ~isempty(spiketimes)
        [data] = waveform.plotTriggered(spiketimes,width,offset,...
            'Overdraw',showOverdraw,'ErrorBar',showErrorbar);
        figh.(field) = gcf;
        figh.(field).Tag = tag;
        title(sprintf('%s: %s',tag,original.ChanTitle));
        local_addVerticalLines(figh.(field),spikewindow)
        
        xdata = data.t;
        ydata = data.mean;
        return
    end
else
    if ~isempty(spiketimes)
        [~,~,~,~,~,data] = waveform.plotTriggered(spiketimes,width,offset,...
            'Average','off','Overdraw','off','ErrorBar','off');
        xdata = data.t;
        ydata = data.mean;
        return
    end
end

end

%--------------------------------------------------------------------------

function figh = local_BW(fwrbuanan,spikeEvent,wave,spiketimes,spikewindow,figh,plotBW)

if strcmpi(plotBW,'on')
    if ~isempty(spiketimes)
        local_plotBWforNaN(fwrbuanan.Data,spikeEvent,fwrbuanan.SRate,fwrbuanan.ChanTitle);
        figh.BW = gcf;
        figh.BW.Tag = 'Black & White';
        title(sprintf('%s, [%.2f,%.2f]',wave.ChanTitle,spikewindow(1),spikewindow(2)))
    end
end

end

%--------------------------------------------------------------------------

function indRep = local_getIndicesOfReplaced(n,spikeEvent,befp,aftp)
% indRep      non-redundant numeric indices for replaced

spkInd = find(spikeEvent);

repStart = spkInd - befp;
repEnd   = spkInd + aftp;

indRep = startsends2ind(repStart,repEnd,n);

end

%--------------------------------------------------------------------------







