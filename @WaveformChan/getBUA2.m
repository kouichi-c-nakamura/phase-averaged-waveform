function [fwrbua,spkfree,fwrbuanan,spiketimes,figh] = getBUA2(wave,from,varargin)
% getBUA2 differs from fromgetBUA in preparation of BUA.
% Median filtering is followed by spike removal and then by high-pass
% filtering.
%
% getBUA returns continuous BUA signal out of wideband LFP signal
%
%   [fwrbua,spkfree,fwrbuanan,spiketimes,figh] = getBUA(wave,'wideband',order,spiketimes)
%
%   [fwrbua,spkfree,fwrbuanan,spiketimes,figh] = getBUA(wave,'wideband thresholding',order,spikeThreshold)
%
%   [fwrbua,spkfree,fwrbuanan,spiketimes,figh] = getBUA(_____, spikewindow)
%   [fwrbua,spkfree,fwrbuanan,spiketimes,figh] = getBUA(_____, ParamName,ParamValue)
%
% INPUT ARGUEMENTS
% wave        A WaveformChan object
%
% from       'wideband'| 'wideband thresholding' | 'highpass'
%             If 'wideband' is chosen, getBUA first high-pass filter the
%             data and then create BUA. 'wideband thresholding' carries out
%             further spike detection with XXXX using findpeaks. If 'highpass' is chosen,
%             the filtering step is skipped.
%
% spiketimes  A column of spike times or a MetaEventChan object with
%             matching Start, SRate and Length
%
% order       A positive integer for the order of butterworth filter.
%             For example, when sampling rate is 17 KHz, order = 9 is good.
%             Use fvtool(b,a,'Fs',samplingrate) for validation
%
% spikeThredhold
%             A zero or positive value for detecting peak locations with
%             findpeaks' 'MinPeakHeight' option.
%
% OPTION
% spikewindow      Two element row vector of zero or positive numbers.
%             [before_ms, after_ms]
%             where before_ms and after_ms are durations relative to a
%             spike time in millisec during which data points are replaced
%             with randomly chosen background activity. See ceateBUA for
%             more details.
%             Default [0.5, 2.5]
%
% OPTIONAL PARAMETER/VALUE PAIRS
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
% 'plotOverdraw'
%             'on' | 'off' (default)
%             Show waveform overdraw for the original waveform data and FWR
%             BUA signal relative to spiketimes
%
% 'plotAverages'
%             'on' | 'off' (default)
%             Only to show overlaid average waveforms. This is overridden
%             by plotOverdraw
%
% OUTPUT ARGUMENTS
% fwrbua      A WaveformChan object holding FWR BUA (full width rectified
%             background unit activity) of the original data.
%             In order to retrieve the envelope of background spiking
%             activity at a low frequency range, you ought to down sample
%             the FWR BUA data or low-pass filter the FWR BUA data.
%             ChanTitle is *_FWRBUA where * is the original ChanTitle.
%
% spkfree     A WaveformChan object holding spike-free data after subtraction 
%             of  moving median filter
%
% fwrbuanan   A WaveformChan object holding FWR BUA of the original data.
%             Instead replacing fragments around spike event with
%             background, they were padded with NaN. They can be useful in
%             order to assess the effect of spike removal.
%
% Algorhythm
% 1. Median filtering (window 0.1 sec) to remove low-frequency fluctuations
% of base line ... medfilt1(x, Fs*0.1)
%
% 2. Spike removal and replacement with background by createBUA() function
%
% 3. High-pass filtering > 300 Hz with butter() and filtfilt()
%
% 4. Rectification and mean subtraction .... FWR BUA signals
%
% 5. Low-pass filtering < 300 Hz with butter() and filtfilt() to retrieve
% envelope of BUA
%
%
% See also
% createBUA (by Izhar Bar-Gad), eventTriggeredAverage, normalizedfreq, timestamps2binned
% getBUA

%TODO overlaid average waveforms

[wave,from,n,plotBW,plotFindpeaks,plotOverdraw,plotAverages,spiketimes,spikewindow,spkthre,...
    fvtoolmode] = local_parse(wave,from,varargin{:});

%% Subtraction of the median filtering of the broadband (LFP) signal,

figh.fvtool = [];
figh.overdraw_original = [];
figh.overdraw_spikeremoved = [];
figh.overdraw_highpassfiltered = [];
figh.overdraw_FWRBUA = [];
figh.overdraw_FWRBUANaN = [];
figh.BW = [];


medsubtructed_data = wave.Data - medfilt1(wave.Data,wave.SRate*0.02); % remove 0-10 Hz oscillations

if strcmp(from,'wideband thresholding')
    
    wh = which('findpeaks');
    if ~isempty(wh) && isempty(strfind(wh, fullfile('toolbox','signal','signal','findpeaks.m')))
        rmpath(fileparts(wh));
    end
    
    warning off 'signal:findpeaks:largeMinPeakHeight'
    [~,locs] = findpeaks(medsubtructed_data,'MinPeakHeight',spkthre);
    t = wave.time;
    warning on 'signal:findpeaks:largeMinPeakHeight'
    
    spiketimes = t(locs);
    
    addpath(fileparts(wh));
end

if strcmp(plotFindpeaks,'on')
    if ~isempty(spiketimes)
        findpeaks(medsubtructed_data,wave.SRate,'MinPeakHeight',spkthre);
        set(gca,'TickDir','out','Box','off')
        xlabel('Time (sec)')
        ylabel('Amplitude (mV)')
        zoom xon, pan xon;
        title(['Magnituide Response (dB):',wave.ChanTitle])

        figh.fvtool = gcf;
        figh.fvtool.Tag = 'findpeaks';
    end

end

if strcmp(plotOverdraw, 'on') || strcmp(plotAverages, 'on')
    if ~isempty(spiketimes)
        if strcmp(plotOverdraw, 'on') 
            wave.plotTriggered(spiketimes,0.020,0.010);
        else
            wave.plotTriggered(spiketimes,0.020,0.010,...
                'Overdraw','off');
        end
        figh.overdraw_original = gcf;
        figh.overdraw_original.Tag = 'original';

        hold on
        if strcmp(from,'wideband thresholding')
            plot(xlim,[spkthre,spkthre],'k--','Tag','spkthre');
        end
        delete(findobj(gcf,'Type','legend'));
        
        if strcmp(plotOverdraw, 'on')
            
            local_plotspikewindow(gca,spikewindow);
            warning off 'MATLAB:legend:IgnoringExtraEntries'
            leg = legend([findobj(gcf,'Tag','Average Waveform'),...
                findobj(gcf,'Tag','Dummy'),...
                findobj(gcf,'Tag','spkthre')],...
                {'Average','Overdraw','Thershold'},...
                'Box','off');
            warning on 'MATLAB:legend:IgnoringExtraEntries'
            
            title(['Original: ',wave.ChanTitle])
        end

        medsub = wave;
        medsub.ChanTitle = [wave.ChanTitle,'_mediansubtracted'];
        medsub.Data = medsubtructed_data;
        if strcmp(plotOverdraw, 'on') 
            medsub.plotTriggered(spiketimes,0.020,0.010);
        else
            medsub.plotTriggered(spiketimes,0.020,0.010,...
                'Overdraw','off');
        end
        figh.overdraw_medsub = gcf;
        figh.overdraw_medsub.Tag = 'median subtracted';
        title(['Median subtracted: ',wave.ChanTitle])

    end
end


%% Replace the high amplitude units with data from elsewhere in the signal (see Moran et al, 2006)

if isa(spiketimes,'MetaEventChan')
    spkevent = spiketimes.Data;
    spiketimes = spiketimes.TimeStamps;
else
    [spkevent, spiketimes] = timestamps2binned(spiketimes,wave.Start,wave.MaxTime,wave.SRate);
end

if isempty(spiketimes)
    spkfree_data = medsubtructed_data;
    if strcmp(plotOverdraw, 'on') || strcmpi(plotBW,'on') || strcmp(plotFindpeaks,'on')...
            || strcmp(plotAverages, 'on')
        fprintf('No spike over %f has been detected\n',spkthre)
    end
else
    spkfree_data = createBUA(medsubtructed_data, spiketimes, wave.SRate, 0, spikewindow)';
end

spkfree = wave;
spkfree.Data = spkfree_data;
spkfree.ChanTitle = [spkfree.ChanTitle, '_spikefree'];

if strcmp(plotOverdraw, 'on') || strcmp(plotAverages, 'on')
    if ~isempty(spiketimes)
        if strcmp(plotOverdraw, 'on')
            spkfree.plotTriggered(spiketimes,0.020,0.010);
        else
            spkfree.plotTriggered(spiketimes,0.020,0.010,...
                'Overdraw','off');
        end
        local_plotspikewindow(gca,spikewindow);

        figh.overdraw_spikeremoved = gcf;
        figh.overdraw_spikeremoved.Tag = 'spike removed';
        title(['spike removed: ',wave.ChanTitle])
    end
end


%% High-pass filtering (300 Hz)

Wn = normalizedfreq(300,spkfree.SRate);
[b,a] = butter(n, Wn,'high'); % the dimension n can be empirically explored
assert(isstable(b,a));

if strcmpi(fvtoolmode,'on')
    fvtool(b,a,'Fs',spkfree.SRate);
    xlim([0, 0.5]);
end
hpfiltered = filtfilt(b, a, spkfree.Data);

if strcmp(plotOverdraw, 'on') || strcmp(plotAverages, 'on')
    if ~isempty(spiketimes)
        hpf = wave;
        hpf.Data = hpfiltered;
        hpf.ChanTitle = [wave.ChanTitle, '_highpass'];
        if strcmp(plotOverdraw, 'on') 
            hpf.plotTriggered(spiketimes,0.020,0.010);
        else
            hpf.plotTriggered(spiketimes,0.020,0.010,...
                'Overdraw','off');
        end
        local_plotspikewindow(gca,spikewindow);

        figh.overdraw_highpassfiltered = gcf;
        figh.overdraw_highpassfiltered.Tag = 'High-pass filtered';
        title(['High-pass filtered: ',wave.ChanTitle])
    end
end

%% Rectify and mean subtraction

FWRBUA = abs(hpfiltered) - mean(hpfiltered);

fwrbua = wave;
fwrbua.Data = FWRBUA;
fwrbua.ChanTitle = [fwrbua.ChanTitle, '_FWRBUA'];

if strcmp(plotOverdraw, 'on') || strcmp(plotAverages, 'on')
    if ~isempty(spiketimes)
        if strcmp(plotOverdraw, 'on') 
            fwrbua.plotTriggered(spiketimes,0.020,0.010);
        else
            fwrbua.plotTriggered(spiketimes,0.020,0.010,...
                'Overdraw','off');
        end
        local_plotspikewindow(gca,spikewindow);

        figh.overdraw_FWRBUA = gcf;
        figh.overdraw_FWRBUA.Tag = 'FWR BUA';
        title(sprintf('%s, [%.2f,%.2f]',wave.ChanTitle,spikewindow(1),spikewindow(2)))
    end
    
end

FWRBUANaN = replacewithnan(hpfiltered,spkevent,wave.SRate,spikewindow);
fwrbuanan = wave;
fwrbuanan.Data = FWRBUANaN;
fwrbuanan.ChanTitle = [wave.ChanTitle, '_FWRBUANaN'];

if strcmpi(plotBW,'on')
    if ~isempty(spiketimes)
        plotBWforNaN(fwrbuanan.Data,spkevent,fwrbuanan.SRate,fwrbuanan.ChanTitle);
        figh.BW = gcf;
        figh.BW.Tag = 'Black & White';
        title(sprintf('%s, [%.2f,%.2f]',wave.ChanTitle,spikewindow(1),spikewindow(2)))
    end
    
end

if strcmp(plotOverdraw, 'on') || strcmp(plotAverages, 'on')
    if ~isempty(spiketimes)
        if strcmp(plotOverdraw, 'on')
            fwrbuanan.plotTriggered(spiketimes,0.020,0.010);
        else
            fwrbuanan.plotTriggered(spiketimes,0.020,0.010,...
                'Overdraw','off');
        end
        local_plotspikewindow(gca,spikewindow);

        figh.overdraw_FWRBUANaN = gcf;
        figh.overdraw_FWRBUANaN.Tag = 'FWR BUA NaN';
        title(sprintf('%s, [%.2f,%.2f]',wave.ChanTitle,spikewindow(1),spikewindow(2)))
    end
end
   
if strcmp(plotOverdraw, 'on') || strcmp(plotAverages, 'on')
    if ~isempty(spiketimes)

        figh.overdraw_all = figure;
        newh = copyobj(findobj([figh.overdraw_original,...
            figh.overdraw_medsub,...
            figh.overdraw_spikeremoved,...
            figh.overdraw_highpassfiltered,...
            figh.overdraw_FWRBUA,...
            figh.overdraw_FWRBUANaN],...
            'Tag','Average Waveform'),axes);
        
        for i = 1:length(newh)
            newh(i).Color = defaultPlotColors(i);
        end
        local_plotspikewindow(gca,spikewindow);

        legend(newh,{'Original',...
            'Median filtered',...
            'Spike removed',...
            'High-pass filtered',...
            'FWR BUA',...
            'FWR BUA NaN'},...
            'Box','off')
        ylabel('Amplitude (mV)')
        xlabel('Time relative to event (msec)')
        title('Overlaid','Interpreter','none')
        set(gca,'Box','off','TickDir','out')
        clear i newh
        
        if strcmp(plotAverages, 'on') && ~strcmp(plotOverdraw, 'on')
            close([figh.overdraw_original,...
                figh.overdraw_medsub,...
                figh.overdraw_spikeremoved,...
                figh.overdraw_highpassfiltered,...
                figh.overdraw_FWRBUA,...
                figh.overdraw_FWRBUANaN]);
        end
    end
end

end

%--------------------------------------------------------------------------

function [FWRBUAnan] = replacewithnan(WBHF,spkdata,Fs,spikeLen)

spikeTrain = find(spkdata);
N = length(spkdata);

befp = round(Fs/1000*spikeLen(1));
aftp = round(Fs/1000*spikeLen(2));

left = spikeTrain - befp;
right = spikeTrain + aftp;

ind = arrayfun(@(x,y) x:y, left, right, 'UniformOutput',false);
ind2 = cellfun(@(x) getridofexcess(x,N), ind,'UniformOutput',false);

index =[ind2{:}];

BUAnan = WBHF;
BUAnan(index) = NaN;

FWRBUAnan = abs(BUAnan) - nanmean(BUAnan);

end

%--------------------------------------------------------------------------

function plotBWforNaN(nandata,spkdata,Fs,chantitle)

[trigT,~,~,~,seg] = eventTriggeredAverage(single(nandata),spkdata,Fs,0.1,0.05);

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

function X = getridofexcess(x,n)

x(x < 1) = [];
x(x > n) = [];

X = x;

end

%--------------------------------------------------------------------------

function hptch = local_plotspikewindow(axh,spikewindow)

currentYLim = ylim(axh);

axes(axh);
hold on
hptch = patch([-spikewindow(1),-spikewindow(1),spikewindow(2),spikewindow(2)],...
    [-5,5,5,-5],'c','FaceAlpha',0.1,'EdgeColor','none','Tag','spike window');

axh.Children = axh.Children([2,1,3:end]);
ylim(axh,currentYLim);

end

%--------------------------------------------------------------------------

function [wave,from,n,plotBW,plotFindpeaks,plotOverdraw,plotAverages,spiketimes,...
    spikewindow,spkthre,fvtoolmode] = local_parse(wave,from,varargin)

p = inputParser;
assert(isa(wave, 'WaveformChan'));

from = lower(from);
assert(ismember(from,{'wideband','wideband thresholding'}));

vfeventchan = @(x) isscalar(x) && isa(x,'MetaEventChan') && wave.SRate == x.SRate ...
    && wave.Start == x.Start && wave.Length == x.Length;

switch from
    %     case 'highpass'
    %         % [fwrbua,spkfree,fwrbuanan] = getBUA(wave,'highpass',spiketimes)
    %
    %         p.addRequired('spiketimes', @(x) (iscolumn(x) && isreal(x) ...
    %             && min(x) >= wave.StartTime && max(x) <= wave.MaxTime )...
    %             || vfeventchan(x));
    
    case 'wideband'
        % [fwrbua,spkfree,fwrbuanan] = getBUA(wave,'wideband',order,spiketimes)
        p.addRequired('order',@(x) (ischar(x) && isrow(x) && strcmpi(x,'')) ...
            || isscalar(x) && fix(x) == x && x > 0);
        
        p.addRequired('spiketimes', @(x) (iscolumn(x) && isreal(x) ...
            && min(x) >= wave.Start && max(x) <= wave.MaxTime )...
            || vfeventchan(x));
        
    case 'wideband thresholding'
        % [fwrbua,spkfree,fwrbuanan] = getBUA(wave,'wideband thresholding',order,spikeThreshold)
        p.addRequired('order',@(x) (ischar(x) && isrow(x) && strcmpi(x,'')) ...
            || isscalar(x) && fix(x) == x && x > 0);
        
        p.addRequired('spikeThreshold',@(x) isscalar(x) && isreal(x) && x >= 0)
        
        
    otherwise
        error('Value for "from" must be either ''%s'', ''%s'', or ''%s''',...
            'wideband','wideband thresholding')
end

p.addOptional('spikewindow',[0.5,2.5], @(x) isrow(x) && numel(x) == 2 && isreal(x) ...
    && all(x >= 0));

p.addParameter('fvtool','off',@(x) ismember(lower(x),{'on','off'}));

p.addParameter('plotBW','off',@(x) ismember(lower(x),{'on','off'}));

p.addParameter('plotFindpeaks','off',@(x) ischar(x) && isrow(x) &&...
    ismember(lower(x),{'on','off'}));

p.addParameter('plotOverdraw','off',@(x) ischar(x) && isrow(x) &&...
    ismember(lower(x),{'on','off'}));

p.addParameter('plotAverages','off',@(x) ischar(x) && isrow(x) &&...
    ismember(lower(x),{'on','off'}));

p.parse(varargin{:});

switch from
    %     case 'highpass'
    %         % [fwrbua,spkfree,fwrbuanan] = getBUA(wave,'highpass',spiketimes)
    %
    %         spiketimes = p.Results.spiketimes;
    %         n = [];
    %         spkthre = [];
    
    case 'wideband'
        % [fwrbua,spkfree,fwrbuanan] = getBUA(wave,'wideband',order,spiketimes)
        n = p.Results.order;
        spiketimes = p.Results.spiketimes;
        spkthre = [];
        
    case 'wideband thresholding'
        % [fwrbua,spkfree,fwrbuanan] = getBUA(wave,'wideband thresholding',order,spikeThreshold)
        n = p.Results.order;
        spkthre = p.Results.spikeThreshold;
        spiketimes = [];
end

spikewindow = p.Results.spikewindow;
fvtoolmode = p.Results.fvtool;
plotBW = p.Results.plotBW;
plotFindpeaks = p.Results.plotFindpeaks;
plotOverdraw = p.Results.plotOverdraw;
plotAverages = p.Results.plotAverages;

end
