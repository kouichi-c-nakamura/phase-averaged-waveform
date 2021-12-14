function [h, out] = plotPowerSpectrum(obj, varargin)
%[h, out] = plotPowerSpectrum(obj)
%[h, out] = plotPowerSpectrum(obj, newRate, window, nfft)
%[h, out] = plotPowerSpectrum(obj, ______, 'Param', value, ...)
%[h, out] = plotPowerSpectrum(obj, axh, ______) %TODO
%
% INPUT ARGUMENTS
%
% newRate      [Hz]. obj.Data is to be resamled to this sampling frequency
%         	   (default) 1024 (obj.SRate > 1024) or 256
%
% window       The size of Hamming window. A positive even integer. The size
%              of window in points after resampling. 
%
% nfft         the number of fft. A positive interger, must be
%              power of 2. Minimum is 256. Frequency resolution
%              of the resultant power spectrum is newRate/nfft.
%
% OPTIONAL PARAMETER/VALUE PAIRS
%
% 'PlotType'     'line'   line drawing
%                'hist'   histogram
%                'none'   no plot produced and only returns out and empty h.
%
% 'XLim'         [left, right] in Hz. Default is [0 100]
%
% 'YLim'         [bottom, top]. Default is [0 1e-3]
%
% 'noverlap'     The number of overlap.
%                (default) window/2
%
% 'Normalize'    [] (default) | [lowHz, highHz]#
%                If empty, absolute power will be plotted. When specified,
%                noramlized power with normalizepower function will be
%                shown.
%
% OUTPUT ARGUMENTS
%
% h            Graphic handles in struct with fields
%                     h.figure
%                     h.axes
%                     h.line
%                     h.text
%
% out          A strucut with power spectram data
%                     out.Pxx      Power
%                     out.f        Frequency vector
%                     out.window   window for pwelch
%                     out.noverlap noverlap for pwelch
%                     out.nfft     nfft for pwelch
%                     out.newRate  newRate for pwelch
%                     out.frequencyRange         = [0, newRate/2]
%                     out.frequencyResolution    = newRate/nfft
%
% See also
% pwelch, resample, WaveformChan.plotCohere,
% textforpwelch


%% parse inputs

[obj, newRate, window, noverlap, nfft, xrange, yrange, plotType, axh, normalizeLowHigh] = ...
    local_parse(obj,varargin{:});

%% job
x = resample(obj.Data, newRate, obj.SRate);
target_title = obj.ChanTitle;

[Pxx,f] = pwelch(x, window, noverlap, nfft, newRate);

if ~isempty(normalizeLowHigh)
    Pxxn = normalizepower(Pxx,f,normalizeLowHigh(1),normalizeLowHigh(2))*100;
else
    Pxxn = [];
end
%% plot

switch plotType
    case {'line', 'hist'}
        if isempty(axh)
            fig1 = figure;
            axh = axes;
        else
            axes(axh);
       
            fig1 = local_recursiveparentfigure(axh);
        end
        
        switch plotType
            case 'line'
                 if isempty(normalizeLowHigh)
                     h1 = line(axh, f, Pxx,'color',defaultPlotColors(1));

                 else
                     h1 = line(axh, f, Pxxn,'color',defaultPlotColors(1));

                 end
            case 'hist'
                if isempty(normalizeLowHigh)
                    h1 = bar(axh, f, Pxx, 'BarWidth', 1);

                else
                    h1 = bar(axh, f, Pxxn, 'BarWidth', 1);

                end
        end
        
        xlim(axh,xrange);
        ylim(axh,yrange);
        xlabel(axh,'Frequency (Hz)');
        if normalizeLowHigh
            ylabel(axh,'Relative Power (%)');
        else
            ylabel(axh,'Power (mV^2)');
        end
        set(axh, 'TickDir', 'out', 'Box', 'off');
        title(axh,sprintf('Power spectrum: %s', target_title),'Interpreter','none');
        
        %% customize YTickLabel
%         ytick = get(axh, 'YTick');
%         ytickc = cell(1, length(ytick));
%         decimalp = -floor(log10(max(yrange)))+1;
%         for i = 1:length(ytick)
%             if ytick(i) == 0
%                 ytickc{i} = '0';
%             else
%                 ytickc{i} = num2str(ytick(i), ['%.', num2str(round(decimalp)),'f']);
%             end
%         end
%         set(axh, 'YTick', ytick, 'YTickLabel', ytickc);
        set(axh,'TickDir', 'out', 'Box', 'off');
        texth = textforpwelch(axh, window, noverlap, nfft, newRate);
        
    case 'none'
        h1 = [];
        fig1 = [];
        axh = [];
        texth = [];
end


%% output

out.Pxx = Pxx;
out.f = f;

out.window =window;
out.noverlap = noverlap;
out.nfft = nfft;
out.newRate = newRate;
out.frequencyRange = [0, newRate/2];
out.frequencyResolution = newRate/nfft; % [Hz]



h.figure = fig1;
h.axes = axh;
h.line = h1;
h.text = texth;

end

%--------------------------------------------------------------------------

function [obj, newRate, window, noverlap, nfft, xrange, yrange, plotType, ...
    axh, normalizeLowHigh] = local_parse(obj, varargin)
% [obj, newRate, window, noverlap, nfft, xrange, yrange, plotType, axh] ...
%     = local_parse(obj, varargin)
%
% See also
% plotPowerSpectra_parse_test.m

p = inputParser;
p.addRequired('obj');

if ~isempty(varargin) && isscalar(varargin{1}) && ishandle(varargin{1})
    
    axh = varargin{1};
    assert(isequal(axh.Type,'axes'));
    varargin = varargin(2:end);
    
else
    axh = [];
end

if ~isempty(varargin) && isnumeric(varargin{1})
    assert(~ischar(varargin{2}) &&  ~ischar(varargin{3}),...
        'plotPowerSpectra:windownfft:missing',...
        'The syntax plotPowerSpectrum(obj, newRate, window, nfft) requires all  newRate, window and nfft at the same time');
end

vfscpos = @(x) ~isempty(x) && isscalar(x) && isnumeric(x) && x > 0;
p.addOptional('newRate', [], vfscpos);

vfscposint = @(x) ~isempty(x) && isscalar(x) && isnumeric(x) && ...
    x > 0 && fix(x) == x;

p.addOptional('window', [], vfscposint);
p.addOptional('nfft',   [], vfscposint);

p.addParameter('plotType', 'line', ...
    @(x) ~isempty(x) && ismember(lower(x),{'line','hist','none'}));

vflim = @(x) ~isempty(x) && isreal(x) && isrow(x) && ...
    (numel(x) == 2) && x(1) >= 0 && x(2) > 0 && x(1) < x(2);

p.addParameter('xlim', [0 100], vflim);

p.addParameter('ylim', [0 1e-4], vflim);

p.addParameter('noverlap', [], vfscposint);

p.addParameter('Normalize', [], @(x) numel(x) ==2 && isrow(x) && isreal(x));


p.parse(obj,varargin{:});

newRate = p.Results.newRate;
window = p.Results.window;
nfft = p.Results.nfft;
plotType = p.Results.plotType;
xrange = p.Results.xlim;
yrange = p.Results.ylim;
noverlap = p.Results.noverlap;
normalizeLowHigh = p.Results.Normalize;

if isempty(window)
    if obj.SRate >= 1024
        newRate = 1024;
        window = 1024;
        nfft = 1024;
    else
        newRate = obj.SRate;
        window = 256;
        nfft = 256;
    end
end

if isempty(noverlap)
    noverlap = round(window/2);
end

if ~isempty(normalizeLowHigh)
    
    if  isequal(yrange,[0 1e-4])
        
        nbin = diff(xrange)/(newRate/nfft);
        yrange = [0 100/nbin*5]; % default for normailized
    end
    
end

end
%--------------------------------------------------------------------------

function out = local_recursiveparentfigure(obj)

if isgraphics(obj.Parent,'figure')
    out = obj.Parent;
else
    out = local_recursiveparentfigure(obj.Parent);
end
  
end
