function [h, out] = plotCohere(obj, wf2, varargin)
% plotCohere is a wrapper of mscoher
% [h, out] = plotCohere(obj, wf2)
% [h, out] = plotCohere(obj, wf2, newRate, window, nfft)
% [h, out] = plotCohere(obj, axh,_____)
% [h, out] = plotCohere(obj, _____, 'Param',value,...)
%

% INPUT ARGUMENTS
%
% newRate      [Hz]. obj.Data is to be resamled to this sampling frequency
%         	   (default) 1024 (obj.SRate > 1024) or 256
%
% window       The size of Hamming window. A positive even integer. The size of window in
%              points after resampling. The 50% of window will
%              be overlapped.
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
%                     out.Cxy      Coherence
%                     out.f        Frequency vector
%                     out.window   window for pwelch
%                     out.noverlap noverlap for pwelch
%                     out.nfft     nfft for pwelch
%                     out.newRate  newRate for pwelch
%                     out.frequencyRange         = [0, newRate/2]
%                     out.frequencyResolution    = newRate/nfft
%
% See also
% pwelch, resample, WaveformChan.plotPowerSpectrum,
% plotCohere, textforpwelch


[obj, wf2, newRate, window, noverlap, nfft, xrange, yrange, plotType, axh] ...
    = local_parse(obj, wf2, varargin{:});

%% job

x = resample(obj.Data, newRate, obj.SRate);
y = resample(wf2.Data, newRate, wf2.SRate);

first_title = obj.ChanTitle;
second_title = wf2.ChanTitle;

[Cxy,f] = mscohere(x, y, window, noverlap, nfft, newRate);


%% plot

switch plotType
    case {'line','hist'}
        if isempty(axh)
            fig1 = figure;
            axh = axes;
        else
            axes(axh);
            fig1 = gcf;
        end
        
        switch plotType
            case 'line'
                h1 = plot(axh, f, Cxy);
            case 'hist'
                h1 = bar(axh, f, Cxy, 'BarWidth', 1);
        end
        
        xlim(xrange);
        ylim(yrange);
        xlabel('Frequency (Hz)');
        ylabel('Coherence');
        set(axh, 'TickDir', 'out', 'Box', 'off');
        title(sprintf('Coherence spectrum: %s VS %s', first_title, second_title),...
            'Interpreter','none');
        
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
%         set(axh, 'YTick', ytick, 'YTickLabel', ytickc)
        set(axh,'TickDir', 'out', 'Box', 'off');
        texth = textforpwelch(gca, window, noverlap, nfft, newRate);

    case 'none'
        fig1 = [];
        axh = [];
        h1 = [];
        texth = [];
end


%% output

out.Cxy = Cxy;
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

function [obj, wf2, newRate, window, noverlap, nfft, xrange, yrange, plotType, axh] ...
    = local_parse(obj, wf2, varargin)
% [obj, wf2, newRate, window, noverlap, nfft, xrange, yrange, plotType, axh] ...
%     = local_parse(obj, wf2, varargin)
%
% See also
% plotPowerSpectra_parse_test.m

p = inputParser;
p.addRequired('obj');
p.addRequired('wf2');


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

p.addParameter('ylim', [0 1], vflim);

p.addParameter('noverlap', [], vfscposint);

p.parse(obj,wf2,varargin{:});


newRate = p.Results.newRate;
window = p.Results.window;
nfft = p.Results.nfft;
plotType = p.Results.plotType;
xrange = p.Results.xlim;
yrange = p.Results.ylim;
noverlap = p.Results.noverlap;


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

end