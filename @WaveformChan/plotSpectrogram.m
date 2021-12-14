function [out] = plotSpectrogram(obj, newRate, window, nfft, varargin)
%[out] = plotSpectrogram(obj, newRate, window, nfft, varargin)
%
% INPUT ARGUMENTS
% obj         A WaveformChan object
%
% newRate     [Hz]. obj.Data is to be resamled to this sampling frequency
%
% window      The size of Hamming window. A positive even integer. The size of window in
%             points after resampling. The 50% of window will
%             be overlapped (noverlap = window/2).
%
% nfft        the number of fft. A positive interger, must be
%             power of 2. Minimum is 256. Frequency resolution
%             of the resultant power spectrum is newRate/nfft
%
% OPTIONAL PARAMETER/VALUE PAIRS
% 'PlotType'  'normal'  (default) | '1/f'
%             With '1/f' option, the power in spectrogram is devided by 1/f
%
% 'PlotFun'   'imagesc' (default) | 'pcolor'
%             'pcolor' uses pcolor with shading interp
%
% 'Fspecial'  {} (default) | {'gaussian', hsize, sigma}  
%             Cell array format. By default, or with {}, no filter will
%             beused. With {'gaussian', hsize, sigma} option, you can apply
%             a rotationally symmetric Gaussian lowpass filter of size
%             hsize with standard deviation sigma (positive). hsize can be
%             a vector specifying the number of rows and columns in h, or
%             it can be a scalar, in which case h is a square matrix. The
%             default value for hsize is [3 3]; the default value for sigma
%             is 0.5. Type 'help fspecial' for more information.
%
% 'CLim'         [left, right] in Hz. Default is [0 100]
%
% 'YLim'         [bottom, top]. Default is [0 1e-3]
%
% 'DoPlot'       true (default) | false 
%
% OUTPUT ARGUMENTS
% out         structure with the following fields
% 
%                 'S'
%                 'F'
%                 'T'
%                 'P'
%                 'P_dB'       power in dB
%                 'PoverF_dB'  power / f in dB
%                 'P_final'    depends on 'PlotType' and 'FSpecial' options
%                 'window'
%                 'noverlap'
%                 'nfft'
%                 'newRate'
%                 'frequencyRange'
%                 'bins'
%                 'frequencyResolution'
%                 'fftdur'
%                 'clim'
%                 'colorbar'
%                 'h1'
%                 'ax1'
%                 'fig1'
% 
%
% See also
% spectrogram
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% kouichi.c.nakamura@gmail.com
% 05-Nov-2016 05:10:03


%% parse inputs

[obj,newRate,window,nfft,plot1overFcorrection,crange,yrange,dopcolor,...
    dofilter,doplot,filtparam] = local_parse(obj,newRate,window,nfft,varargin{:});

%% job
x = resample(obj.Data, newRate, round(obj.SRate));
target_title = obj.ChanTitle;

noverlap = window/2;
[S, F, T, P] = spectrogram(x,window,noverlap,nfft,newRate);
P_dB = 10*log10(abs(P)); % 10*log10(abs(P)) see  http://www.mathworksco.uk/help/signal/ref/spectrogram.html
PoverF = P.*repmat(F,1,size(P,2));
PoverF_dB =10*log10(abs(PoverF)); % in [dB]

switch plot1overFcorrection
    case false
        P_final = P_dB;
    case true
        P_final = PoverF_dB;
end

% igonore if the first row is inf (1/f mode)
if any(any(isinf(P_final)))
    % hide inf or -inf from plot
    if all(isinf(P_final(1,:))) && ...
            ~any(any(isinf(P_final(2:end,:))))
        % only the first row (F == 0) need to be ignored
        F = F(2:end);
        P_final = P_final(2:end,:);
    else
        warning('K:plotSpectrogram:inf', ...
            'infinite values found in P_final in other than the first row.');
    end
end

% image filtering
if dofilter
    imagefilter = fspecial(filtparam{:});
    P_final = imfilter(P_final, imagefilter, 'replicate'); %TODO any shift?
end

T = T + obj.Start; % support when obj.Start ~= 0

if doplot
    fig1 = figure;
    
    switch dopcolor
        case true
            h1 = pcolor(T,F,P_final);
            shading(gca, 'interp');
        case false
            h1 = imagesc(T,F,P_final);
    end
    
    set(gca, 'CLim', crange);
    set(gca, 'YDir', 'normal', 'TickDir', 'out', 'Box', 'off');
    title(sprintf('Spectrogram: %s', target_title));
    ylim(yrange);
    xlabel('Time (sec)');
    ylabel('Frequency (Hz)');
    ax1 = gca;
    colb = colorbar;
    % ylabel(colorbar, '[dB]'); % need to use text() instead
    % zoom xon; pan xon;
    % set(colorbar, 'TickDir', 'out',  'Box', 'off'); % doesn't work properly
    
else
    h1 = [];
    ax1 =[];
    fig1 = [];
    colb = [];
end


out.S = S;
out.F = F;
out.T = T;
out.P = P;
out.P_dB = P_dB;
out.PoverF_dB = PoverF_dB;
out.P_final = P_final;

out.window =window;
out.noverlap = noverlap;
out.nfft = nfft;
out.newRate = newRate;
out.frequencyRange = [0, newRate/2];
out.bins = []; %TODO
out.frequencyResolution = newRate/nfft; % [Hz]
out.fftdur = nfft*1/newRate; % [sec]
out.clim = crange;
out.colorbar = colb;

out.h1 = h1; 
out.ax1 = ax1;
out.fig1 = fig1;

end

%--------------------------------------------------------------------------

function [obj,newRate,window,nfft,plot1overFcorrection,crange,yrange,...
    dopcolor,dofilter,doplot,filtparam] = local_parse(obj,newRate,window,...
    nfft,varargin)


if nargin == 1
    if obj.SRate > 1024
        newRate = 1024;
        window = 512;
        nfft = 512;
    else
        newRate = obj.SRate;
        window = 256;
        nfft = 256;
    end
    
else
    narginchk(4, inf);
end

p = inputParser;

vf_newRate = @(x) isnumeric(x) &&...
    isscalar(x) &&...
    (x > 0) && ...
    (x <= obj.SRate);

p.addRequired('newRate', vf_newRate);

vf_window = @(x) isnumeric(x) &&...
    fix(x) == x &&...
    x > 0 &&...
    rem(x, 2) == 0; % positive even integer

p.addRequired('window', vf_window);

vf_nfft = @(x) isnumeric(x) &&...
    fix(x) == x &&...
    x >= 256 &&...
    rem(x, 2) == 0 && ... % positive even integer
    fix(log2(x)) == log2(x) ; % power of 2

p.addRequired('nfft', vf_nfft);

parse(p, newRate, window, nfft);


if window * 10 > obj.Length
    warning('K:WaveformChan:plotSpectrogram:window',...
        'Data length is shorter than 10*window. The FFT result may not be reliable.');
end



%% optional input agruments
plot1overFcorrection = false;
crange = [-80 -30];
yrange = [0 100];
dopcolor = false;
dofilter = false;
doplot = true;
filtparam = [];

ni = length(varargin); % ni >= 1
PNVStart = 1;
for i=PNVStart:2:ni
    % Set each Property Name/Value pair in turn.
    Property = varargin{i};
    if i+1>ni
        error('K:WaveformChan:plotSpectrogram:options:pvsetNoValue', 'Value is missing')
    else
        Value = varargin{i+1};
    end
    
    % Perform assignment
    switch lower(Property)
        case 'plottype'
            if ~isempty(Value) && ischar(Value) && isrow(Value)
                Value = validatestring(Value, {'normal', '1/f'});
                if strcmpi(Value, 'normal')
                    plot1overFcorrection = false;
                elseif strcmpi(Value, '1/f')
                    plot1overFcorrection = true;
                else
                    error('K:WaveformChan:plotSpectrogram:options:PlotType:invalid', 'PlotType value invalid')
                end
            else
                error('K:WaveformChan:plotSpectrogram:options:PlotType:invalid', 'PlotType value invalid')
            end
        case 'plotfun'
            if ~isempty(Value) && ischar(Value) && isrow(Value)
                Value = validatestring(Value, {'imagesc', 'pcolor'});
                if strcmpi(Value, 'imagesc')
                    dopcolor = false;
                elseif strcmpi(Value, 'pcolor')
                    dopcolor = true;
                else
                    error('K:WaveformChan:plotSpectrogram:options:PlotFun:invalid', 'plotfun value invalid')
                end
            else
                error('K:WaveformChan:plotSpectrogram:options:PlotFun:invalid', 'plotfun value invalid')
            end
        case 'fspecial'
            if ~isempty(Value) && iscell(Value) && isvector(Value) && ischar(Value{1})
                dofilter = true;
                filtparam = Value;
                
            elseif ~isempty(Value) && iscell(Value) && isvector(Value) && ~ischar(Value{1})
                error('K:WaveformChan:plotSpectrogram:options:Fspecial:invalid', 'fspecial value invalid')
            elseif isempty(Value)
                dofilter = false;
            else
                error('K:WaveformChan:plotSpectrogram:options:Fspecial:invalid', 'fspecial value invalid')
            end
            
        case 'clim'
            if ~isempty(Value) && isnumeric(Value) && isrow(Value) && (numel(Value) == 2) && ...
                    Value(1) <= Value(2) 
                crange = Value;
            else
                error('K:WaveformChan:plotSpectrogram:options:CLim:invalid', 'CLim value invalid')
            end
        case 'ylim'
            if ~isempty(Value) && isnumeric(Value) && isrow(Value) && (numel(Value) == 2) && ...
                    Value(1) >= 0 && Value(2) > 0
                yrange = Value;
            else
                error('K:WaveformChan:plotSpectrogram:options:YLim:invalid', 'YLim value invalid')
            end
        case 'doplot'
            if isscalar(Value) && islogical(Value) || isnumeric(Value) &&  ...
                    (Value == 0) || (Value == 1)
                doplot = Value;
            else
                error('K:WaveformChan:plotSpectrogram:options:doplot:invalid', 'doplot value invalid')
            end
        otherwise
            error('K:WaveformChan:plotSpectrogram:option:pvsetInvalid', 'Parameter and values are invalid')
    end % switch
end % for

if plot1overFcorrection && ~any(strcmpi('clim', varargin))
    crange = [-45 -30];
end


end

