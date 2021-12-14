function [val, out] = powerRatio(obj, varargin)
% [val, out] = powerRatio(obj, frange, wholefrange , newRate, window, nfft)
% [val, out] = powerRatio(obj, _____, 'Preset', 'band')
%
%
% INPUT ARGUMENTS
% A           integer | char
%             Description about A comes here.
%
%
% B           0 (default) | non-negative integers
%             (Optional) Description about B comes here.
%
% frange      [low, high] frequnency range of interest as numerator
%
% wholefrange [low, high] frequnency range of interest as denominator
%             wholefrange must include frange.
%
% newRate     [Hz]. obj.Data is to be resamled to this sampling frequency
%
% window      The size of Hamming window. A positive even integer. The size of window in
%             points after resampling. The 50% of window will
%             be overlapped.
% nfft        the number of fft. A positive interger, must be
%             power of 2. Minimum is 256. Frequency resolution
%             of the resultant power spectrum is newRate/nfft.
%
%             examples:
%             newRate = 1024;
%             window = 512;
%             nfft = 512;
%
% OPTIONAL PARAMETER/VALUE PAIRS
% 'PlotType'  'line'   line drawing
%             'hist'   histogram
% 'XLim'      [left, right] in Hz. Default is [0 100]
%
% 'YLim'      [bottom, top]. Default is [0 1e-3]
%
% 'Preset'    'slow', 'spindles', 'beta', 'gamma'
%             override frange and wholefrange settings
%
% OUTPUT ARGUMENTS
% val         ratio of power (0 to 1)
%
% out         structure
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% kouichi.c.nakamura@gmail.com
% 03-Nov-2016 08:38:20




%% parse inputs

narginchk(3, 8); %TODO

[frange, wholefrange, newRate, window, nfft] = local_parse(obj, nargin, varargin{:});



%% job
x = resample(obj.Data, newRate, obj.SRate);
% target_title = obj.ChanTitle;

noverlap = window/2; %TODO

[Pxx,f] = pwelch(x, window, noverlap, nfft, newRate);

numerind = f >= frange(1) & f <= frange(2);

denomind = f >= wholefrange(1) & f <= wholefrange(2);

val = sum(Pxx(numerind)) / sum(Pxx(denomind));


%% output

out.Pxx = Pxx;
out.f = f;

out.window =window;
out.noverlap = noverlap;
out.nfft = nfft;
out.newRate = newRate;
out.frequencyRange = [0, newRate/2];
out.bins = []; %TODO
out.frequencyResolution = newRate/nfft; % [Hz]


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [frange, wholefrange, newRate, window, nfft] = local_parse(obj, nargin1, varargin)

if nargin1 > 3
        
    assert(nargin1 >= 6 && nargin1 <= 8, 'WaveformChan:powerRatio:nargin','nargin must be 6 to 8.');
    
    
    p = inputParser;
    
    vf1 = @(x) ~isempty(x) &&...
        isnumeric(x) &&...
        isrow(x) &&...
        numel(x) == 2 && ...
        all (x >= 0) &&...
        x(1) < x(2);
    
    addRequired(p, 'frange', vf1);
    
    vf2 = @(x) ~isempty(x) &&...
        isnumeric(x) &&...
        isrow(x) &&...
        numel(x) == 2 && ...
        all (x >= 0) &&...
        x(1) < x(2) &&...
        x(1) <= varargin{1}(1) && varargin{1}(2) <= x(2);
    
    addRequired(p, 'wholefrange', vf2);
    
    vf3 = @(x) ~isempty(x) &&...
        isnumeric(x) &&...
        isscalar(x) &&...
        (x > 0) && ...
        (x < obj.SRate);
    
    addRequired(p, 'newRate', vf3);
    
    vf4 = @(x) ~isempty(x) &&...
        isnumeric(x) &&...
        fix(x) == x &&...
        x > 0 &&...
        rem(x, 2) == 0; % positive even integer
    
    addRequired(p, 'window', vf4);
    
    vf5 = @(x) ~isempty(x) &&...
        isnumeric(x) &&...
        fix(x) == x &&...
        x >= 256 &&...
        rem(x, 2) == 0 && ...
        fix(log2(x)) == log2(x) ; % positive even integer, power of 2
    
    addRequired(p, 'nfft', vf5);
    
    parse(p, varargin{1:5});
    
    frange = varargin{1};
    wholefrange = varargin{2};
    newRate = varargin{3};
    window = varargin{4};
    nfft = varargin{5};
    
    
    if window * 10 > obj.Length
        warning('K:WaveformChan:powerRatio:window',...
            'Data length is shorter than 10*window. The FFT result may not be reliable.');
    end
    
end



%% Optional input arguments

ni = length(varargin); % ni >= 1
DataInputs = 0;
PNVStart = 0; % index of the first parameter in varargin
while DataInputs<ni && PNVStart==0
    nextarg = varargin{DataInputs+1};
    if ischar(nextarg) && isvector(nextarg)
        PNVStart = DataInputs+1;
    else
        DataInputs = DataInputs+1;
    end
end

preset = [];
if PNVStart > 0
    for i=PNVStart:2:ni
        % Set each Property Name/Value pair in turn.
        Property = varargin{i};
        if i+1>ni
            error('K:WaveformChan:powerRatio:options:pvsetNoValue', 'Value is missing')
        else
            Value = varargin{i+1};
        end
        
        % Perform assignment
        switch lower(Property)
            case 'preset'
                %% Assign the value
                if ~isempty(Value) && ischar(Value) && isrow(Value)
                    Value = validatestring(Value, {'slow', 'spindle', 'beta', 'gamma'});
                    preset = Value;
                else
                    error('K:WaveformChan:powerRatio:preset:pvsetInvalid', 'preset is invalid')
                end
                
            otherwise
                error('K:WaveformChan:powerRatio:option:pvsetInvalid', 'Parameter and values are invalid')
        end % switch
    end % for
    
    
    %% Preset parameters
    
    if ~isempty(preset)
        wholefrange = [0, 100];
        switch preset
            case 'slow'
                frange = [0.4, 1.6];
                newRate = 1024;
                window = 8192;
                nfft = 8192;
            case 'spindle'
                frange = [7, 12];
                newRate = 1024;
                window = 1024;% 512
                nfft = 1024;
                error('not implemented yet')
            case 'beta'
                frange = [15, 30];
                newRate = 1024;
                window = 1024;% 512
                nfft = 1024;
            case 'gamma'
                frange = [30, 60];
                newRate = 1024;
                window = 512;
                nfft = 512;
                error('not implemented yet')
        end
        
        if window * 10 > obj.Length
            warning('K:WaveformChan:powerRatio:window',...
                'Data length is shorter than 10*window. The FFT result may not be reliable.');
        end
    end
    
end

end


