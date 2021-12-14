function [ds, h] = detectUPandDOWNwavelet(obj, varargin)
% ds = detectUPandDOWNwavelet(obj)
% ds = detectUPandDOWNwavelet(obj, lowpass, phase_threshold, time_threshold)
% ds = detectUPandDOWNwavelet(obj, ______, doplot)
%
% INPUTS
% lowpass               [Hz] Cutoff frequency for low pass (default 2)
% phase_threshold       between -1 and 1 (dafault is 0);
% time_threshold        [msec] for states (dafault is 250 msec);
%
% Param/Value pairs
% 'doplot'              true of false (default)
%
% OUTPUTs
% ds                    dataset containing thre retuls in second
% h                     handle for graphics


narginchk(1,6)


%% PNVStart
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


%% Parse inputs

lowpass =2;
phase_threshold =0;
time_threshold = 250;
doplot = false;

if nargin == 1 || PNVStart == 1
elseif PNVStart >= 2
    lowpass = varargin{1};
    if PNVStart >= 3
        phase_threshold =varargin{2}; 
        if PNVStart == 4
            time_threshold = varargin{3};
        end
    end
end

p = inputParser;

vf1 = @(x) ~isempty(x) &&...
    isnumeric(x) &&...
    isscalar(x) &&...
    x > 0;

addRequired(p, 'lowpass', vf1);

vf2 = @(x) ~isempty(x) &&...
    isnumeric(x) &&...
    isscalar(x) &&...
    -1 <= x && x <= 1;

addRequired(p, 'phase_threshold', vf2);

vf3 = @(x) ~isempty(x) &&...
    isnumeric(x) &&...
    isscalar(x) &&...
    x > 0;

addRequired(p, 'time_threshold', vf3);

parse(p, lowpass, phase_threshold, time_threshold);


%% Param/Value pairs

if PNVStart > 0
    for i=PNVStart:2:ni
        % Set each Property Name/Value pair in turn.
        Property = varargin{i};
        if i+1>ni
            error('id:pvsetInvalid', 'message')
        else
            Value = varargin{i+1};
        end
        % Perform assignment
        switch lower(Property)
            case 'doplot'
                %% Assign the value
                if ~isempty(Value) && isscalar(Value) && ...
                        islogical(Value) || ...
                        isnumeric(Value) && Value == 1 || Value ==0;
                    % Name has been specified
                    doplot = Value;
                else
                    error('K:WaveformChan:detectUPandDOWNwavelet',...
                        'Value for doplot parameter must be logical or numeric 0 or 1');
                end
                
            otherwise
                error('id:pvsetInvalid', 'message')
        end % switch
    end % for
end

%% job

W = resample(obj.Data, 1000, round(obj.SRate));

states= K_UPDOWNwavelet(W, lowpass, phase_threshold, time_threshold);

ds = mat2dataset(states/1000, 'VarNames', {'UPstart','UPend','UPdur','DOWNstart','DOWNend','DOWNdur', 'skip'});


h.fig1 =[];h.l1=[];h.ax1=[];
if doplot
    t = (1:length(W))/1000;
    h.fig1 = figure; 
    hold on;
    h.l1 = plot(t, W, 'g');
    h.ax1 = gca;
    
    box off;set(gca, 'TickDir', 'out');
    xlabel('Time [sec]');
    pan xon; zoom xon;
    
    ylim([-1, 1]);
    ylabel(obj.ChanTitle);
    
    for i = 1:size(states, 1)
        plot([ds.DOWNstart(i), ds.DOWNend(i)], [-0.5, -0.5], 'b');
        plot([ds.DOWNstart(i), ds.DOWNstart(i)], [-0.5, 0.5], 'Color', [0.5, 0.5, 0.5]);
        plot([ds.DOWNend(i), ds.DOWNend(i)], [-0.5, 0.5], 'Color', [0.5, 0.5, 0.5]);
        
        plot([ds.UPstart(i), ds.UPend(i)], [0.5, 0.5], 'r');
        plot([ds.UPstart(i), ds.UPstart(i)], [-0.5, 0.5], 'Color', [0.5, 0.5, 0.5]);
        plot([ds.UPend(i), ds.UPend(i)], [-0.5, 0.5], 'Color', [0.5, 0.5, 0.5]);
    end
    
end

end