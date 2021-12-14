function h = plot(obj, varargin)
% h = plot(obj)
% h = plot(obj, 'Param', Value)
% h = plot(obj, ax, 'Param', Value, ...)
%      'Param', Value pairs for line function

narginchk(1, inf);

if nargin == 1
    fig = figure;
    ax = axes;
    
elseif nargin >=2
    if isgraphics(varargin{1},'axes')
        ax = varargin{1};
        if isscalar(ax) && ishandle(ax)
        
        if ~strcmp('axes', get(ax, 'Type'))
            error('K:EventChan:plot:ax:invalid',...
                'not valid axes handle.');
        end
        
        % axes(ax); %TODO SLOW
        fig = ancestor(ax,'figure');
        PVset = varargin(2:end);
        
        else
            error('K:EventChan:plot:ax:invalid',...
                'not valid axes handle.');
        end
    else
        fig = figure;
        ax = axes;
        PVset = varargin;
    end
end

xlim(ax, [obj.Start, obj.MaxTime]);
% custimization for EEG/LFPs
if any(strfind(lower(obj.ChanTitle), 'eeg')) ||...
        any(strfind(lower(obj.ChanTitle), 'ecog'))||...
        any(strfind(lower(obj.ChanTitle), 'lfp'))
    ylim(ax,[-0.5, 0.5]); 
end
set(ax, 'TickDir', 'out', 'Box', 'off');
xlabel(ax, 'Time (sec)');
ylabel(ax, sprintf('%s (%s)', obj.ChanTitle, obj.DataUnit),'Interpreter','none');

hold(ax,'on')
l1 = plot(ax, obj.time, obj.Data,'Tag','Waveform');
hold(ax,'off')

if exist('PVset', 'var') && ~isempty(PVset)
    set(l1, PVset{:});
end

h.fig = fig;
h.ax = ax;
h.l1 = l1;


end
