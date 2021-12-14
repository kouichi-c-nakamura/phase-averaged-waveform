function [data,S] = plotTriggered(obj,varargin)
%
%   [data,S]  = plotTriggered(obj,trigger,window_sec,offset_sec)
%   [data,S]  = plotTriggered(obj,ax,trigger,window_sec,offset_sec)
%   [data,S]  = plotTriggered(_____,ParamName,ParamValue)
%
% INPUT ARGUMETNS
% obj             A WaveformChan object
%
% trigger         A MetaEventChan object with the matching Start, Srate and
%                 Length with obj, or a column vector of 0 and 1 for
%                 trigger events, or a column vector of event time stamps
%
% window_sec      Window in seconds. Must be 0 or positive.
%
% offset_sec      Offset in seconds. Must be 0 or postive.
%
% ax              axes object
%
% OPTIONAL PARAMETER/VALUE PAIRS
% 'Average'       'on' (default) | 'off'
%                 for event channel this is equivalent to create PSTH, but
%                 you need to set histbin
%
% 'AverageColor'   Colorsepc for average waveform
%
% 'Overdraw'      'on' (default) | 'off'
%
% 'OverdrawColor'
%                 coloSpec for overdrawn waveforms
%
% 'OverdrawAlpha' Scalar double [0, 1]
%                 Default 0.2
%                 Alpha for overdrawn waveforms
%
% 'ErrorBar'      'off' (default) | 'std' | 'sem'
%
% 'ErrorBarColor' Colorspec
%
% 'ErrorBarAlpha' Scalar double [0, 1]
%                 Default 0.2
%
% 'Legend'       'on' (default) | 'off' %TODO
%
%
% OUTPUT ARGUMENTS
% axh        Handle of axes
% h_mean     Handle of line for mean
% h_err      Handle of patch for error
% h_overdr   Handle of lines for overdraw
% h_leg      Handle of legend
%
% data       data.t
%            data.mean
%            data.std
%            data.sem
%            data.segment   for overdraw
%
%            For the sake of memory and speed, if 'ErrorBar' is 'off',
%            data.std and data.sem will be empty.
%
% S          structure
%            With fields holding graphic handle objects:
%                 axh           axes
%                 line_mean     mean
%                 line_err      error
%                 line_overdr   overdraw
%                 leg      	    legend
%
% See also
% eventTriggeredAverage, linealpha, MarkerChan.plotTriggered

[ax,trigger,window_sec,offset_sec,average,averagecolor,errorbarmode,errobarcolor,...
    errorbaralpha,overdraw,overdrawcolor,overdrawalpha,showlegend] = parse(obj,varargin{:});

%% Job

Fs = obj.SRate;

if isempty(trigger)
    triggerevent = false(size(obj.Data)); % no trigger event
    
elseif isa(trigger,'MetaEventChan')
    triggerevent = trigger.Data;
elseif vf0or1(trigger)
    assert(length(trigger) == obj.Length);
    triggerevent = trigger;
    
elseif vftimestamps(trigger)
    assert(max(trigger) <= obj.MaxTime);
    assert(min(trigger) >= obj.Start);
    
    triggerevent = timestamps2binned(trigger,obj.Start,obj.MaxTime,obj.SRate,'ignore');

    %NOTE keep this comment
    %     ws = warning('error','K:timestamp2vec:yyNotBinary');
    %     try
    %         triggerevent = timestamps2binned(trigger,obj.Start,obj.MaxTime,obj.SRate);
    %     catch mexc
    %         if strcmp(mexc.identifier,'K:timestamp2vec:yyNotBinary')
    %             warning(mexc.message);
    %             disp('ignore option of timestamps2binned is used.')
    %             triggerevent = timestamps2binned(trigger,obj.Start,obj.MaxTime,obj.SRate,'ignore');
    %         else
    %             throw(mexc);
    %         end
    %     end
    %     warning(ws);
    
end

if strcmpi(errorbarmode,'off') && strcmpi(overdraw,'off') && strcmpi(average,'off')

    [t, outmean] = eventTriggeredAverage(obj.Data,...
        triggerevent,Fs,window_sec,offset_sec);
    %NOTE this often edns up in memory issue
    
    outstd = [];
    outsem = [];
    outseg = [];
    
else
    
    [t, outmean, outstd, outsem, outseg] = eventTriggeredAverage(obj.Data,...
        triggerevent,Fs,window_sec,offset_sec);
    %NOTE this often edns up in memory issue

end


if window_sec < 1
    T = t*1000;
else
    T = t;
end

%% Plot
data.t = T;
data.mean = outmean;
data.std = outstd;
data.sem = outsem;
data.segment = outseg;

if strcmpi(average,'off') && strcmpi(overdraw,'off') && strcmpi(errorbarmode,'off')
   % no plot to make
   ax = [];
   h_mean = [];
   h_err = [];
   h_overdr = [];
   h_leg = [];
else

    if ismember(lower(errorbarmode),{'std','sem'}) || strcmpi(overdraw,'on')...
            || strcmpi(average,'on')
        if ~isempty(ax)
            axes(ax);
        else
            figure;
        end
        
        if isa(trigger,'MetaEventChan')
            triggername = trigger.ChanTitle;
        else
            triggername = 'event';
        end
        
        if window_sec < 1
            xunit = 'ms';
        else
            xunit = 's';
        end
        
        xlabel(sprintf('Time relative to %s (%s)', triggername,xunit));
        ylabel(sprintf('Amplitude (%s)', obj.DataUnit));
        set(gca,'TickDir','out');
    end
    
    switch errorbarmode
        case 'std'
            hold on
            err = outstd;
            h_err = shadederrorbar(T,outmean,err,errobarcolor,...
                'FaceAlpha',errorbaralpha);
            set(h_err,'Tag','Shaded Error Bar STD');
        case 'sem'
            hold on
            err = outsem;
            h_err = shadederrorbar(T,outmean,err,errobarcolor,...
                'FaceAlpha',errorbaralpha);
            set(h_err,'Tag','Shaded Error Bar SEM');
        case 'off'
            h_err = [];
    end
    
    
    if strcmpi(overdraw,'on')
        
        segments = squeeze(outseg);
        if verLessThan('matlab','8.4.0')
            hold on
            h_overdr = zeros(size(segments,2),1);
        else
            hold on
            h_overdr = gobjects(size(segments,2),1);
        end
        
        for i = 1:size(segments,2)
            % transparent 'patchline' ... only works for hign magnification
            % h_overdr(i) = linealpha(T,segments(:,i),overdrawcolor,overdrawalpha);
            h_overdr(i) = line(gca,T,segments(:,i),'Color',[overdrawcolor,overdrawalpha]);
            
        end
        
        set(h_overdr,'Tag','Overdrawn Waveforms');
    else
        h_overdr = [];
    end
    
    if strcmpi(average,'on')
        hold on
        h_mean = plot(T,outmean,'Color',averagecolor,'Tag','Average Waveform');
    else
        h_mean = [];
    end
    
    hold off
    
    title(sprintf('%s',obj.ChanTitle),'Interpreter','none');
    
    %% legend
    
    if verLessThan('matlab','8.4.0')
        legtargets = [];
    else
        legtargets = gobjects(1,0);
    end
    legstr = {};
    
    if ~isempty(h_mean)
        legtargets = [legtargets,h_mean];
        legstr = [legstr,{'Average'}];
    end
    
    if ~isempty(h_err)
        legtargets = [legtargets,h_err(1)];
        legstr = [legstr,{sprintf('%s %s',char(177),upper(errorbarmode))}];
    end
    
    if ~isempty(h_overdr)
        dummy = line(0,0,'Color',overdrawcolor,'Visible','off','Tag','Dummy');
        legtargets = [legtargets,dummy];
        legstr = [legstr,{'Overdraw'}];
    end
    
    if verLessThan('matlab','9.2') %R2017a
        h_leg = legend(legtargets,legstr,'Location','NorthEast','Box','off');
    else
        h_leg = legend(legtargets,legstr,'Location','NorthEast','Box','off',...
            'AutoUpdate','off');
    end
    
    switch showlegend
        case 'off'
            h_leg.Visible = 'off';
    end
    
    ax = gca;

end

S.ax =          ax;
S.line_mean =   h_mean;
S.line_err =    h_err;
S.line_overdr = h_overdr;
S.leg =         h_leg;

end


%--------------------------------------------------------------------------

function [ax,trigger,window_sec,offset_sec,average,averagecolor,errorbarmode,errobarcolor,...
    errorbaralpha,overdraw,overdrawcolor,overdrawalpha,showlegend] = parse(obj,varargin)

narginchk(4,inf);
p = inputParser;
p.addRequired('obj',@(x) isa(x,'WaveformChan'));

if isscalar(varargin{1}) && ishandle(varargin{1}) && strcmp(get(varargin{1},'Type'),'axes')
    ax = varargin{1};
    args = varargin(2:end);
    
else
    ax = [];
    args = varargin(1:end);
end


p.addRequired('trigger',@(x) (isa(x,'MetaEventChan') && isscalar(x))...
    || vf0or1(x) || vftimestamps(x));

p.addRequired('window_sec',@(x) isscalar(x) && isreal(x) && x > 0);
p.addRequired('offset_sec',@(x) isscalar(x) && isreal(x) && x >= 0);

p.addParameter('Average','on',@(x) ismember(x,{'on','off'}));
p.addParameter('AverageColor','r',@iscolorspec);

p.addParameter('ErrorBar','off',@(x) ismember(x,{'std','sem','off'}));
p.addParameter('ErrorBarColor','k',@iscolorspec);
p.addParameter('ErrorBarAlpha',0.2,@(x) isreal(x) && isscalar(x) && x >=0 && x <= 1);

p.addParameter('Overdraw','on',@(x) ismember(x,{'on','off'}));
p.addParameter('OverdrawColor',[0.5 0.5 0.5],@iscolorspec);
p.addParameter('OverdrawAlpha',0.2,@(x) isreal(x) && isscalar(x) && x >=0 && x <= 1);

p.addParameter('Legend','on',@(x) ismember(x,{'on','off'}));

p.parse(obj,args{:});

trigger = p.Results.trigger;
window_sec = p.Results.window_sec;
offset_sec = p.Results.offset_sec;

average = p.Results.Average;
averagecolor = p.Results.AverageColor;
if ischar(averagecolor)
    averagecolor = lower(averagecolor);
end

errorbarmode = p.Results.ErrorBar;
errobarcolor = p.Results.ErrorBarColor;
if ischar(errobarcolor)
    errobarcolor = lower(errobarcolor);
end
errorbaralpha = p.Results.ErrorBarAlpha;


overdraw = p.Results.Overdraw;
overdrawcolor = p.Results.OverdrawColor;
if ischar(overdrawcolor)
    overdrawcolor = lower(overdrawcolor);
end
overdrawalpha = p.Results.OverdrawAlpha;

showlegend = p.Results.Legend;

end

%--------------------------------------------------------------------------

function out = vf0or1(x)
% vf0or1 = @(x) iscolumn(x) && all(x == 0 | x == 1);

out = iscolumn(x) && all(x == 0 | x == 1);

end

%--------------------------------------------------------------------------

function out = vftimestamps(x)
% vftimestamps = @(x) iscolumn(x) && isreal(x) && all(diff(x) > 0);

out = iscolumn(x) && isreal(x) && all(diff(x) > 0);

end