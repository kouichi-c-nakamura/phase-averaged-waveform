function obj2 = extractTime(obj, StartTime, EndTime, mode, time)
% Chan.extractTime() is a method to extract fragment of data between the
% specified time window.
%
% SYNTAX
% obj2 = extracttime(obj, StartTime, EndTime)
% obj2 = extracttime(obj, StartTime, EndTime, mode)
% obj2 = extracttime(obj, StartTime, EndTime, time)
%
% INPUT ARGUMENTS
% StartTime, EndTime      
%             in second
%
% mode        'normal' (default) | 'extend' 
%             (OPTIONAL) 'extend' mode accepts StartTime and EndTime
%             outside of the range of Time vector and fill the gap with NaN
%             or zeros.
%
% time        column vector of time stamps
%             (OPTIONAL) time stamps for all the data points. To be used to
%             avoid repeated call of Chan.time method, which is slow.
%
%
% OUTPUT ARGUMENTS
% obj2        an object of the same class as obj 
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% kouichi.c.nakamura@gmail.com
% 21-Dec-2018 16:35:04
%
% See also
% Chan, WaveformChan, EventChan

%% parse

arguments
    
    obj
    StartTime
    EndTime
    mode (1,:) char {mustBeMember(mode,{'normal','extend'})} = 'normal';
    time (:,1) = zeros(1, 0)
end

if strcmpi(mode, 'normal')
    if obj.Start > StartTime || obj.MaxTime < EndTime
        error('K:Chan:extractTime:needToExtend',...
            'You need ''extend'' option to accept StartTime or EndTime outside of the Time.');
    end
end

p = inputParser;

switch mode
    case 'normal'
        
        vf_StartTime = @(x) ~isempty(x) &&...
            isscalar(x) && ...
            isreal(x) && ...
            isfinite(x) && ...
            ~isnan(x) && ...
            x <= EndTime && ...
            x >= obj.Start &&...
            x <= obj.MaxTime;
        addRequired(p, 'StartTime', vf_StartTime);
        
        vf_EndTime = @(x) ~isempty(x) &&...
            isscalar(x) && ...
            isreal(x) && ...
            isfinite(x) && ...
            ~isnan(x) && ...           
            x >= obj.Start &&...
            x <= obj.MaxTime;
        addRequired(p, 'EndTime', vf_EndTime);
        
    case 'extend'
        vf_StartTime = @(x) ~isempty(x) &&...
            isscalar(x) && ...
            isreal(x) && ...
            isfinite(x) && ...
            ~isnan(x) && ...
            x <= EndTime;
        addRequired(p, 'StartTime', vf_StartTime);
        
        vf_EndTime = @(x) ~isempty(x) &&...
            isscalar(x) && ...
            isreal(x) && ...
            isfinite(x) && ...
            ~isnan(x);
        addRequired(p, 'EndTime', vf_EndTime);
        
end

parse(p, StartTime, EndTime);

L = round((EndTime - StartTime)/obj.SInterval) + 1;

if isempty(time)
    time = obj.time; %SLOW
else
    assert(obj.Length == length(time))
    assert(obj.Start == time(1))
    assert(obj.MaxTime == time(end))    
end

%% job
switch lower(mode)
    case 'normal'
        % find the Time indices
        
        
        startind = find(time <= StartTime, 1, 'last');
        
        endind = startind + L - 1;% find(time >= EndTime, 1, 'first');
        
        obj2 = obj;
        if obj2.Length < endind
            endind = obj2.Length; 
        end
        
        obj2.Data = obj2.Data(startind:endind, 1);
        obj2.Start = time(startind);
        obj2.Header = obj.Header;
        
    case 'extend'
        
        if StartTime < obj.Start
            startind = 1;
            if obj.Start < EndTime
                startpad = round((obj.Start - StartTime)/obj.SInterval) + 1;
            else
                startpad = round((EndTime - StartTime)/obj.SInterval) + 1;
            end
            
            obj_bef = obj;
            if isa(obj_bef, 'WaveformChan')
                obj_bef.Start = 0;
                obj_bef.Data = NaN(startpad, 1);                
            elseif isa(obj_bef, 'EventChan')
                 obj_bef.Start = 0;
                 obj_bef.Data = zeros(startpad, 1);
            elseif isa(obj_bef, 'MarkerChan')
                 obj_bef.Start = 0;
                 obj_bef.Data = zeros(startpad, 1); 
                 obj_bef.TextMark = cell(0,1);
                 obj_bef.MarkerCodes = zeros(0,4);
            else
                error('K:Chan:extractTime:extend:isa',...
                    'unexpected class for obj %s', inputname(1));
            end
            
        else
            startind = find(time <= StartTime, 1, 'last');
            startpad = 0;
        end
        
        if EndTime > obj.MaxTime
            endind = obj.Length;
            endpad = round((EndTime - obj.MaxTime)/obj.SInterval);
            
            obj_aft = obj;
            if isa(obj_aft, 'WaveformChan')
                obj_aft.Start = 0;
                obj_aft.Data = NaN(endpad, 1);                
            elseif isa(obj_aft, 'MetaEventChan')
                 obj_aft.Start = 0;
                 obj_aft.Data = zeros(endpad, 1);
            else
                error('K:Chan:extractTime:extend:isa',...
                    'unexpected class for obj %s', inputname(1));
            end
            
        else
            if obj.Start < EndTime 
                if obj.Start <= StartTime
                    endind = startind + L -1;
                else
                    endind = startind + L -1 - startpad;                    
                end
                
            else
                endind = [];                
            end
            endpad = 0;
        end

        obj2 = obj; % copy parameters
        obj2.Data = obj2.Data(startind:endind, 1);
        obj2.Start = time(startind);
        obj2.Header = obj.Header;
        
        if startpad > 0
            obj2 = [obj_bef; obj2];
            obj2.Start = StartTime;
        end
        
        if endpad > 0
            obj2 = [obj2; obj_aft];
        end
                
end


end
