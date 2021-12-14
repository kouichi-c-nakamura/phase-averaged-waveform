function [ind,time, varargout] = crossThreshold(obj,direction,threshold,varargin)
% WaveformChan.crossThreshold allows you to find data points where the
% curve crosses the threshold by rising or falling or both. Thsi method is
% a wrapper of the crossThreshold function. crossThreshold mimics the the
% behaviours of Spike2 active cursor's rising/falling threshold search.
%
% SYNTAX
% [ind,time] = crossThreshold(obj,direction,threshold)
% [ind,time] = crossThreshold(____,'Param',value)
% [ind,time,evt] = crossThreshold(____)
%
%
% INPUT ARGUMENTS
% obj         WaveformChan object
%
% direction   'rising' | 'falling' | 'any'
%
%
% threshold   scalar real value
%             Threshold value
%
% OPTIONAL PARAMETER/VALUE PAIRS
%
% 'Fs'        (Default) 1
%             Sampling frequency. Ingored if wave is a WaveformChan object.
%
% 'minDelay'  Minimum duration of corrsing in seconds when Fs is provided,
%             or in data points when Fs is omitted.
%
% 'minCross'  Minimum amount of crossing in terms of amplitude
%
%
% OUTPUT ARGUMENTS
% ind         Indices of crossing points. Crossing points are the points
%             immediately after the ideal crossing poits. 
%
% time        Time stamps for crossing points. Only meaningful when you set
%             Fs properly (without Fs, time is identical to ind).
%
% evt         EventChan object
%             An EventChan object representing the rising or falling data
%             points with ChanTitle 'rising' or 'falling' and
%             Header.threshold for threshold
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% kouichi.c.nakamura@gmail.com
% 24-Nov-2018 18:11:47
%
% See also
% WaveformChan, crossthreshold

[ind,time] = crossthreshold(direction,obj,threshold,varargin{:});

if nargout == 3 
    
    data = false(obj.Length,1);
    data(ind) = true;
        
    evt = EventChan(data,obj.Start, obj.SRate,direction);
    
    evt.Header.threshold = threshold;
    
    varargout{1} = evt;
    
elseif nargout > 3
    
    error('Too many output arguments')

end
   

end

