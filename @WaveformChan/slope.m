function slp = slope(obj,psec)
% WaveformChan.slope returns a WaveformChan object that holds a slope of
% the obj
%
% SYNTAX
% slp = slope(obj,psec)
%
% INPUT ARGUMENTS
% obj         WaveformChan object
%
% psec        a positive scalar
%             Time constant in seconds
%
%
% OUTPUT ARGUMENTS
% slp         WaveformChan object
%             Holding a slope of obj. slp.DataUnit is [obj.DataUnit,'/s']
%
% from Spike 2 version 8 documentation
%
% This process has one argument, a time period in seconds, p. The slope at
% time t is calculated using an equal weighting of the points from time t-p
% to t+p. This is done by calculating the mean of the points ahead and the
% mean of the points behind each data point, and the slope is taken from
% the line through the centre of the points behind to the centre of the
% points ahead. This calculation is equivalent to applying an FIR filter,
% it is not a least squares fit through the points.
%
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% kouichi.c.nakamura@gmail.com
% 24-Nov-2018 07:43:47
%
% See also
% WaveformChan


w = obj.Data;
Fs = obj.SRate;
n = obj.Length;

ppt = round(psec*Fs);

SLP = NaN(size(w));

for t = 1:n
    
    if all((t - ppt):t > 0) && all((t - ppt):t <= n)
        befmean = mean(w((t - ppt):t));
    else
        befmean = NaN;
    end
    
    if all(t:(t + ppt) > 0) && all(t:(t + ppt) <= n)
        aftmean = mean(w(t:(t + ppt)));
    else
        aftmean = NaN;
    end
    
    SLP(t) = (aftmean - befmean)/ (ppt*1/Fs);
    
end

slp = WaveformChan(SLP,obj.Start,obj.SRate,[obj.ChanTitle,' slope']);
slp.DataUnit = [obj.DataUnit,'/s'];


end