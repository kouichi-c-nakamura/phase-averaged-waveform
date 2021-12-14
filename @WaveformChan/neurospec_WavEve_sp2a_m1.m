function [f,t,cl,sc] = neurospec_WavEve_sp2a_m1(obj, eventchan, varargin)
%[f,t,cl,sc] = neurospec_WavEve_sp2a_m1(obj, eventchan, varargin)
%requires Neurospec 2.0


narginchk(2, inf);


p = inputParser;
vf = @(x) isa(x, 'EventChan');
addRequired(p, 'eventchan', vf);
parse(p, eventchan);

try
    [f,t,cl,sc] = sp2a_m1(0,eventchan.TimeStamps,obj.Data,varargin);
catch ME1
    disp('Have you installed Neurospec 2.0?');
    throw(ME1);
end

end