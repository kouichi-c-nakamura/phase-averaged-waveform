function [f,t,cl,sc] = neurospec_WavWav_sp2a2_m1(obj, waveform, varargin)
%[f,t,cl,sc] = neurospec_WavWav_sp2a2_m1(obj, waveform, varargin)
%requires Neurospec 2.0

narginchk(2, inf);
    
p = inputParser;
vf = @(x) isa(x, 'WaveformChan');
addRequired(p, 'waveform', vf);
parse(p, waveform);


samp_rate = obj.SRate;
if obj.SRate ~= waveform.SRate
    error('SRate mismatch');
end   

seg_pwr = 10;
opt_str = 'h';

try
    [f,t,cl,sc] = sp2a2_m1(0, obj.Data, waveform.Data, samp_rate, seg_pwr, opt_str);
catch ME1
    disp('Have you installed Neurospec 2.0?');
    throw(ME1);
end

end