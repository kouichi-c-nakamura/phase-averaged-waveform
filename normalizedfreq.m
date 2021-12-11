function nF = normalizedfreq(freq, Fs)
% nF = normalizedfreq(freq, Fs)
%
% Note:   This toolbox (Signal Processing Toolbox) uses the convention that unit frequency 
% is the Nyquist frequency, defined as half the sampling frequency. The cutoff frequency
% parameter for all basic filter design functions is normalized by the Nyquist frequency.
%
% For a system with a 1000 Hz sampling frequency, for example, 300 Hz is 300/500 = 0.6. 
%
% To convert normalized frequency to angular frequency around the unit circle, multiply by
% pi. To convert normalized frequency back to hertz, multiply by half the sample frequency.
% http://mathworks.com/help/signal/ug/frequency-response.html
%
% INPUT ARGUMENTS
% freq      cycles/sec or Hz; scaler or 2 element vector
% Fs        sampling frequency [Hz]
%
% OUTPUT ARGUMENT
% nF        normalized frequency [0 1] in "half-cycle per sample"
%           where Nyquist frequency [pi radian/sample] corresponds to 1
%
% ANONYMUS FUNCTION FORM
%    getNF = @(freq, Fs) freq/(Fs/2);
%
%
%% cf. OTHER DEFINITIONS
% nF = freq*2*pi/Fs; % [radians/sample] expresion
%     (cycles/second) / (samples/second) = cycles/sample. 
%
% nF = freq/Fs % [cycle per sample] expression
%
%
% NOTE: There are other defintions of normalized frequency!!! DO NOT CONFUSE
%
% "Some programs (such as MATLAB) that design filters with real-valued 
% coefficients use the Nyquist frequency (\scriptstyle f_s/2) as the normalization
% constant. The resultant normalized frequency has units of [half-cycles/sample] or
% equivalently [cycles per 2 samples]."
% https://en.wikipedia.org/wiki/Normalized_frequency_(unit)
%
%    omega = 2 * pi * f/Fs
%
% where omega is normalized frequency, f is physical frequency [Hz], and Fs is sampling frequency [Hz}.
% http://www.mathworks.co.uk/help/signal/ug/spectral-analysis.html
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% kouichi.c.nakamura@gmail.com
% 09-Dec-2016 20:27:47




narginchk(2,2)

p = inputParser;
p.addRequired('freq', @(x) all(isreal(x)) && all(x >= 0) && ...
 isscalar(x) || numel(x) == 2 && isrow(x) ) % frequency of interest
p.addRequired('Fs',  @(x) isscalar(x) && isreal(x) && x >= 0) % sampgling frequency
p.parse(freq, Fs);

%% Job

nF = freq/(Fs/2);    % == 2*freq/Fs  [half-cycle per sample]
% == freq/NiquistFreq
