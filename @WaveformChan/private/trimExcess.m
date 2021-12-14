function X = trimExcess(x,n)
% X = trimExcess(x,n)
%
% x(x < 1) = [];
% x(x > n) = [];
% 
% X = x;
%
% See also
% WaveformChan.getBUA, startsends2ind

x(x < 1) = [];
x(x > n) = [];

X = x;

end