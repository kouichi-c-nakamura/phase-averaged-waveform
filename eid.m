function eid = eid(varargin)
% eid = eid()
% eid = eid(string)
% Local function that generates error id that begins with K:
% NOTE: onlyt to be used as a local function.
%
% input argument
% str (Optional) string in char type (row vector)
%
% output argument
% eid an error id composed of 'K:(functionname):str'
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% kouichi.c.nakamura@gmail.com
% 12-Dec-2016 16:17:33



narginchk(0, 1);
p = inputParser;
p.addOptional('str', '', @(x) isempty(x) || ischar(x) && isrow(x));
p.parse(varargin{:});
str = p.Results.str;

if isempty(str)
    str = '';
else
    str = [':', str];
end

[~,m,~] = fileparts(mfilename('fullpath'));

eid = ['K:', m, str];


end