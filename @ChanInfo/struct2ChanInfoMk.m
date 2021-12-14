function obj = struct2ChanInfoMk(obj, S, Sbin, varname1, varname2)
% obj = struct2ChanInfoMk(obj, S, Sbin, varname1, varname2)
%
% protected function
% used for the constructor ChanInfo()
% construct obj from struct input S
% Sbin is struct from a binned channel for reference

%% parse

narginchk(3,5);

if nargin == 3
    varname1 = '';
    varname2 = '';
end

if nargin == 4
    error('K:ChanInfo:MarkerChan:struct:nargin',...
        'nargin must be 3 or 5');
end

p = inputParser;
vf_S = @(x) isstruct(x) && ...
    ~isempty(fieldnames(x));
vf_varname = @(x) ischar(x) &&...
    isrow(x);
addRequired(p, 'S', vf_S);
addRequired(p, 'Sbin', vf_S);

addRequired(p, 'varname1', vf_varname);
addRequired(p, 'varname2', vf_varname);

parse(p, S, Sbin, varname1, varname2);


%% job

% obj = ChanInfo(); %TODO this doesn't allow subclass

assert(ChanInfo.vf_structMarker(S),...
    'K:ChanInfo:MarkerChan:struct:structMk',...
    ['To get some information on a marker channel, you need ',...
    'structMk (not binned) for the marker channel and structBin as a reference ',...
    'channel in binned format.']);

assert(ChanInfo.vf_structWaveform(Sbin) || ChanInfo.vf_structEvent(Sbin), ...
    'K:ChanInfo:MarkerChan:struct:structBin',...
    ['To get some information on a marker channel, you need ',...
    'structMk (not binned) for the marker channel and structBin as a reference ',...
    'channel in binned format.']);

obj.ChanTitle = S.title;
obj.ChanStructVarName = varname1;

obj.ChanTitleRef = Sbin.title;
obj.ChanStructVarNameRef = varname2;

obj.DataUnit = '';
obj.Start = Sbin.start;
obj.SRate = 1/Sbin.interval;
obj.Length = Sbin.length; %TODO need check


obj.Header = S;
if isfield(S, 'channumber')
    obj.ChanNumber = S.channumber;
else
    obj.ChanNumber = [];
end
end