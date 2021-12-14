function obj = struct2ChanInfo(obj, S, varname)
% obj = Struct2ChanInfo(obj, S, varname)
%
% protected function
% used for the constructor ChanInfo()
% construct obj from struct input S
% varname is the variable name in workspace
% EventChan and WaveformChan

%% parse

narginchk(2, 3);

if nargin == 2
    varname = '';
end

p = inputParser;
vf_S = @(x) isstruct(x) && ...
    ~isempty(fieldnames(x));
vf_varname = @(x) ischar(x) &&...
    isrow(x);
addRequired(p, 'S', vf_S);
addRequired(p, 'varname', vf_varname);
parse(p, S, varname);

%% job

% obj = ChanInfo(); %TODO this prevent subclass from calling superclass constructor

if ChanInfo.vf_structWaveform(S)
    
    obj.DataUnit = S.units;
    
elseif ChanInfo.vf_structEvent(S)
    
    obj.DataUnit = '';
    
elseif ChanInfo.vf_structMarker(S)
    error('K:ChanInfo:MarkerChan:struct:structBin',...
        ['To get some information on a marker channel, you need ',...
        'structMk (not binned) for the marker channel and structBin as a reference ',...
        'channel in binned format.']);
else
    error('K:ChanInfo:struct:invalid',...
        'strucut doesn''t seem to be in Spike2 format.');
end

obj.Start = S.start;
obj.SRate = 1/S.interval;
obj.ChanTitle = S.title;
obj.ChanStructVarName = varname;
obj.Header = rmfield(S, 'values');
obj.Length = S.length;

if isfield(S, 'channumber')
    obj.ChanNumber = S.channumber;
else
    obj.ChanNumber = 0;
end

end