function tf = vf_structWaveform(x)
% tf = vf_structWaveform(x)
%
% tf     scalar logical
%
% validation function for field names of structure for event channel
% produced by Spike2

tf = false;

finamesW ={ 'title',...
    'comment',...
    'interval',...
    'scale',...
    'offset',...
    'units',...
    'start',...
    'length',...
    'values'};


if isstruct(x) && isscalar(x) && all(isfield(x, finamesW)) 
    
    if iscolumn(x.values)
        if isa(x.values, 'double') || isa(x.values, 'single') || isa(x.values, 'int16')
            tf = true;
        end
    end
end


