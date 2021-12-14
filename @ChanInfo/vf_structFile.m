function tf = vf_structFile(x)
% tf = vf_structFile(x)
%
% tf     scalar logical
%
% validation function for field names of structure for event channel
% produced by Spike2
%
% Spike2 version 8 MATLAB export include a structure named file
tf = false;

finamesM ={'name',...
            'start'};
        
if isstruct(x) && isscalar(x) &&...
            all(isfield(x, finamesM))
    tf = true;
end