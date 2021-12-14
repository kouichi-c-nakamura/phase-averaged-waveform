function ind = startsends2ind(startInd,endInd,n)
% startsends2ind returns a column vector of numeric indices for eposchs whose
% onset and offset specified by startInd and endInd, respectively, against
% a vector of n length. In case, those epochs overlap with each other, the
% epochs will be joined up.
%
% ind = startsends2ind(startInd,endInd,n)
%
% startInd     column vector of index numbers of an array defining the
%              onset of specific periods
%
% endInd       Similar to the above but for offset
%
% OUTPUT ARGUMENTS
% ind          Column vector of indices of epoches defined by startInd 
%              and endInd. Excess indices (ind < 1 or ind >N) will be removed.
%
% See also
% WaveformChan.getBUA, trimExcess


vfscposint = @(x) isscalar(x) && x > 0 && all(fix(x) == x);
vfint = @(x) all(fix(x) == x);

p = inputParser;
p.addRequired('startInd',vfint);
p.addRequired('BendInd',vfint);
p.addRequired('n',vfscposint);
p.parse(startInd,endInd,n);

%% Job

C1 = arrayfun(@(x,y) x:y, startInd, endInd, 'UniformOutput',false)';
C2 = cellfun(@(x) trimExcess(x,n), C1,'UniformOutput',false);

ind = unique([C2{:}]'); % avoid repeative referencing by numeric indices

end
