function status = exportCSV(obj, varargin)
%
% status = exportCSV(obj)
% status = exportCSV(obj, start, end)
% status = exportCSV(_______, 'Param', value )

if nargin == 2
    error('K:WaveformChan:exportCSV', 'start needs to be provided with end')
end

p = inputParser;
p.addRequired('obj');
p.addOptional('start', 0, @(x) isscalar(x) && isreal(x));
p.addOptional('end', 0, @(x) isscalar(x) && isreal(x));
p.parse(obj, varargin{:});

startT = p.Results.start;
endT = p.Results.end;

%TODO startT and endT as second, or should it be index of Data?



fid = openfile();

% when invoked by Record, you need to skip this
if skipfileinfo
    fprinf(fid, ['###\n', ...
        'filename, "%s"\n',...
        'SampligRate, %f',...
        'start, %f\n',...
        'maxtime, %f\n'],...
        ...
        'New',...
        obj.SRate,...
        obj.Start,...
        obj.MaxTime);
end
     
fprintf(fid, ['##',...
    'chantitle, "%s"' ,...
    'chankind, "%s"',...
    'scale, %f',...
    'offset, %f'],...
    ...
    obj.ChanTitle,...
    'event',...
    obj.Scale,...
    obj.Offset);


% Int16toW = @(x) double(x).*obj.Scale + obj.Offset;

WtoInt16 = @(x) int16((x - obj.Offset)/obj.Scale);



data = obj.Data(startT:endT);

fprintf(fid, ['#',...
    '%f\n',...
    ...
    WtoInt16(obj.Data)]);
    % cf. ChanWriteWave() in Spike2


end

