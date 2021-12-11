function h = shadederrorbar(varargin)
%
% shadederrorbar draws a shaded error range using a patch object.
%
% h = shadederrorbar(x,y,err)
% h = shadederrorbar(x,y,err,facecolor)
% h = shadederrorbar(x,y,err,facecolor,edgecolor)
% h = shadederrorbar(x,y,{neg,pos},_____)
% h = shadederrorbar(ax, _____)
% h = shadederrorbar(______, 'Name',Value)
%
% INPUT ARGUMENTS
% x              x axis values in a vector
%
% y              Typically the mean values for y axis in a vector format.
%                As long as NaNs appear at the same data points in y
%                and err, they can be accepted.
%
% err            a vecotr of non-negative numbers | {neg, pos}
%                The error range for y axis in a vector, typically standard
%                deviation or standard error of means.
%
%                Alternatively, you can use a cell vector {neg, pos} as
%                err, where neg and pos are the length of y error range (in
%                the same format as err as vector), in order to specify
%                different error ranges for postive and negative values.
%                For example, you can use medians for y, and the distance
%                from medians to the 25 and 75 percentiles for neg and pos
%                to show interquarter range (IQR).
%
% The size of x, y, and err (if err is cell array, the size of neg and pos) 
% must be identical.
%
%
% facecolor      (Optional) FaceColor property of the patch object.
%                By default it uses color that is specified by ColorOrder and
%                ColorOrderIndex propertis of axes (R2014b or later), or 'k'
%                (earlier releases).
%
% edgecolor      (Optional) EdgeColor property of the patch object. By 
%                default it is set to 'none'
%
%
% OPTINAL PARAMETER/VALUE PAIRS
%
% 'FaceAlpha'    FaceColor property of the patch object.
%                Default is 0.5.
%
% 'Params'       A row vector of cell array containing property names and values 
%                for patch object.
%
%                e.g.   'params', {'FaceAlpha',0.3,'EdgeAlpha',0.3}
%
%
% OUTPUT ARGUEMENTS
% h          Graphic handle of the patch object
%
% Writteny by
% Kouichi C. Nakamura, Ph.D
% Dept Morhological Brain Science
% Kyoto University, Kyoto, Japan
% kouichi.c.nakamura@gmail.com
%
% 9 Nov 2015
%
% See also
% errorbar, patch, shadederrorbar_demo
% https://uk.mathworks.com/help/matlab/creating_plots/_mw_322212e1-ae64-481a-bce6-b9fd5120b5fb.html
% (Line Plot with Confidence Bounds)

%% Parse

[ax,x,y,neg,pos,facecolor,edgecolor,facealpha,params] = parse(varargin{:});

%% Job

if isempty(facecolor)
    if ~verLessThan('matlab','8.4.0')
        colors = get(gca,'ColorOrder');
        indexBefore = get(gca,'ColorOrderIndex');

        ind = mod(indexBefore,size(colors,1));

        facecolor = colors(ind,:);
    else
        facecolor = 'k';
    end
end


neg(neg == 0) = eps; %NOTE Workround for weird shape related to err == 0
pos(pos == 0) = eps; %NOTE Workround for weird shape related to err == 0


if ~isempty(ax)
    axes(ax)
else
    ax = gca;
end

[~ ,negC,~ ] = splitNumbers(y,neg,x);
[yC,posC,xC] = splitNumbers(y,pos,x);

hold(ax,'on')

h = gobjects(length(yC),1);

for i = 1:length(yC)
    h(i) = fill(ax,[xC{i}, fliplr(xC{i})], ...
        [yC{i} + posC{i}, fliplr(yC{i} - negC{i})], ...
        facecolor, 'EdgeColor', edgecolor,'FaceAlpha', facealpha, 'linestyle', 'none',...
        'Tag','Shaded Error Bar'); %TODO is patch better than fill?
end


if ~isempty(params)

    set(h, params{:});

end

if isempty(facecolor)
    if ~verLessThan('matlab','8.4.0')
        set(ax,'ColorOrderIndex',indexBefore + 1);
    end
end


end

%--------------------------------------------------------------------------

function [ax,x,y,neg,pos,facecolor,edgecolor,facealpha,params] = parse(varargin)

p = inputParser;

if isscalar(varargin{1}) && ishghandle(varargin{1}) && ...
        isgraphics(varargin{1},'Axes')
    ax      = varargin{1};
    x        = varargin{2};
    y        = varargin{3};
    err      = varargin{4};
    paramval = varargin(5:end);

else
    ax      = [];
    x        = varargin{1};
    y        = varargin{2};
    err      = varargin{3};
    paramval = varargin(4:end);

end


vf1 = @(x) isvector(x) && isreal(x);
p.addRequired('x', vf1);
p.addRequired('y', vf1);
p.addRequired('err', @(x) isvector(x) && isreal(x) || ...
    ( iscell(x) && isrow(x) && numel(x) && all(cellfun(vf1,x)))) ;

vfcolor = @(x) iscolorspec(x) || ischar(x) && isrow(x) && strcmpi(x,'none');
p.addOptional('facecolor','',vfcolor);
p.addOptional('edgecolor','none',vfcolor);
p.addParameter('Params',{},@(x) iscell(x) && isrow(x) && rem(numel(x),2) == 0);
p.addParameter('FaceAlpha',0.5,@(x) isscalar(x) && x >= 0 && x <= 1);
p.parse(x,y,err,paramval{:});

facecolor = p.Results.facecolor;
edgecolor = p.Results.edgecolor;
facealpha = p.Results.FaceAlpha;
params    = p.Results.Params;


assert(isequal(size(x),size(y)) || isequal(size(x),size(y')));

if isequal(size(x),size(y'))
    y = y';
else
    assert(isequal(size(x),size(y)));
end

if iscell(err)
    neg = err{1};
    if isequal(size(x),size(neg'))
        neg = neg';
    else
        assert(isequal(size(x),size(neg)),'Size of x and err{1} does not match');
    end

    pos = err{2};
    if isequal(size(x),size(pos'))
        pos = pos';
    else
        assert(isequal(size(x),size(pos)),'Size of x and err{2} does not match');
    end

else
    if isequal(size(x),size(err'))
        neg = err';
        pos = neg;
    else

        assert(isequal(size(x),size(err)),'Size of x and err does not match');
        neg = err;
        pos = err;
        
    end
end


if iscolumn(x)
    x = x';
    y = y';
    pos = pos';
    neg = neg';
end
end

%--------------------------------------------------------------------------

function [yC,errC,xC] = splitNumbers(y,err,x)
%
% [yC,errC,xC] = splitNumbers(y,err,x)
%
% splitNumbers takes out NaNs and enables splitted shaded error bars for data
% including NaN values.
%
% splitNumbers checks if y or err contains NaN, and if one of them
% does, take out fragments of non-NaN numbers and put them in cell vector.
%
% See also
% isnan

tfnan_m = isnan(y);
tfnan_e = isnan(err);

assert(isequal(tfnan_m,tfnan_e),...
    'shadederrorbar:splitNumbers:isnan',...
    'NaN appears in different parts in y and err')

if any(tfnan_m)
    % split

    diffvec = diff(tfnan_m);

    endofNum = find(diffvec == 1);
    startofNum = find(diffvec == -1) +1;

    if isempty(endofNum) && isempty(startofNum)
        % all NaN
        assert(all(tfnan_m))
        yC = {NaN};
        errC = {NaN};
        xC = {NaN};
    elseif isempty(endofNum)
        % starts with NaNs and ends with numbers
        assert(length(startofNum) == 1)

        yC = {y(startofNum:end)};
        errC = {err(startofNum:end)};
        xC = {x(startofNum:end)};

    elseif isempty(startofNum)
         % starts with numbers and ends with NaNs
        assert(length(endofNum) == 1)

        yC = {y(1:endofNum)};
        errC = {err(1:endofNum)};
        xC = {x(1:endofNum)};

    else

        if endofNum(1) < startofNum(1)
            startofNum = [1,startofNum];
        end

        if endofNum(end) < startofNum(end)
            endofNum = [endofNum,length(y)];
        end

        assert(length(startofNum) == length(endofNum),...
            'shadederrorbar:splitNumbers:length',...
            'startofNum (%d) and endofNum (%d) are different in length',...
            length(startofNum),length(endofNum))

        yC = cell(1,length(startofNum));
        errC     = yC;
        xC       = yC;

        for i = 1:length(startofNum)
            ind = startofNum(i):endofNum(i);
            yC{i} = y(ind);
            errC{i} = err(ind);
            xC{i} = x(ind);
        end
    end

else
    yC = {y};
    errC = {err};
    xC = {x};
end


end
