function [h, outdata] = plotLinearPhase(axh,x,y,yerr,analysismode,ColorSpec,xgrid,...
    plottype,titlestr,ylabelstr,errorBar)
%
% plotLinearPhase is a low-level plotting subroutine for
% K_plotLinearPhaseHist and K_plotLinearPhaseWave to draw a plot for a
% linear phase histogram
%
% [h, outdata] = pvt_plotLinearPhase(axh,x,y,yerr,analysismode,ColorSpec,xgrid,...
%     plottype,titlestr,ylabelstr,errorBar)
%
% INPUT ARGUMENTS
% axh         [] | an Axes object | two Axes objects
%             Axes handle. This can be empty, scalar or two-element array.
%             For scalar axh, this function will delete axh and draw the
%             main and sub axes in place of axh, so that The outmost
%             rectangle of the main and sub axes tallies with Position of
%             axh. If axh has two elements, then axh(1) will be the main
%             axes and axh(2) will be the sub axes.
%
% x,y,yerr    Vectors of the same length for x, y values and y error
%             values.
%
%             Array format of y (rows for phase, columns for samples) is
%             accepted only for plottype = 'surface', as long as size(y,1)
%             == length(x)
%
% analysismode
%             'one'  | 'group'
%
% ColorSpec   colorSpec
%
% xgrid       'on' | 'off'
%             Specifies XGrid property of Axes
%
% plottype    'line' | 'bar' | 'surface' %TODO |'none'
%
% titlestr    String for title
%
% ylabelstr   String for ylabel
%
% errorBar    'std' | 'sem' | 'none'
%
%
% OUTPUT ARGUMENTS
% h           Structure containing graphic handles
%                 h.fig
%                 h.main
%                 h.sub
%                 h.titlepane
%                 h.colorbar
%
% outdata     Structure with the following fields
%                 outdata.x
%                 outdata.y
%                 outdata.yerr
%
% See also
% K_plotLinearPhaseHist, K_plotLinearPhaseWave


%% Parse inputs
narginchk(11,11);

p = inputParser;
p.addRequired('axh', @(X) isempty(X) || ...
    numel(X) <=2 && all(isgraphics(X,'axes')));
p.addRequired('x',   @(X) isvector(X) && isnumeric(X));
p.addRequired('y',   @(X) isnumeric(X) && ...
    (isvector(X) && length(X) == length(x)) || ...
    size(X,1) == length(x));
p.addRequired('yerr',@(X) isempty(X) || isvector(X) && isnumeric(X) &&...
    length(X) == length(x));
p.addRequired('analysismode',@(X) ismember(lower(X),{'one','group'}));
p.addRequired('ColorSpec',@(X) iscolorspec(X));
p.addRequired('xgrid',    @(X) ismember(lower(X),{'on','off'}));
p.addRequired('plottype', @(X) ismember(lower(X),{'line','bar','surface','none'}));
p.addRequired('titlestr', @(X) ischar(X) && isempty(X) || isrow(X));
p.addRequired('ylabelstr',@(X) ischar(X) && isempty(X) || isrow(X));
p.addRequired('errorBar', @(X) ismember(lower(X),{'std','sem','none'}));

p.parse(axh,x,y,yerr,analysismode,ColorSpec,xgrid,plottype,titlestr,...
    ylabelstr,errorBar);

switch plottype
    case {'line','bar','none'}
        assert(isvector(y))
end

if isrow(x)
    x = x';
end

if isrow(y)
    y = y';
end

if isrow(yerr)
    yerr = yerr';
end

rad2ang   = @(X) X/pi*180;
ang2rad   = @(X) X * pi /180;
twocycles = @(X) circshift(repmat(X, 2, 1), round(size(X, 1)/2)); % two cycles
plusalpha = @(X) [X(end,:);X(:,:);X(1,:)];

%% Prepare axes

left1   = 0.14;
bottom1 = 0.12;
width1  = 0.77;
height1 = 0.68;

left2   = left1;
bottom2 = 0.84;
width2  = width1;
height2 = 0.15;

gaph = bottom2 - height1- bottom1;

if isempty(axh) && ~ strcmp(plottype,'none')
    fig =figure;
    main.axh = subplot('Position', [left1, bottom1, width1, height1],'Tag','Main'); hold on;
    sub.axh  = subplot('Position', [left2, bottom2, width2, height2],'Tag','Sub'); hold on;
    
elseif strcmp(plottype,'none')
    fig = [];
    main = [];
    sub = [];
elseif numel(axh) == 1
    %subdivide the axh
    
    origUnits = get(axh,'Units');
    set(axh,'Units','normalized');
    
    origpos = get(axh,'Position');
    axes(axh);
    fig = gcf;
    delete(axh);
    
    pos1 = [origpos(1),...
        origpos(2), ...
        origpos(3),...
        origpos(4) *height1/(gaph+height1+height2)];
    
    pos2 = [origpos(1),...
        origpos(2) + origpos(4) *(height1+gaph)/(gaph+height1+height2), ...
        origpos(3),...
        origpos(4) *height2/(gaph+height1+height2)];
    
    %NOTE on R2017a prerelease, You need to specifiy 'Units', 'normalized'
    main.axh = axes('Units','normalized','Position', pos1,'Tag','Main');
    hold on;
    set(main.axh,'Units',origUnits);
    sub.axh = axes('Units','normalized','Position', pos2,'Tag','Sub');
    hold on;
    set(sub.axh,'Units',origUnits);
    
elseif  numel(axh) == 2
    
    main.axh = axh(1);
    sub.axh  = axh(2);
    axes(axh(1))
    fig = gcf;
    
end

%% Prepare data
xx = round(rad2ang(unwrap(plusalpha(twocycles(x))))); % 2 cylces + alpha
% x = round(rad2ang(unwrap([x(end);x;x(1)]))); % 2 cylces + alpha

switch plottype
    case {'line','bar','none'}
        yy = plusalpha(twocycles(y));  % 2 cylces + alpha
        % y = [y(end); y; y(1)]; % 2 cylces + alpha
    case 'surface'
        yy = plusalpha(twocycles(y)); %twocycles is working
end

% if strcmpi(analysismode,'group')
if isempty(yerr) || strcmpi(errorBar,'none')
    yerryerr = [];
else
    yerryerr = plusalpha(twocycles(yerr)); %TODO *100
end
% end


%% Plot
cbh = [];
switch plottype
    case 'line'
        main.lh1(1) = plot(main.axh, xx, yy, 'LineWidth', 1,...
            'Color', ColorSpec,'Tag','Line Histogram');
        
        axes(main.axh);
        if  isempty(yerryerr) || strcmpi(errorBar,'none')
            main.eh1 = gobjects(1);%TODO need to check
        else
            main.eh1 = shadederrorbar(xx,yy,yerryerr,ColorSpec);
            set(main.eh1,'Tag','Error Range');
            
        end
        
    case 'bar'
        % positive Y
        main.lh1(1) = bar(main.axh, xx(yy>=0), yy(yy>=0),'BarWidth', 1, ...
            'FaceColor', ColorSpec, 'Tag','Bar Histogram Positive');
        
        % negative Y
        if any(yy<0)
            main.lh1(2,1) = bar(main.axh, xx(yy<0), yy(yy<0),'BarWidth', 1, ...
                'FaceColor', [0.5 0.5 0.5] ,'Tag','Bar Histogram Negative');
        end
        
        axes(main.axh);
        if  isempty(yerryerr) || strcmpi(errorBar,'none')
            main.eh1 = gobjects(1);%TODO need to check
        else
            main.eh1 = errorbar(xx,yy,yerryerr);
            set(main.eh1 ,'Color','k','LineStyle','none',...
                'Tag','Error Bar');
            
            uistack(main.eh1, 'bottom');
            
        end
        
        if ispc
            f = gcf;
            f.Renderer = 'painter'; % work around for drawing bug
        end
    case 'surface'
        
        n = size(yy,2);
        
        main.lh1(1) = surface(main.axh, ...
            repmat(xx, 1, n+2),...
            repmat([0,0.5:n-0.5,n],size(yy,1),1),...
            [yy(:,1),yy,yy(:,end)],'Tag','Surface');
        
        %NOTE In order to provide "width" to yy(1,:) and yy(end,:), two
        % bands of 0.5 width have been added at the bottom and top. This
        % solution contains a little bit of misleading, but may be more
        % intuitive than simplely ploting yy, in which case the bottom and
        % top channels do not represent any space.
        
        shading(main.axh,'interp')
        
        if exist('redblue','file') == 2
            colormap(main.axh,redblue);
        else
            colormap(main.axh,parula);
        end
        
        %TODO colorbar?
        main_units = main.axh.Units;
        main.axh.Units = 'normalized';
        pos = main.axh.Position;
        
        cbh = colorbar(main.axh,'eastoutside');
        cbh.Tag = 'PhaseWave Colorbar';
        main.axh.Position = pos; % this works for 'Units','normalized'
        main.axh.Units = main_units;
        
        ylabel(cbh,'Modulation (mV)')
        cbh.TickDirection = 'out';
        
        cbh.Position(3) = cbh.Position(3)/3;
        cbh.Position(1) = pos(1)+pos(3)+0.01;
        
        
        ylim(main.axh,'auto')
        ylabel(main.axh,'')
        
        main.axh.YTick = 0:n;
        main.axh.YTickLabel = [];
        for i = 1:n
            text(main.axh,-10,i-0.5,sprintf('%d',i),'Tag','YTickLabel Text',...
                'HorizontalAlignment','right')
        end
        
    case 'none'
        fig = [];
        main = [];
        sub = [];
        titleh = [];
    otherwise
        error('unexpected plottype %s', plottype);
end

if ~strcmp(plottype,'none')
    ylimhist = ylim(main.axh);
    
    
    xtic = cellfun(@(x) [num2str(x), char(176)], num2cell(0:90:720), 'UniformOutput', false); %  char(176) == degree symbol
    set(main.axh, 'TickDir', 'out', 'Box', 'off',...
        'XTick', 0:90:720, 'XTickLabel', xtic,...
        'XLim', [0, 720]);
    
    if strcmp(xgrid,'on')
        axes(main.axh);
        set(main.axh,'XGrid','on','GridColor',[0.5 0.5 0.5]);
    end
    
    ylim(ylimhist);
    clear ylimhist
    
    xlabel(main.axh,'Phase');
    ylabel(main.axh, ylabelstr);
    
    
    %% Add text
    
    % txth = local_placetext(N,radians,rayleighECDF,main); %TODO
    
    %% subplot for cosine curve
    
    x2 = linspace(0, 720, 1440);
    y2 = cos(ang2rad(x2));
    sub.l1 =plot(sub.axh, x2, y2, 'Color' , 'k', 'LineWidth', 2,...
        'Tag','Cosine Curve');
    set(sub.axh, 'Visible', 'off', 'XLim', [0, 720], 'YLim', [-1.1, 1.1]);
    
    if ~isempty(titlestr)
        
        titleh= title(sub.axh,titlestr,'Interpreter','none');
        titleh.Visible = 'on';
        
    else
        titleh = [];
    end
    
end

h.fig  = fig;
h.main = main;
h.sub  = sub;
h.titlepane = titleh;
h.colorbar = cbh;

outdata.x = xx;
outdata.y = yy;
outdata.yerr = yerryerr;



end