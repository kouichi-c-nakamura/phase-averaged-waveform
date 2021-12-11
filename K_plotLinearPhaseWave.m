function [h, N, outdata] = K_plotLinearPhaseWave(varargin)
% K_plotLinearPhaseWave create a linear histogram out of phase values in
% radians. It works on single sample and group data as well.
%
% [h, N, outdata] = K_plotLinearPhaseWave(S)
% [h, N, outdata] = K_plotLinearPhaseWave(axh, _____)
% [h, N, outdata] = K_plotLinearPhaseWave(_____, 'Param', Value, ...)
%
%
% INPUT ARGUMENTS
%
% S            Scalar (single data) or non-scalar (groupdata) structure 
%              with the fields
%                 
%                     'binmean'
%                     'binstd'
%                     'binsem'
%                     'axrad'
%                     'meanvec'
%                     'bootstrap'
%
% axh         [] | an Axes object | two Axes objects
%             Axes handle. This can be empty, scalar or two-element array.
%             For scalar axh, this function will delete axh and draw the
%             main and sub axes in place of axh, so that The outmost
%             rectangle of the main and sub axes tallies with Position of
%             axh. If axh has two elements, then axh(1) will be the main
%             axes and axh(2) will be the sub axes.
%
% OPTIONAL PARAMETER/VALUE PAIRS
%
% 'Color'     ColorSpec
%
% 'XGrid'     'on' (default) | 'off' (1 or 0)
%
% 'PlotType'  'line' | 'bar' (default) | 'surface' |'none'
%
% 'ErrorBar'  'none' (default) | 'std' | 'sem'
%
% 'AnalysisMode' 
%             'auto' (default) | 'one' | 'group'
%             When S is scalar you can decide which type of analysis you
%             want to perform.
%
% 'Title'     string for title
%
% 'YUnit'     'mv2mv' (default) | 'mv2microv'
%             'mv2microv' will 1000-fold values for Y axis and use
%             microvolt(µV) as unit.
%
% OUTPUT ARGUMENTS
%
% h           structure of handles
%               h.main       main axes
%               h.sub        sub axes for sign curve
%               h.titlepane  axes for title
%
% N           Sample number
%
% outdata     Structure with fields:
%               outdata.x
%               outdata.y
%               outdata.yerr
%
% Supporting group data and shaded error bars, 21/06/2013
%
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% kouichi.c.nakamura@gmail.com
% 24-May-2017 13:13:54
% 
% See also
% K_plotCircPhaseWave_one, K_plotCircPhaseWave_group
% K_plotLinearPhaseHist, K_plotLinearPhase, K_PhaseWave
% K_PhaseHist, K_plotCircPhaseHist_one
% K_plotCircularSingle, K_plotCircularGroup % what's the difference among them?


[axrad,binmean,binerr,errorBar,axh,ColorSpec,xgrid,plottype,titlestr,...
    btstrp,circsh,meanvec,analysismode,yunit] = local_parser(varargin{:});

%%

switch analysismode
    case 'one'

        N = 1;
        x = axrad;
        y = binmean;
        yerr = binerr;
    
    case 'group'
            
        N = size(binmean,2);
        
        assert(all(arrayfun(@(k) all(axrad(:,1) == axrad(:,k)),1:N)),...
            'Assuming axrad is identical among samples (mandatory to compute group mean).')
        
        if isempty(axrad)
            h = [];
            outdata = [];
            disp('No data to plot.')
            return
        end
        x = axrad(:,1);
        
        switch plottype
            case {'line','bar','none'}
                y = nanmean(binmean,2);
                switch errorBar
                    case 'std'
                        yerr = nanstd(binmean,0,2);
                    case 'sem'
                        binstd = nanstd(binmean,0,2);
                        yerr = binstd/sqrt(N);
                        
                    case 'none'
                        yerr = [];
                end
            case 'surface'
                y = binmean;
                yerr = [];
        end
        
end

if strcmpi(yunit,'mv2microv')
    y = y.*1000;
    yerr = yerr.*1000;
end

if strcmp(plottype,'none')
    axh = [];
    ylabelstr = '';
    [h,outdata] = plotLinearPhase(axh,x,y,yerr,analysismode,ColorSpec,xgrid,plottype,...
        titlestr,ylabelstr,errorBar);
    
else

    if isempty(axh)
        figure;
        axh = axes;
    end
    
    switch plottype
        case 'surface'
            ylabelstr = '';
        otherwise
            if strcmpi(yunit,'mv2microv')
                ylabelstr = 'Modulation (µV)';
            else
                ylabelstr = 'Modulation (mV)';
            end
    end
    
    [h,outdata] = plotLinearPhase(axh,x,y,yerr,analysismode,ColorSpec,xgrid,plottype,...
        titlestr,ylabelstr,errorBar);
    
    if strcmpi(yunit,'mv2microv') && strcmpi(plottype,'surface')
           h.colorbar.Label.String = regexprep(h.colorbar.Label.String,'\(mV\)',...
               sprintf('(%sV)',char(181)));
    end
    %% Add text
    
    if ~isempty(btstrp) && ~isempty(meanvec)
        h.txt = local_placetext(N,btstrp,circsh,meanvec,h.main.axh,analysismode); %TODO
    end
    
end
    

end

%--------------------------------------------------------------------------

function [axrad,binmean,binerr,errorBar,axh,ColorSpec,xgrid,plottype,titlestr,...
    btstrp,circsh,meanvec,analysismode,yunit] = local_parser(varargin)


% initialization

axh = [];

% vfaxhrad = @(x) isvector(x) && isnumeric(x) || ...
%     iscell(x) && ...
%     all(cellfun(@(y) isnumeric(y) && iscolumn(y) || isempty(y), x)) ||...
%     all(cellfun(@(y) isnumeric(y) && isrow(y) || isempty(y) , x));
%NOTE cell contents must be all columns or all rows

if ishandle(varargin{1})
    % [h,N] = K_plotLinearPhaseWave(axh, _____)
    
    axh = varargin{1};
    assert(numel(axh) <= 2 && all(isgraphics(axh,'axes')))
    
    varargin = varargin(2:end);
   
end

p = inputParser;

p.addRequired('S',@(x) isstruct(x) && isequal(sort(fieldnames(x)),sort({...
    'binmean';...
    'binstd';...
    'binsem';...
    'axrad';...
    'meanvec';...
    'circshift';...
    'bootstrap'})));

p.addParameter('Color','b',@iscolorspec);

p.addParameter('XGrid','on',@(x) ~isempty(x) &&...
    ischar(x) &&...
    isrow(x) &&...
    ismember(x, {'on','off'}));

p.addParameter('ErrorBar','sem',@(x) ismember(x,{'std','sem','none'}));

p.addParameter('PlotType','line',@(x) ~isempty(x) &&...
    ischar(x) &&...
    isrow(x) &&...
    ismember(x, {'line','bar','surface','none'}));

p.addParameter('Title','',@(x) ~isempty(x) &&...
    ischar(x) &&...
    isrow(x));

p.addParameter('AnalysisMode','auto',@(x) ismember(x,{'one','group','auto'}));

p.addParameter('YUnit','mV2mV',@(x) ismember(lower(x),{'mv2mv','mv2microv'}));

p.parse(varargin{:});

S         = p.Results.S;
ColorSpec = p.Results.Color;
errorBar  = p.Results.ErrorBar;
analysismode = p.Results.AnalysisMode;

xgrid    = p.Results.XGrid;
plottype = p.Results.PlotType;
titlestr = p.Results.Title;
yunit    = lower(p.Results.YUnit);

axrad   = [S(:).axrad]; % S(k).axrad is a column vector
binmean = [S(:).binmean];
btstrp  = [S(:).bootstrap];
circsh  = [S(:).circshift];
meanvec = [S(:).meanvec];

switch errorBar
    case 'std'
        binerr = [S(:).binstd];
    case 'sem'
        binerr = [S(:).binsem];
    case 'none'
        binerr = [];
end

if isscalar(S)
    if isequal(analysismode,'auto')
        analysismode = 'one';
    else
        % user's analysismode is chosen
    end
else
    analysismode = 'group';
end

end

%--------------------------------------------------------------------------

function txth = local_placetext(N,btstrp,circsh,meanvec,mainaxh,analysismode)


switch analysismode
    case 'one'
        if (isempty(btstrp) && isempty(circsh)) || ...
                (isempty(btstrp.p_lessthan) && isempty(circsh.p_lessthan))
            pstr = '';
            
        else
             if ~isempty(circsh.p_lessthan)
                pstr = sprintf('Circshift test: p < %f\n', ...
                    circsh.p_lessthan);
                
             elseif ~isempty(btstrp.p_lessthan)
            
                pstr = sprintf('Bootstrap test: p < %f\n', ...
                    btstrp.p_lessthan);

            end
        end
        
        thestr = sprintf(['%s',...
            'Phase of mean vector: %.3f%s\n',...
            'Length of mean vector: %.3e mV\n'], ...
            pstr,...
            meanvec.degree,char(176),...
            meanvec.length);
        
    case 'group'
        
        % circmean of phase values and mean length
        
        samplevecrad = [meanvec(:).radian]';
        groupradmean = circ_mean(samplevecrad);
        groupradstd = circ_std(samplevecrad);
        
        sampleveclen = [meanvec(:).length];
        groupveclenmean = mean(sampleveclen);
        groupveclenstd = std(sampleveclen);
        
        if isempty([circsh(:).p_lessthan]) && isempty([circsh(:).p_lessthan])
            pstr = '';
            
        else
            if ~isempty([circsh(:).p_lessthan])
                p = [circsh(:).p_lessthan]';
                
                n1 = nnz(p < 0.001);
                n2 = nnz(p < 0.01);
                n3 = nnz(p < 0.05);
                n4 = nnz(0>= 0.05);
                
                sigN = n3;
                
                pstr = sprintf('Circshift test (p < 0.05): %d of %d\n', ...
                    sigN,N);
                
            elseif ~isempty([circsh(:).p_lessthan])
                
                p = [btstrp(:).p_lessthan]';
                
                n1 = nnz(p < 0.001);
                n2 = nnz(p < 0.01);
                n3 = nnz(p < 0.05);
                n4 = nnz(0>= 0.05);
                
                sigN = n3;
                
                pstr = sprintf('Bootstrap test (p < 0.05): %d of %d\n', ...
                    sigN,N);

            end
        end
        
        thestr = sprintf(['n = %d\n',...
            '%s',...
            'Group cicrular mean %s std: %.3f %s %.3f%s\n', ...
            'Vector length mean %s std: %3f %s %3f mV\n'], ...
            N,...
            pstr,...
            char(177),rad2deg(groupradmean),char(177),rad2deg(groupradstd),char(176),...
            char(177),groupveclenmean,char(177),groupveclenstd);
end

axes(mainaxh);

txth = text(0.95, 0.95, thestr,...
    'HorizontalAlignment', 'right', 'VerticalAlignment','top',...
    'Units', 'normalized','Tag','Annotation Text');

end
