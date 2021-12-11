function h = K_plotCircPhaseHist_one(varargin)
% K_plotCircPhaseHist_one create the following out of phase values in radians.
% It works on single sample only.
%
%    1. a circular histogram (rose plot) 
%    2. a single vector representing circular mean of phase values
%    3. small circles that indicate each phase value 
%
%   h = K_plotCircPhaseHist_one(radians)
%   h = K_plotCircPhaseHist_one(radians,nbin)
%   h = K_plotCircPhaseHist_one(axh, ______)
%   h = K_plotCircPhaseHist_one(_____, 'Param', Value, ...)
% 
% INPUT ARGUMENTS
% 
% radians          A vector of phase values in radian
%
% nbin             The number of histogram bins
%
% axh              Axes handle (optional)
% 
% OPTIONAL PARAMETER/VALUE PAIRS
% 'Color'          ColorSpec 
%
% 'ZeroPos'        String. 
%                  'top' (default) | 'right' | 'left' |'bottom'
%
% 'Direction'      String. 
%                  'clockwise' (default) | 'anti'     
%
% 'MeanVector'     'on' (default) | 'off'
%
% 'SmallCircle'    'on' (default) | 'off'
%
% 'SmallCircleSize'
%                   Positive number. Default is 10. 
%
% 'Histogram'       'on' (default) | 'off'
%
% 'HistLimPercent'  Option for nomalized histogram rather than counts. You
%                   can specifiy the outerlimit of the circular histogram in
%                   percent. If HistLimPercent is 0, the outerlimit is
%                   automatically set.If Histbin is not provided,
%                   HistLimPercent will be ignored.
%
% 'rayleighECDF'    Test results (P value) of Rayleigh's test (cirt_rtest) in
%                   vector. The length must be the same as that of the sample
%                   number of radians (1 if radians is numeric vector, 
%                   numel(radians) if cell vector).
%
% 'RadiusForHistScale'
%                   [0.60, 1.10] (default) | [inner, outer]
%                   Determines the positions of texts indicating scales
%
% 'RadiusForAngles'
%                   [1.2, 1.2] (default) | [r0, r90]
%                   Determines the positions of texts indicating angles.
%                   r0 is for 0 and 180 degrees, whereas r90 is for 90 and 270.
%
% 'RadiusForScales'
%                   [0.65,1.2] (default) | [inner, outer] 
%                   Determines the positions of texts indicating scales
%
% 'RadianForScales'
%                   pi*3/8 (default) | scalar
%                   Phase of the positions of texts indicating scales in radian
%
% OUTPUT ARGUMENTS
%
% h                 A structure of handles for graphic objects.
%
%
% See Also 
% K_PhaseHist
% K_plotLinearPhaseHist 
% K_plotLinearPhaseHist_S  (to directly use output of K_PhaseHist as an input argument)
% K_PhaseHist_histlabel    (add summary text to the plots made by K_PhaseHist)
% K_plotCircPhaseHist_one, K_plotCircPhaseHist_group (for circular plots)
% K_plotColorPhaseHist     (heatmap representation of phase coupling)
% K_ECDFforRayleigh        (ECDF-based correction for Rayleigh's uniformity test)
% K_PhaseHist_test         (UnitTest)
% pvt_K_plotCircPhaseHist_parseInputs
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% kouichi.c.nakamura@gmail.com
% 05-May-2017 09:31:33



%% parse input arguments

[axh, radians, nbin, zeropos, plotdir, meanvector, smallcircle, smallcirclesize,...
    histogram, histlimpercent, ColorSpec, rayleighECDF, radiusForHistScale, ...
    radiusForAngles, radiusForScales, radianForScales] = local_parser(varargin{:});

%% Job

if isempty(radians)
    cmean = NaN;
else
    cmean = circ_mean(radians);
end
veclen = circ_r(radians);

% convert a radian to a complex number for polar coodinate
rad2cmp = @(x) exp(1i * x);

% zz = rad2cmp(linspace(0, 2*pi, 3600));

z1 = rad2cmp(cmean);
z2 = rad2cmp(radians);

%%

h = plotCircleGrids(axh,'RadiusForAngles',radiusForAngles,...
    'RadiusForScales',radiusForScales,'RadianForScales',radianForScales);


%% histogram

if ~isempty(nbin)
    fig2 = figure('Visible', 'off');
    histline = rose(radians, nbin);
    lim = max(get(gca, 'XLim'));
    histx = get(histline, 'XData'); 
    histy = get(histline, 'YData');
    % four succesive values form a triangle starting (the first of four) and ending (the last of four) at the origin 
    
    delete(fig2);

    axes(h.axh);
    
    if ~isempty(histlimpercent)
       
        if histlimpercent == 0
            edges = linspace(-pi, pi, nbin +1);
            bincount = histc(radians, edges);
            bincount = [bincount(1:end-2);...
                (bincount(end-1) +bincount(end))];
           histlimpercent = ceil( max(bincount)/numel(radians) *110);                      
        end
        
        
        num = numel(radians);
        % defaultouterlimpercent = lim/num *100;
        
        scale = 1 /num *100 /histlimpercent;
        
        X = histx * scale;
        Y = histy * scale;
        
        histstr1 = [num2str(histlimpercent/2), '%'];
        histstr2 = [num2str(histlimpercent), '%'];
        
    else
        X = histx/lim;
        Y = histy/lim;   
        
        histstr1 = num2str(lim/2);
        histstr2 = num2str(lim);
    end
    
    h.histline = plot(h.axh, X, Y, 'Color', 'k',...
        'Tag','Petal Edge');
    
    % four succesive values of X and Y form a triangle starting (the first
    % of four) and ending (the last of four) at the origin
    
    assert(rem(length(X), 4) == 0)
    h.histpatch = zeros(length(X)/4, 1);

    for i = 1:length(X)/4 
        h.histpatch(i) = fill(X(i*4-3:i*4), Y(i*4-3:i*4), [0.7, 0.7, 0.7],...
            'Tag','Petal');
    end
    hg = hggroup('Tag','Rose Histogram');
    set(h.histpatch,'Parent',hg);
    set(h.histline,'Parent',hg);
    
    axis equal
    
    tz2 =rad2cmp(pi*1/8);

    h.txt_scalehist(1) = text(real(tz2)*radiusForHistScale(1), ...
        imag(tz2)*radiusForHistScale(1), histstr1,'Tag','Histogram Scale Full');
    
    h.txt_scalehist(2) = text(real(tz2)*radiusForHistScale(2), ...
        imag(tz2)*radiusForHistScale(2), histstr2,'Tag','Histogram Scale Half');
    
    set(h.txt_scalehist, 'HorizontalAlignment', 'center');
    
    if strcmp(histogram,'off')
        set(hg,'Visible','off');
    end 

end

%% draw small open circles for radians on perimeter
h.circles = plot(h.axh, real(z2), imag(z2),...
    'Marker', 'o', 'MarkerSize', smallcirclesize,'Linestyle', 'none', ...
    'Color', ColorSpec, 'Tag','Small Circles on Perimeter');

if strcmp(smallcircle,'off')
   set(h.circles,'Visible','off'); 
end
    

%% draw a vector and X mark on perimeter for circular mean

h.vector = plot(h.axh, [0 real(z1)*veclen], [0 imag(z1).*veclen],...
    'Color', ColorSpec , 'LineWidth', 5,'Tag','Vector');
axis equal;

if strcmp(meanvector,'off')
   set(h.vector,'Visible','off'); 
end

h.mark = plot(h.axh, real(z1), imag(z1),...
    'Marker', '*', 'MarkerSize', smallcirclesize,'Linestyle', 'none', ...
    'Color', ColorSpec, 'Tag','Mark on Perimeter per Group');
hold off


%% add text labels

hold off

txth = local_placetext(radians,rayleighECDF,h.axh);

%% rotation
switch plotdir
    case 'clockwise'
        switch zeropos
            case 'right'
                view(0, -90)
            case 'left'
                view(180, -90)
            case 'top' % default
                view(90, -90)
            case 'bottom'
                view(270, -90)
        end
    case 'anti'
        switch zeropos
            case 'right'
                view(0, 90)
            case 'left'
                view(180, 90)
            case 'top'
                view(270, 90)
            case 'bottom'
                view(90, 90)
        end
end
end

%--------------------------------------------------------------------------
function [axh, radians, nbin, zeropos, plotdir, meanvector, smallcircle, ...
    smallcirclesize, histogram, histlimpercent, ColorSpec, rayleighECDF,...
    radiusForHistScale, radiusForAngles, radiusForScales, radianForScales] ...
    = local_parser(varargin)

vfrad = @(x) isreal(x) && isvector(x);

% initialize
axh = [];
ColorSpec  = 'b';
nbin = 36;
zeropos = 'top';
plotdir = 'clockwise';
histlimpercent = [];
rayleighECDF = [];
meanvector = 'on';
smallcircle = 'on';
histogram = 'on';
smallcirclesize = 10;

radiusForHistScale = [0.60, 1.10];
radiusForAngles = [1.2,1.2];
radiusForScales = [0.65,1.2];
radianForScales =  pi*3/8;

% vfrad = @(x) isnumeric(x) &&...
%     isvector(x);

if nargin == 1
    % h = K_plotCircPhaseHist_one(radians)

    p = inputParser;
    
    addRequired(p, 'radians', vfrad);
    parse(p, varargin{:});
    
    radians = varargin{1};
    
elseif nargin >= 2
    % h = K_plotCircPhaseHist_one(radians, nbin) %TODO
    % h = K_plotCircPhaseHist_one(axh, ______)
    % h = K_plotCircPhaseHist_one(_____, 'Param', Value, ...)
    
    if ishandle(varargin{1})
        % h = K_plotCircPhaseHist_one(axh, ______)

        axh = varargin{1};
        radians = varargin{2};
        
        if length(varargin) >=3 && ~ischar(varargin{3})
            % h = K_plotCircPhaseHist_one(axh,radians,nbin)

            nbin = varargin{3};
            PNVStart = 4;
        else
            % h = K_plotCircPhaseHist_one(axh,radians)

            PNVStart = 3;
        end
        
    else
        radians = varargin{1};
        
        if ~ischar(varargin{2})
            nbin = varargin{2};
            PNVStart = 3;
        else
            PNVStart = 2;
        end
        
    end
    
    PNV = varargin(PNVStart:end);
    
    p = inputParser;
    
    p.addRequired('radians', vfrad);
    
    vf2 = @(x) isempty(x) ||...
        isnumeric(x) && isscalar(x) &&...
        fix(x) == x && x >=0;
    p.addRequired('nbin', vf2);
    
    p.addParameter('Color',ColorSpec,@iscolorspec);
    
    p.addParameter('MeanVector',meanvector,@(x) ~isempty(x) && ischar(x) && isrow(x) ...
        && ismember(lower(x),{'on','off'}));

    p.addParameter('ZeroPos',zeropos,@(x) ~isempty(x) && ischar(x) ...
        && isrow(x) && ismember(lower(x),{'top','left','right','bottom'}));
    
    p.addParameter('SmallCircle',smallcircle,@(x) ~isempty(x) && ischar(x) && isrow(x) ...
        && ismember(lower(x),{'on','off'}));

    p.addParameter('SmallCircleSize',smallcirclesize,@(x) isscalar(x) && isreal(x) && x > 0)
    
    p.addParameter('Direction',plotdir,@(x) ~isempty(x) && ischar(x) && isrow(x) ...
        && ismember(lower(x),{'clockwise','anti'}));
    
    p.addParameter('Histogram',histogram,@(x) ~isempty(x) && ischar(x) && isrow(x) ...
        && ismember(lower(x),{'on','off'}));

    p.addParameter('HistLimPercent',histlimpercent,@(x) isempty(x) ...
        || isnumeric(x) && isscalar(x) && x >= 0);

    p.addParameter('rayleighECDF',rayleighECDF,@(x) isvector(x) && isreal(x));
    
    
    p.addParameter('RadiusForHistScale',radiusForHistScale, @(x) isrow(x) && numel(x) == 2 ...
        && isnumeric(x) && all(x >=0) && x(1) < x(2))
    
    p.addParameter('RadiusForAngles', radiusForAngles, @(x) isrow(x) && numel(x) == 2 ...
        && isnumeric(x) && all(x >=0))
    
    p.addParameter('RadiusForScales', radiusForScales, @(x) isrow(x) && numel(x) == 2 ...
        && isnumeric(x) && all(x >=0) && x(1) < x(2))
    
    p.addParameter('RadianForScales',radianForScales, @(x) isscalar(x) && isreal(x))

    %NOTE see also plotCircleGrids
    
    
    
    p.parse(radians, nbin, PNV{:});
    
    if isempty(nbin) || nbin == 0
        nbin = 36; % default value
    end
    
    ColorSpec = p.Results.Color;
    meanvector = lower(p.Results.MeanVector);
    histogram = lower(p.Results.Histogram);
    smallcircle = lower(p.Results.SmallCircle);
    histlimpercent = p.Results.HistLimPercent;
    zeropos = lower(p.Results.ZeroPos);
    smallcirclesize = p.Results.SmallCircleSize;
    plotdir = lower(p.Results.Direction);
    rayleighECDF = p.Results.rayleighECDF;
    
    radiusForHistScale = p.Results.RadiusForHistScale;
    radiusForAngles = p.Results.RadiusForAngles;
    radiusForScales = p.Results.RadiusForScales;
    radianForScales = p.Results.RadianForScales;

end


end

%--------------------------------------------------------------------------

function txth = local_placetext(radians,rayleighECDF,axh)


if isempty(rayleighECDF)
    raylstr = '';
    
else
    raylstr = sprintf('Rayleigh'' test: p = %.3e\n', ...
        rayleighECDF);
end

thestr = sprintf(['%s',...
    'Cicrular mean %s std: %.3f %s %.3f%s\n',...
    'Vector length: %.3f\n'], ...
    raylstr,...
    char(177),rad2deg(circ_mean(radians)),char(177),rad2deg(circ_std(radians)),char(176),...
    circ_r(radians));
    
axes(axh);

txth = text(1.3, -0.15, thestr,...
    'HorizontalAlignment', 'right', 'VerticalAlignment','bottom',...
    'Units', 'normalized','Tag','Annotation Text',...
    'FontSize',9);

end
%--------------------------------------------------------------------------

function eid = eid(varargin)
% eid = eid()
% eid = eid(string)
% Local function that generates error id that begins with K:
%
%
% input argument
% str (Optional) string in char type (row vector)
%
% output argument
% eid an error id composed of 'K:(functionname):str'

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
            
    
 



 


