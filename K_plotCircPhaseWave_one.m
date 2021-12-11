function h = K_plotCircPhaseWave_one(varargin)
% K_plotCircPhaseWave_one create the following out of phase values in radians.
% It works on single sample only.
%
%    1. a circular representation of average potentials (patch plot) 
%    2. a single vector representing the synthetic mean of phase-potential
%    vectors
%
%   h = K_plotCircPhaseWave_one(S) 
%   h = K_plotCircPhaseWave_one(axh, ______)
%   h = K_plotCircPhaseWave_one(_____, 'Param', Value, ...)
% 
% INPUT ARGUMENTS
% 
% S              Structure with the fields
%                 
%                     'binmean'
%                     'binstd'
%                     'binsem'
%                     'axrad'
%                     'meanvec'
%                     'bootstrap'
%                     'circshift'
%
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
% 'ShowMeanVector' true (default) | false
%
% 'ShowPatch'      true (default) | false
%
% 'ShowBootstrapVectors'
%                   true | false (default)
%                   If true plots sample vectors for each iteration of
%                   Bootstraping as a small dot. They are usually very
%                   short and not quite visible.
%
% 'ErrorRange'     'none' (default) | 'std' | 'sem'
%
% 'Radius'         'auto' (default) | scalar number
%                  Specify the radius of the circular axes.
%
% 'YUnit'        'mv2mv' (default) | 'mv2microv'
%                'mv2microv' will 1000-fold values for Y axis and use
%                microvolt(µV) as unit.
%
% OUTPUT ARGUMENTS
%
% h                 A structure of handles for graphic objects.
%
%
% See Also 
% K_plotCircPhaseWave_group, K_plotLinearPhaseWave,
% K_plotCircPhaseHist_group, K_plotCircPhaseHist_one, K_PhaseHist,
% K_PhaseHist_test, pvt_K_plotCircPhaseHist_parseInputs

%% parse input arguments

[axrad,binmean,binerr,axh,zeropos, plotdir, showmeanvector, showpatch,...
    showbtstrpvecs,ColorSpec,radius,meanvec,btstrp,titlestr,yunit] = ...
    local_parser(varargin{:});

%% Job

% convert a radian to a complex number for polar coodinate
plusalpha = @(X) [X(end);X;X(1)]; % to complete the circle

theta = plusalpha(axrad);
rho   = plusalpha(binmean);

if isequal(radius,'auto')
    if ~isempty(binerr)
        rho_ = plusalpha(binerr);
        RHO = arrayfun(@(x,y) x+sign(x)*y, rho,rho_);
        maxradius = max(xlimautoval([0,max(abs(RHO))]));
        
    else
        maxradius = max(xlimautoval([0,max(abs(rho))]));
    end
else
    maxradius = radius;
end

h = plotCircleGrids(axh,maxradius);

if strcmpi(yunit,'mv2microv') %TODO
        
    h.txt_scalevec(1).String = [num2str(str2double(h.txt_scalevec(1).String)*1000), ...
        sprintf(' %sV',char_greekmu)];
    h.txt_scalevec(2).String = [num2str(str2double(h.txt_scalevec(2).String)*1000), ...
        sprintf(' %sV',char_greekmu)];
    
else
    h.txt_scalevec(1).String = [h.txt_scalevec(1).String, ' mV'];
    h.txt_scalevec(2).String = [h.txt_scalevec(2).String, ' mV'];
    
    
end

h.txt_scalevec(1).String = [h.txt_scalevec(1).String, ' mV'];
h.txt_scalevec(2).String = [h.txt_scalevec(2).String, ' mV'];

if ~isempty(binerr)
    h2_ = local_plotCircData(theta,RHO,ColorSpec);
    set(h2_.Children,'EdgeColor','none','FaceAlpha',0.2)
    
end

h2 = local_plotCircData(theta,rho,ColorSpec);



if ~showpatch
    h2.Visible = 'off';
end

%% Add mark and text

hold on

L = meanvec.length;

xmark = maxradius/L*real(meanvec.vec);
ymark = maxradius/L*imag(meanvec.vec);
handles.circ.hvec = plot(xmark,ymark,'+','MarkerSize',30,'Color','k');

text(xmark*1.1,ymark*1.1,sprintf('%.1f%s',...
    meanvec.degree,char(176)));

%% add group vectors

h3.vec = plot([0,real(meanvec.vec)],[0,imag(meanvec.vec)],'Color',ColorSpec,...
    'LineWidth',5,'DisplayName','Group Mean');
hold off

if ~showmeanvector
    h3.vec.Visible = 'off';
end
%% plot mean vectors for each bootstrap trial
if showbtstrpvecs
    axes(h2.Parent)
    hold on
    
    btvecmean = rad2cmp([btstrp(:).radian]').*[btstrp(:).length]';
    
    h3.bt = plot(h2.Parent,real(btvecmean),imag(btvecmean),'Color','k',...
        'Marker','.','LineStyle','none','Tag','Bootstrap mean',...
        'DisplayName','Bootstrap Mean'); % usually invisible
    hold off
end

%% add text labels

if ~isempty(titlestr)
    titleh = title(titlestr,'Interpreter','none','Visible','on','Units','normalized');
    titleh.Position = [0.5 1.15 0];
end

txth = local_placetext(btstrp,meanvec,h.axh);

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

function [axrad,binmean,binerr,axh,zeropos, plotdir, showmeanvector,showpatch,...
    showbtstrpvecs,ColorSpec, radius,meanvec,btstrp,titlestr,yunit] = local_parser(varargin)

axh = [];

if ishandle(varargin{1})
    % [h,N] = K_plotLinearPhaseWave(axh, _____)
    
    axh = varargin{1};
    assert(isscalar(axh) && strcmpi(axh.Type,'axes'))
    
    varargin = varargin(2:end);
   
end

% initialize
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

p.addParameter('Title','',@(x) ~isempty(x) &&...
    ischar(x) &&...
    isrow(x));

p.addParameter('ErrorRange','none',@(x) ismember(x,{'std','sem','none'}));

p.addParameter('ShowMeanVector',true,@(x) x ==0 || x==1);

p.addParameter('ShowPatch',true,@(x) x ==0 || x==1);

p.addParameter('ShowBootstrapVectors',false,@(x) x ==0 || x==1);

p.addParameter('Radius','auto',@(x) isequal(x,'auto') || isscalar(x) &&...
    x > 0);

p.addParameter('ZeroPos','top',@(x) ~isempty(x) && ischar(x) ...
    && isrow(x) && ismember(lower(x),{'top','left','right','bottom'}));

p.addParameter('Direction','clockwise',@(x) ~isempty(x) && ischar(x) && isrow(x) ...
    && ismember(lower(x),{'clockwise','anti'}));

p.addParameter('YUnit','mV2mV',@(x) ismember(lower(x),{'mv2mv','mv2microv'}));

p.parse(varargin{:});
    

S          = p.Results.S;
ColorSpec  = p.Results.Color;
errorrange = p.Results.ErrorRange;
axrad      = S.axrad;
binmean    = S.binmean;
    
showmeanvector = p.Results.ShowMeanVector;
showbtstrpvecs = p.Results.ShowBootstrapVectors;

showpatch      = p.Results.ShowPatch;
radius         = p.Results.Radius;


zeropos  = p.Results.ZeroPos;
plotdir  = p.Results.Direction;
titlestr = p.Results.Title;
btstrp   = S.bootstrap;
meanvec  = S.meanvec;
yunit    = lower(p.Results.YUnit);


switch errorrange
    case 'std'
        binerr = S.binstd;
    case 'sem'
        binerr = S.binsem;
    case 'none'
        binerr = [];
end


end

%--------------------------------------------------------------------------

function txth = local_placetext(btstrp,meanvec,axh)


if isempty(btstrp) || isempty(btstrp.p_lessthan)
    btstrpstr = '';
    
else
    btstrpstr = sprintf('Bootstrap test: p < %f\n', ...
        btstrp.p_lessthan);
end

thestr = sprintf(['%s',...
    'Phase of mean vector: %.3f%s\n',...
    'Length of mean vector: %.3e mV\n'], ...
    btstrpstr,...
    meanvec.degree,char(176),...
    meanvec.length);

axes(axh);

txth = text(1.2, -0.15, thestr,...
    'HorizontalAlignment', 'right', 'VerticalAlignment','bottom',...
    'Units', 'normalized','Tag','Annotation Text');

end

%--------------------------------------------------------------------------

function h = local_plotCircData(theta,rho,color)
%
% h = local_plotCircData(theta,rho,color)
%
% See also
% scr2016_07_04_223502_K_PhaseWaveWave_draft

% linh = gobjects(2,length(rho));

%% Negative rho

rhoneg = rho;
rhoneg(rho > 0) = 0; % to make it go through (0,0)

binmeanZ__ = rad2cmp(theta).* -rhoneg; % negative values are flipped over
x = real(binmeanZ__);
y = imag(binmeanZ__);

% X = [zeros(size(x)),x]';
% Y = [zeros(size(y)),y]';

ptc(1,1) = patch(x,y,'k',...
    'EdgeColor','k','LineWidth',1,'FaceAlpha',0.3,'EdgeColor','k',...
    'Tag','negative rho');

% linh(1,:) = line(X,Y,'LineStyle','-','Color','k','LineWidth',2,...
%     'Tag','negative rho','DisplayName','Negative rho')';

%% Positive rho
rhopos = rho;
rhopos(rho < 0) = 0;  % to make it go through (0,0)

binmeanZ = rad2cmp(theta).* rhopos;
x = real(binmeanZ);
y = imag(binmeanZ);

% X = [zeros(size(x)),x]';
% Y = [zeros(size(y)),y]';

ptc(2,1) = patch(x,y,color,...
    'EdgeColor',color,'LineWidth',1,'FaceAlpha',0.3,'EdgeColor',color,...
    'Tag','positive rho');

% linh(2,:) = line(X,Y,'LineStyle','-','Color',color,'LineWidth',2,...
%     'Tag','positive rho','DisplayName','Positive rho');

h = hggroup('Tag','Circular Data','DisplayName','Circular Data');

set(ptc,'Parent',h);

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

%--------------------------------------------------------------------------

function y = rad2cmp(x)

y = exp(1i .* x);
end            
    
 



 


