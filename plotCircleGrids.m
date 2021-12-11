function h = plotCircleGrids(varargin)
% h = plotCircleGrids()
% h = plotCircleGrids(axh)
% h = plotCircleGrids(axh,radius)
%
% INPUT ARGUMENTS
% 'axh'       (Optional) Axes object
%
% 'radius'    1 (default) | scalar
%             (Optional) Radius of the circle
%
% OPTIONAL PARAMETER/VALUE PAIRS
%
% 'RadiusForAngles'
%             [1.2, 1.2] (default) | [r0, r90]
%             Determines the positions of texts indicating angles.
%             r0 is for 0 and 180 degrees, whereas r90 is for 90 and 270.
%
% 'RadiusForScales'
%             [0.65,1.2] (default) | [inner, outer] 
%             Determines the positions of texts indicating scales
%
% 'RadianForScales'
%             pi*3/8 (default) | scalar
%             Phase of the positions of texts indicating scales in radian
%
%
% OUTPUT ARGUMENTS
% h           structure 
%
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% kouichi.c.nakamura@gmail.com
% 23-May-2017 14:23:38
%
% See also
% K_plotCircPhaseHist_one, K_plotCircPhaseHist_group, K_plotCircPhaseWave_one


p = inputParser;
p.addOptional('axh',[],@(x) isempty(x) || ishandle(x) && strcmpi(x.Type,'axes'))
p.addOptional('radius',1,@(x) isscalar(x) && isnumeric(x) && x >=0)

p.addParameter('RadiusForAngles',  [1.2,1.2], @(x) isrow(x) && numel(x) == 2 ...
    && isnumeric(x) && all(x >=0))
p.addParameter('RadiusForScales', [0.65,1.2], @(x) isrow(x) && numel(x) == 2 ...
    && isnumeric(x) && all(x >=0) && x(1) < x(2))
p.addParameter('RadianForScales', pi*3/8, @(x) isscalar(x) && isreal(x))

p.parse(varargin{:})

radius = p.Results.radius;
axh    = p.Results.axh;

radiusForAngles = p.Results.RadiusForAngles;
radiusForScales = p.Results.RadiusForScales;
radianForScales = p.Results.RadianForScales;


if isempty(axh) 
    
    h.fig = figure('Color',[1 1 1]); hold on;
    h.axh = gca;
    axh = gca;
    
    axh.Position = [0.1300    0.16    0.7    0.7];

else
    axes(axh);

    h.fig = gcf;
    set(gcf, 'Color', [1 1 1]);
    h.axh = axh;
    
end    

set(h.axh, ...
    'Visible','off',...
    'Box', 'off',...
    'PlotBoxAspectRatioMode', 'manual',...
    'Units', 'normalized');

rad2cmp = @(x) exp(1i * x);
zz = rad2cmp(linspace(0, 2*pi, 3600));

view(90,-90);

hold on
axis tight
axis equal

% draw concentric circles
h.outercirc = plot(axh, radius*real(zz),     radius*imag(zz),...
    'Color', 'k', 'LineWidth',2 ,'Tag','Outer Circle');
h.innercirc = plot(axh, radius*1/2*real(zz), radius*1/2*imag(zz),...
    'Color', [0.5 0.5 0.5],      'Tag','Inner Circle');

% draw vertical and holizontal lines
h.vertline = plot(axh, radius*[0 0],radius*[-1 1],'Color', [0.5 0.5 0.5],...
    'Tag','Vertical Line');
h.horzline = plot(axh, radius*[-1 1],radius*[0 0],'Color', [0.5 0.5 0.5],...
    'Tag','Horizontal Line');

hg = hggroup('Tag','Circular Grid');
set([h.outercirc,h.innercirc,h.vertline,h.horzline],'Parent',hg);

%%
h.t(1) = text(radius*radiusForAngles(1),  radius*0,    ...
    ['0',   char(176)],'Tag','Phase 0'); % char(176) == degree symbol

h.t(2) = text(radius*0,               radius*radiusForAngles(2),  ...
    ['90',  char(176)],'Tag','Phase 90');

h.t(3) = text(radius*-radiusForAngles(1), radius*0,    ...
    ['180', char(176)],'Tag','Phase 180');

h.t(4) = text(radius*0,               radius*-radiusForAngles(2), ...
    ['270', char(176)],'Tag','Phase 270');


hg = hggroup('Tag','Phase 0 90 180 270');
set(h.t,'Parent',hg);

set(h.t, 'HorizontalAlignment', 'center', 'FontSize', 12);

tz1 =rad2cmp(radianForScales);

h.txt_scalevec(1) = text(radius*real(tz1)*radiusForScales(1), ...
    radius*imag(tz1)*radiusForScales(1),...
    num2str(radius/2),...
    'Color', 'k','Tag','Scale Half');
h.txt_scalevec(2) = text(radius*real(tz1)*radiusForScales(2), ...
    radius*imag(tz1)*radiusForScales(2), ...
    num2str(radius),...
    'Color', 'k','Tag','Scale Full');

set(h.txt_scalevec, 'HorizontalAlignment', 'center');

axh.Visible = 'off';
h.fig.Color = 'w';

end