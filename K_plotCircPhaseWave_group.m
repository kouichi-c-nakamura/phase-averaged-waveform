function h = K_plotCircPhaseWave_group(varargin)
% K_plotCircPhaseWave_group create the following out of phase values in
% radians. It works on single sample only.
%
%    1. Multiple vectors representing circular mean and vector length of
%    one sample.
%
%    2. Small circles at the perimeter that indicate circular mean of one
%    sample.
%
%    3. Also create a group vector with 'Tag' property set to 'Group Vector
%    CircMean', which represent circular mean of phases of each samples.
%    However, this is by default hidden by 'Visible','off'. You can see it
%    by:
%
%    set(findobj(gca,'Tag','Group Vector CircMean'),'Visible','on')
%
%    4. Also create a group vector with 'Tag' property set to 'Group Vector
%    Synthetic', which represent mean of sample vectors (computed as mean
%    of complex numbers). However, this is by default hidden by
%    'Visible','off'. You can see it by:
%
%    set(findobj(gca,'Tag','Group Vector Synthetic'),'Visible','on')
%
% h = K_plotCircPhaseWave_group(S)
% h = K_plotCircPhaseWave_group(axh, ______)
% h = K_plotCircPhaseWave_group(_____, 'Param', Value, ...)
% 
% INPUT ARGUMENTS
% 
% radians          A cell vector containing a vector of phase values in             
%                  radian. Each cell represent one sample.
%
% axh              Axes handle (optional)
% 
% OPTIONS (paramter/value pairs)
% 'Color'          ColorSpec 
%
% 'ZeroPos'        'top' (default) | 'right' | 'left' | 'bottom'
%
% 'Direction'      'clockwise' (default) | 'anti'
%
% 'Radius'         'auto' (default) | scalar number
%                  Specify the radius of the circular axes for sample
%                  vectors in mV.
%
% 'ShowGroupCirc'  true (default) | false
%
% 'ShowGroupSyn'   true (default) | false
%
% 'ShowSampleVec'  true (default) | false
%                  Shows sample vectors
%
% 'YUnit'          'mv2mv' (default) | 'mv2microv'
%                  'mv2microv' will 1000-fold values for Y axis and use
%                  microvolt(µV) as unit.
%
% 'RadiusForAngles'
%                  [1.2, 1.2] (default) | [r0, r90]
%                  Determines the positions of texts indicating angles.
%                  r0 is for 0 and 180 degrees, whereas r90 is for 90 and 270.
%
% 'RadianForStandardScales'
%                  pi*3/8 (default) | scalar
%                  Phase of the positions of texts indicating scales in radian
%
% 'RadianForAmplitudeScales'
%                  pi*1/8 (default) | scalar
%                  Phase of the positions of texts indicating scales in radian
%
% 'RadiusForScales'
%                  [0.65,1.2] (default) | [inner, outer] 
%                  Determines the positions of texts indicating scales
%
% OUTPUT ARGUMENTS
%
% h                A structure of handles for graphic objects.
%
%
% EXAMPLES
% ah2 = CircularPlotSingle_KCN(formated.cmean_group, formated.vec_group,...
%           formated.cmean, colorsp);
%
% See Also 
% K_plotCircPhaseWave_one, K_plotCircPhaseHist_one, K_PhaseHist,
% K_PhaseHist_test pvt_K_plotCircPhaseHist_parseInputs


%% parse input arguments


[axh,zeropos,plotdir,ColorSpec,radius,showsamplevec,showgroupcircmean,...
    showgroupsyn,titlestr,btstrp,circsh,meanvec,yunit,...
    radiusForAngles,radiusForScales,radianForStandardScales,radianForAmplitudeScales]...
    = local_parser(varargin{:});

%% Job

N = numel(meanvec);

if N == 0
   disp('No data to plot.') 
   h = [];
   return 
end


if isempty(axh) 
    h.fig = figure('Color',[1 1 1]); hold on;

    h.axh = gca;
    
    h.axh.Position = [0.1300    0.16    0.7    0.7];

else
    axes(axh);

    h.fig = gcf;
    set(gcf, 'Color', [1 1 1]);
    h.axh = axh;
end  

axis tight
axis equal;


radians  = vertcat(meanvec(:).radian);
veclen = vertcat(meanvec(:).length);

h = K_plotCircPhaseHist_one(h.axh,radians,...
    'Color',ColorSpec,'ZeroPos',zeropos,'Direction',plotdir,...
    ... 'RadiusForHistScale',radiusForHistScale, ...
    'RadiusForAngles',radiusForAngles, 'RadiusForScales',radiusForScales,...
    'RadianForScales',radianForStandardScales);


delete(findobj(h.axh,'Tag','Rose Histogram',...
    '-or','Tag','Histogram Scale Full',...
    '-or','Tag','Histogram Scale Half',...
    '-or','Tag','Annotation Text'))

set(findobj(h.axh,'Tag','Vector'),'Tag','Vector per Group');



%% Add sample vectors
axes(h.axh);
hold on

if isequal(radius,'auto')
    maxradius = max(xlimautoval([0,max(veclen)]));
else
    maxradius = radius;
end


h.vector = gobjects(N,1);
for i = 1:N

    z1 = rad2cmp(radians(i));
    
    h.vector(i) = plot(h.axh, [0, real(z1)*veclen(i)/maxradius], ...
        [0, imag(z1).*veclen(i)/maxradius],...
        'Color', ColorSpec , 'LineWidth', 1,'Tag','Vector per Sample');

end


tz1 =rad2cmp(radianForAmplitudeScales);


 if strcmpi(yunit,'mv2mv')
     h.txt_scalevec(3) = text(real(tz1)*radiusForScales(1), ...
         imag(tz1)*radiusForScales(1),...
         [num2str(maxradius/2),' mV'],...
         'Color', 'k','Tag','Scale Half (Sample)');
     
     h.txt_scalevec(4) = text(real(tz1)*radiusForScales(2), ...
         imag(tz1)*radiusForScales(2), ...
         [num2str(maxradius),' mV'],...
         'Color', 'k','Tag','Scale Full (Sample)');

 elseif strcmpi(yunit,'mv2microv')
     h.txt_scalevec(3) = text(real(tz1)*radiusForScales(1), ...
         imag(tz1)*radiusForScales(1),...
         [num2str(maxradius/2*1000),sprintf(' %sV',char_greekmu)],...
         'Color', 'k','Tag','Scale Half (Sample)');
     
     h.txt_scalevec(4) = text(real(tz1)*radiusForScales(2), ...
         imag(tz1)*radiusForScales(2), ...
         [num2str(maxradius*1000),sprintf(' %sV',char_greekmu)],...
         'Color', 'k','Tag','Scale Full (Sample)');
     
 end


if ~showsamplevec
    set(h.vector,'Visible','off');
    set(h.txt_scalevec(3:4),'Visible','off');
end

%% Add a synthetic group vector

cmp = rad2cmp(radians).*veclen;

synthvec = mean(cmp);

hold on
h.synthvec = plot(h.axh,[0,real(synthvec)/maxradius],...
    [0,imag(synthvec)/maxradius],...
    'Color','r','LineWidth',5,'Tag','Group Vector Synthetic');

h.synthmk = plot(h.axh,[0,real(synthvec)/abs(synthvec)],...
    [0,imag(synthvec)/abs(synthvec)],...
    'Color','r','LineStyle','none','Marker','*','MarkerSize',10,...
    'Tag','Mark on Perimeter for Group Vector Synthetic');
hold off

if ~showgroupsyn
    h.synthvec.Visible = 'off';
end

% rename tag
set(findobj(h.axh,'Tag','Mark on Perimeter per Group'),...
    'Tag','Mark on Perimeter for Group Vector CircMean');

set(findobj(h.axh,'Tag','Vector per Group'),...
    'Tag','Group Vector CircMean');


%% Appearance of Group vector
if ~showgroupcircmean
    set(findobj(h.axh,'Tag','Group Vector CircMean'),'Visible','off');
    set(findobj(h.axh,'Tag','Mark on Perimeter for Group Vector CircMean'),'Visible','off'); 
end

if ~showgroupsyn
    set(h.synthmk,'Visible','off');

    set(h.synthvec,'Visible','off');
    
end

uistack(findobj(gca,'Tag','Vector per Group'),'top')
uistack(h.synthvec,'top')


txth = local_placetext(btstrp,circsh,meanvec,h.axh);

if ~isempty(titlestr)
   title(titlestr,'Interpreter','none');   
end

end

%% LOCAL FUNCTIONS

%--------------------------------------------------------------------------

function [axh, zeropos, plotdir, ColorSpec,radius,showsamplevec,showgroupcircmean,...
    showgroupsyn,titlestr,btstrp,circsh,meanvec,yunit,...
    radiusForAngles,radiusForScales,radianForScales,radianForAmplitudeScales]...
    = local_parser(varargin)


if ishandle(varargin{1})
    
    axh = varargin{1};
    varargin = varargin(2:end);
    
else
   
    axh = [];
    
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

p.addParameter('ShowSampleVec',true,@(x) x ==0 || x==1);

p.addParameter('ShowGroupCirc',true,@(x) x == 1 || x == 0);

p.addParameter('ShowGroupSyn',true,@(x) x == 1 || x == 0);

p.addParameter('ZeroPos','top',@(x) ~isempty(x) && ischar(x) ...
    && isrow(x) && ismember(lower(x),{'top','left','right','bottom'}));

p.addParameter('Direction','clockwise',@(x) ~isempty(x) && ischar(x) && isrow(x) ...
    && ismember(lower(x),{'clockwise','anti'}));

p.addParameter('Radius','auto',@(x) isequal(x,'auto') || isscalar(x) &&...
    x > 0);

p.addParameter('YUnit','mV2mV',@(x) ismember(lower(x),{'mv2mv','mv2microv'}));

p.addParameter('RadiusForAngles',  [1.2,1.2], @(x) isrow(x) && numel(x) == 2 ...
    && isnumeric(x) && all(x >=0))

p.addParameter('RadiusForScales', [0.65,1.2], @(x) isrow(x) && numel(x) == 2 ...
    && isnumeric(x) && all(x >=0) && x(1) < x(2))

p.addParameter('RadianForScales', pi*3/8, @(x) isscalar(x) && isreal(x))

p.addParameter('RadianForAmplitudeScales', pi*1/8, @(x) isscalar(x) && isreal(x));

p.parse(varargin{:});
    

S          = p.Results.S;
ColorSpec  = p.Results.Color;
errorrange = p.Results.ErrorRange;
    
showsamplevec   = p.Results.ShowSampleVec;
showgroupcircmean  = p.Results.ShowGroupCirc;
showgroupsyn = p.Results.ShowGroupSyn;


zeropos  = p.Results.ZeroPos;
plotdir  = p.Results.Direction;
titlestr = p.Results.Title;
radius   = p.Results.Radius;
btstrp   = [S(:).bootstrap];
circsh   = [S(:).circshift];
meanvec  = [S(:).meanvec];
yunit    = lower(p.Results.YUnit);

radiusForAngles = p.Results.RadiusForAngles;
radiusForScales = p.Results.RadiusForScales;
radianForScales = p.Results.RadianForScales;
radianForAmplitudeScales = p.Results.RadianForAmplitudeScales;


switch errorrange
    case 'std'
        binerr = [S(:).binstd];
    case 'sem'
        binerr = [S(:).binsem];
    case 'none'
        binerr = [];
end


end

%--------------------------------------------------------------------------

function txth = local_placetext(btstrp,circsh,meanvec,axh)
%
% See also
% K_plotCircPhaseHist_group/local_placetext
N = length(meanvec);

% circmean of phase values and mean length

samplevecrad = [meanvec(:).radian]';
groupradmean = circ_mean(samplevecrad);
groupradstd = circ_std(samplevecrad);

groupveclen = circ_r(samplevecrad);

sampleveclen = [meanvec(:).length];
sampleveclenmean = mean(sampleveclen);
sampleveclenstd = std(sampleveclen);

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
    'Group vector length: %.3f\n', ...
    'Sample Vector length mean %s std: %3f %s %3f mV\n'], ...
    N,...
    pstr,...
    char(177),rad2deg(groupradmean),char(177),rad2deg(groupradstd),char(176),...
    groupveclen,...
    char(177),sampleveclenmean,char(177),sampleveclenstd);

axes(axh);  

txth = text(1.3, -0.15, thestr,...
    'HorizontalAlignment', 'right', 'VerticalAlignment','bottom',...
    'Units', 'normalized','Tag','Annotation Text',...
    'FontSize',9);

end

%--------------------------------------------------------------------------

function cmp = rad2cmp(rad)

cmp = exp(1i .* rad);

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