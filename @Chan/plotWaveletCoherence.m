function [h,out] = plotWaveletCoherence(obj,obj2,varargin)
% plotWaveletScalogram draws wavelet scalogram with wcoherence
%
% [h,out] = plotWaveletCoherence(obj,obj2)
% [h,out] = plotWaveletCoherence(____,)
% [h,out] = plotWaveletCoherence(axh,____)
%
% INPUT ARGUMENTS
% obj, obj2   Chan objects that share sampling frequency SRate and data
%             length Length
%
% OPTIONAL PARAMETER/VALUE PAIRS
% 'PlotType'  'log' (default) | 'linear' | 'none'
%             'log' will use cwt's plotting capability, which plots image
%             with log Y axis. 'linear' will plot surface with linearly
%             spaced Y axis. 'none' will only export data without producing
%             a plot.
%
% 'ShowCOI'   true (default) | false | 1 | 0
%             Whether to show Cone of Influence, which indicates where edge
%             effects occur in the CWT
%
% 'ShowArrows' 
%             true (default) | false | 1 | 0
%             Whether to show arrows for phase lag of obj2 with respect to obj
%
%
% OUTPUT ARGUMENTS
% h           structure for graphic objects
%
% out         structure for cwt outputs
%             out.wt
%             out.f
%             out.coi
%
% See also
% wcoherence, Chan.plotWaveletScalogram
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% kouichi.c.nakamura@gmail.com
% 21-Apr-2017 16:08:06




if ~isempty(varargin) && isscalar(varargin{1}) && isgraphics(varargin{1},'axes')
    axh = varargin{1};
    axes(axh);
    varargin = varargin(2:end);
else
    axh = [];
end


p = inputParser;
p.addRequired('obj');
p.addRequired('obj2',@(x) isa(x,'Chan'));

p.addParameter('PlotType','log',@(x) ischar(x) && ismember(x,{'log','linear','none'}));
p.addParameter('ShowCOI',true,@(x)isscalar(x) && x == 0 || x == 1);
p.addParameter('ShowArrows',true,@(x)isscalar(x) && x == 0 || x == 1);

p.parse(obj,obj2,varargin{:});

PlotType = p.Results.PlotType;
showCOI = p.Results.ShowCOI;
showArrrows = p.Results.ShowArrows;

S = p.Unmatched;

if ~isempty(S)
    fn = fieldnames(S)';
    vals = struct2cell(S)';
    
    C = cell(1,length(fn)*2);
    for i = 1:length(fn)
        
        C((i-1)*2+1) = fn(i);
        C((i-1)*2+2) = vals(i);
        
    end
end


if string(PlotType) == 'none'
    axh = [];
else
    axh = gca;
end

assert(obj.SRate == obj2.SRate,...
    'Sampling rate (SRate property) of obj2 is different from that of obj')


[wcoh,wcs,f,coi] = wcoherence(obj.Data,obj2.Data,obj.SRate,C{:});

out.wcoh = wcoh;
out.wcs  = wcs;
out.f = f;
out.coi = coi;

switch PlotType
    case 'log'
        axes(axh)
        wcoherence(obj.Data,obj2.Data,obj.SRate,C{:});
        
        h.axh = axh;
        h.imagesc = findobj(h.axh,'Type','image');
        h.surface = [];
        h.line = findobj(h.axh,'Type','line');
       
        if ~showCOI
            h.line.Visible = 'off';
        end
        
        if ~showArrrows
            h.line.Visible = 'off';
        end
        
        set(gca,'TickDir','out')
        
    case 'linear'

        sf = surf(axh,obj.time,f,wcoh);
        axis tight
        view(0,90);
        sf.LineStyle = 'none';
        cbh = colorbar;
        ylimit = ylim(axh);
        
        lh = line(axh,obj.time,coi,  ...
            'LineStyle','--','Color','w','LineWidth',2);
        ylim(axh,ylimit);
        
        if ~showCOI
            lh.Visible = 'off';
        end
        
        xlabel('Time (sec)');
        ylabel('Frequencuy (Hz)')
        axh.TickDir = 'out';
        
        cbh.Label.String = 'Magnituide';
        cbh.TickDirection = 'out';
        
        h.axh = axh;
        h.imagesc = [];
        h.surface = sf;
        h.line = lh;
       

    otherwise
        
        h.axh = axh;
        h.imagesc = [];
        h.surface = [];
        h.line = [];
        
end

end