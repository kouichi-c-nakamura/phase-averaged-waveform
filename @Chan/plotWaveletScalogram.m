function [h,out] = plotWaveletScalogram(obj,varargin)
% plotWaveletScalogram draws wavelet scalogram with cwt
%
% [h,out] = plotWaveletScalogram(obj)
% [h,out] = plotWaveletScalogram(obj,wname)
% [h,out] = plotWaveletScalogram(____,)
% [h,out] = plotWaveletScalogram(axh,____)
%
% INPUT ARGUMENTS
% obj         a Chan object
%
%
% wname       'morse' (defualt) | 'amor' | 'bump'
%             (Optional) Wavelet name for cwt
%
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
% cwt, Chan.plotWaveletCoherence
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
p.addOptional('wname','morse',@(x) ismember(x,{'morse', 'amor', 'bump'}));
p.addParameter('PlotType','log',@(x) ischar(x) && ismember(x,{'log','linear','none'}));
p.addParameter('ShowCOI',true,@(x)isscalar(x) && x == 0 || x == 1);

p.parse(obj,varargin{:});

wname = p.Results.wname;
PlotType = p.Results.PlotType;
showCOI = p.Results.ShowCOI;


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



[wt,f,coi] = cwt(obj.Data,wname,obj.SRate,C{:});

out.wt = wt;
out.f = f;
out.coi = coi;

switch PlotType
    case 'log'
        axes(axh)
        cwt(obj.Data,wname,obj.SRate,C{:});
        
        h.axh = axh;
        h.imagesc = findobj(h.axh,'Type','image');
        h.surface = [];
        h.line = findobj(h.axh,'Type','line');
       
        if ~showCOI
            h.line.Visible = 'off';
        end
        
    case 'linear'
        
        sf = surf(axh,obj.time,f,abs(wt));
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