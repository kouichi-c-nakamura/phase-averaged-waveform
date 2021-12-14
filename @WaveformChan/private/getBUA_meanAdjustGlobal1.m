function Wadj = getBUA_meanAdjustGlobal1(W,indRep,befp,aftp,xdata,ydata)
%
% See also
% getBUA, getBUA/local_replacewithnan
% getBUA_meanAdjustGlobal2
 
srdWidth = round(0.003/W.SInterval);

ind_t0 = find(xdata >= 0,1,'first');

N = length(ydata);
acom = @(x) accomodate(x,N);


meanRep = nanmean(ydata(acom(ind_t0 - befp):acom(ind_t0 + aftp)));

meanBef = nanmean(ydata(acom(ind_t0 - befp - srdWidth):acom(ind_t0 - befp -1)));

meanAft = nanmean(ydata(acom(ind_t0 + aftp +1):acom(ind_t0 + aftp + srdWidth)));

meanBefAft = mean([meanBef,meanAft]);

meanDiff = meanBefAft - meanRep;

Wadj = W;
Wadj.Data(indRep) = W.Data(indRep) + meanDiff;

Wadj.ChanTitle = [W.ChanTitle, ' meanadjusted']; % or Header?

end

%--------------------------------------------------------------------------

function x = accomodate(x,N)
if x < 1
    x = 1;
    warning
elseif x > N
    x = N;
    warning
end

end