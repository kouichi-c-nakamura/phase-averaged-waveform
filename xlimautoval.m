function out = xlimautoval(x)
%
% reverse engineering of computation of XLim values by xlim('auto')
%
% out = xlimautoval(x)
%
% INPUT ARGUMENTS
% x           [low high]
%             XLim values
%
%
% OUTPUT ARGUMENTS
% out         [low high]
%             Modified XLim values floored and ceiled to the nearest
%             'simple' numbers.
%
% See also
% xlim,ylim,zlim, xlimautoval_test
% doc Axes Properties
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% kouichi.c.nakamura@gmail.com
% 09-Feb-2017 12:29:16



p = inputParser;
p.addRequired('x',@(x) isnumeric(x) && numel(x) == 2 && isrow(x));
p.parse(x);

y = [0 1];
f = figure('Visible','off');
plot(x,y);
a = gca;
out = a.XLim;
close(f)

end

%--------------------------------------------------------------------------

function dummy
% legacy code

p = inputParser;
p.addRequired('x',@(x) isnumeric(x) && numel(x) == 2 && isrow(x));
p.parse(x);



logint = floor(log10(abs(x)));
logintdiff = floor(log10(abs(diff(x))));


signedexp = zeros(1,2);
signedexp(1) = local_oneIntegerPart(max(logint),x(1)); % signedexp(1) to the power of 10
signedexp(2) = local_oneIntegerPart(max(logint),x(2));

signedexpdiff = local_oneIntegerPart(logintdiff,diff(x));


out = zeros(1,2);
X = zeros(1,2);

% the difference is now a number whose integer part with 1 digit; eg. 7.356
X(1) = x(1)/signedexpdiff; 
X(2) = x(2)/signedexpdiff;

if X(2) - X(1) > 5
    out(1) = floor(X(1))*signedexpdiff;
    out(2) = ceil(X(2))*signedexpdiff;

elseif X(2) - X(1) <= 5 && X(2) - X(1) > 2
    % 0.5 floor/ceil
    array = 0:0.5:10;
    
    ind1 = find( array <= X(1),1,'last');
    out(1) = array(ind1)*signedexpdiff;
    
    ind2 = find( array >= X(2),1,'first');
    out(2) = array(ind2)*signedexpdiff;
    
elseif X(2) - X(1) <= 2 && X(2) - X(1) > 1
      % 0.2 floor/ceil
    array = 0:0.2:10;
    
    ind1 = find( array <= X(1),1,'last');
    out(1) = array(ind1)*signedexpdiff;
    
    ind2 = find( array >= X(2),1,'first');
    out(2) = array(ind2)*signedexpdiff;    
    
elseif X(2) - X(1) <= 1 
    
    array = 0:0.1:10;
    
    ind1 = find( array <= X(1),1,'last');
    out(1) = array(ind1)*signedexpdiff;
    
    ind2 = find( array >= X(2),1,'first');
    out(2) = array(ind2)*signedexpdiff; 
    
end
        


end

%--------------------------------------------------------------------------

function signedexp = local_oneIntegerPart(logint,x)

if x > 0
    signedexp = sign(x)*10^logint;
elseif x < 0
    signedexp = sign(x)*10^logint;
else
    signedexp = 10^0;
end

end