function TF = iscolorspec(colorSpec)
% iscolorspec returns true when colorSpec is a valid colorspec or false if
% not.
%
% TF = iscolorspec(colorSpec)
%
% INPUT ARGUMENTS
% colorSpec   RGB triplet | short name | long name | [R,G,B,alpha] quadruplets
%             iscolorspec also accepts [R,G,B,alpha] quadruplets because it
%             is undocumented but accepted color properties for line
%             objects etc. For further details, see
%             http://undocumentedmatlab.com/blog/plot-line-transparency-and-color-gradient
%
% OUTPUT ARGUMENT
% TF          logical scalar
%
%
% See also
% doc colorspec
% http://uk.mathworks.com/help/matlab/ref/colorspec.html
%
% Written by
% Kouichi C. Nakamura, Ph.D
% Kyoto University
% Nov 2015
% kouichi.c.nakamura@gmail.com

TF = false;

if ~isempty(colorSpec) && isrow(colorSpec)  
    
    if  isnumeric(colorSpec) ...
           && (length(colorSpec) == 3 || length(colorSpec) == 4)...
           &&  all(colorSpec >= 0 & colorSpec <= 1)
       
       TF = true;

    elseif ischar(colorSpec)
        if  ismember(lower(colorSpec), {'y','m','c','r','g','b','w','k'})
            TF = true;
            
        elseif ismember(lower(colorSpec), {'yellow','magenta','cyan','red','green','blue','white','black', 'none'})
            TF = true;
        end
    end

end

end