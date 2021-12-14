function iseq = eq(obj1, obj2)
% iseq = eq(obj1, obj2)
%TODO
%   C:\Program Files\MATLAB\R2012b\toolbox\matlab\timeseries\@timeseries\eq.m

% If one object is empty, return []
if isempty(obj1) && isempty(obj2)
    iseq = [];
    return
end

% If object is scalar do an element-wise comparison against the scalar
% timeseries
if numel(obj1)==1 && ~isempty(obj2)
    iseq = false(size(obj2));
    for k=1:numel(obj2)
        iseq(k) = localeq(obj1,obj2(k));
    end
    return
end

% Check that the sizes agree
if ~isequal(size(obj1), size(obj2))
    error(message('K:Chan:eq:sizemismatch'))
end

% Do an element-wise comparison
iseq = false(size(obj2));
for k=1:numel(obj1)
    iseq(k) = localeq(obj1(k), obj2(k));
end

end


function eqstate = localeq(obj1, obj2)

if ~strcmp(class(obj1), class(obj2)) % different class
    eqstate = false;
else
    
    propnames = properties(obj1);
    eqstate = true;
    
    for i = 1:length(propnames)
        if  ~isequaln(obj1.(propnames{i}), obj2.(propnames{i}))
            eqstate = false;
            
            % disp(propnames{i});
            break
        end
    end
end
end
