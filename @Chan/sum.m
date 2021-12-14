function [ objout ] = sum( obj1, varargin )
%[ objout ] = sum( this, varargin )
%TODO

% C:\Program Files\MATLAB\R2012b\toolbox\matlab\timeseries\@timeseries\sum.m

if ~isa(obj2, class(obj1))
    error('K:Chan:sum:classmismatch',...
        'Class of %s doesn''t match.', inputname(2));
end

if obj1.Length ~= obj2.Length
    error('K:Chan:sum:lengthmismatch',...
        'Length doesn''t match');
end

obj3 = obj;
obj3.Data = obj1.Data + obj2.Data;

%TODO how to accept multiple inputs?

timeseries


end

