function val_ = get(a, property)
% function val_ = get(a, property)
%
% Get method for class analog
%
% Input:
%   a        ... the analog object
%   property ... string name of property to return
%
% Output:
%   val_ ... value of the named property

% Copyright 2005 by Joshua I. Gold
%   University of Pennsylvania

val_ = a.(property);
