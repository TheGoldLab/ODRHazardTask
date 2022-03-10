function sem = rocSE(x,y,n)
% Standard error for area under ROC
% usage: sem = rocSE(x,y,n)
% Third argument is optional.  It controls number of points in the
% ROC. 100 is default.
% x and y can be vectors or matrices.  If matrices, then sem is a row
% vector corresponding to ROCs formed from each pair of corresponding
% columns in x and y.  Note that x and y can be padded with NaNs (see
% also nanroc.m)
%
% Based on Hanley & McNeil, 1982

% 2/22/00 mns wrote it
% 2/23/00 added capacity for matrices padded with NaNs

% test 
% $$$ x = normrnd(0,1,100,1);
% $$$ y = normrnd(1,1,100,1);
% $$$ n = 100;
% $$$ error('comment out test code')

% routine starts here

if nargin < 3
  n = 100;				% default for number of perms
end

if isvector(x)
  x = x(:);
  xvectorflag = 1;
else
  xvectorflag = 0;
  X = x;
end


if isvector(y)
  y = y(:);
  yvectorflag = 1;
else
  yvectorflag = 0;
  Y = y;
end


if xvectorflag & yvectorflag
  % get ROC area
  W = rocN(y,x,n)
  
  nx = length(x);
  ny = length(y);
  
  % Estimate Q1, the probability that two values from y will exceed x
  Y1 = repmat(y,[1 ny]);
  Y2 = repmat(y,[1 ny]);
  for i = 2:ny
    ind = [i:ny 1:i-1]';
    Y2(:,i) = y(ind);
  end
  Y1 = Y1(:);
  Y2 = Y2(:);
  count = 0;
  for i = 1:length(Y1)
    count = count + sum(Y1(i)>x & Y2(i)>x); 
  end
  Q1 = count / (length(Y1) * nx);
  
  % Q2, prob that one value in y will exceed two vals from x
  X1 = repmat(x,[1 nx]);
  X2 = repmat(x,[1 nx]);
  for i = 2:nx
    ind = [i:nx 1:i-1]';
    X2(:,i) = x(ind);
  end
  X1 = X1(:);
  X2 = X2(:);
  count = 0;
  for i = 1:length(X1)
    count = count + sum(X1(i)<y & X2(i)<y); 
  end
  Q2 = count / (length(X1) * ny);
  
  sem = sqrt((W*(1-W) + (ny-1) * (Q1 - W*W) +...
	      (nx-1) * (Q2 - W*W))/(nx * ny));
else
  % dealing with matrices, nan padded
  % go through each column
  if size(X,2) ~= size(Y,2)
    error('Matrices x and y must have same number of columns')
  end
  sem = repmat(nan,[1 size(X,2)]);	% allocate return arg
  for k = 1:size(X,2)
    x = X(:,k);
    x = x(finite(x));
    y = Y(:,k);
    y = y(finite(y));
    
    W = rocN(y,x,n);
  
    nx = length(x);
    ny = length(y);
    
    % Estimate Q1, the probability that two values from y will exceed x
    Y1 = repmat(y,[1 ny]);
    Y2 = repmat(y,[1 ny]);
    for i = 2:ny
      ind = [i:ny 1:i-1]';
      Y2(:,i) = y(ind);
    end
    Y1 = Y1(:);
    Y2 = Y2(:);
    count = 0;
    for i = 1:length(Y1)
      count = count + sum(Y1(i)>x & Y2(i)>x); 
    end
    Q1 = count / (length(Y1) * nx);
    
    % Q2, prob that one value in y will exceed two vals from x
    X1 = repmat(x,[1 nx]);
    X2 = repmat(x,[1 nx]);
    for i = 2:nx
      ind = [i:nx 1:i-1]';
      X2(:,i) = x(ind);
    end
    X1 = X1(:);
    X2 = X2(:);
    count = 0;
    for i = 1:length(X1)
      count = count + sum(X1(i)<y & X2(i)<y); 
    end
    Q2 = count / (length(X1) * ny);
    
    sem(k) = sqrt((W*(1-W) + (ny-1) * (Q1 - W*W) +...
		(nx-1) * (Q2 - W*W))/(nx * ny));
  end
end

