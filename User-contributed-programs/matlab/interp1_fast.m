function v = interp1(x,y,u)
%
% same as interp1, but with error checking stripped out for speed
% assumes input of form interp1_fast(x,y,xi).  Assumes that x goes from
% most to least.
%
%INTERP1 1-D interpolation (table lookup).
%   YI = INTERP1(X,Y,XI) interpolates to find YI, the values of
%   the underlying function Y at the points in the vector XI.
%   The vector X specifies the points at which the data Y is
%   given. If Y is a matrix, then the interpolation is performed
%   for each column of Y and YI will be length(XI)-by-size(Y,2).
%
%   YI = INTERP1(Y,XI) assumes X = 1:N, where N is the length(Y)
%   for vector Y or SIZE(Y,1) for matrix Y.
%
%   Interpolation is the same operation as "table lookup".  Described in
%   "table lookup" terms, the "table" is [X,Y] and INTERP1 "looks-up"
%   the elements of XI in X, and, based upon their location, returns
%   values YI interpolated within the elements of Y.
%
%   YI = INTERP1(X,Y,XI,'method') specifies alternate methods.
%   The default is linear interpolation.  Available methods are:
%
%     'nearest'  - nearest neighbor interpolation
%     'linear'   - linear interpolation
%     'spline'   - piecewise cubic spline interpolation (SPLINE)
%     'pchip'    - piecewise cubic Hermite interpolation (PCHIP)
%     'cubic'    - same as 'pchip'
%     'v5cubic'  - the cubic interpolation from MATLAB 5, which does not
%                  extrapolate and uses 'spline' if X is not equally spaced.
%
%   YI = INTERP1(X,Y,XI,'method','extrap') uses the specified method for
%   extrapolation for any elements of XI outside the interval spanned by X.
%   Alternatively, YI = INTERP1(X,Y,XI,'method',EXTRAPVAL) replaces
%   these values with EXTRAPVAL.  NaN and 0 are often used for EXTRAPVAL.
%   The default extrapolation behavior with four input arguments is 'extrap'
%   for 'spline' and 'pchip' and EXTRAPVAL = NaN for the other methods.
%
%   For example, generate a coarse sine curve and interpolate over a
%   finer abscissa:
%       x = 0:10; y = sin(x); xi = 0:.25:10;
%       yi = interp1(x,y,xi); plot(x,y,'o',xi,yi)
%
%   See also INTERP1Q, INTERPFT, SPLINE, INTERP2, INTERP3, INTERPN.

%   Copyright 1984-2001 The MathWorks, Inc.
%   $Revision: 5.38 $  $Date: 2001/04/15 11:59:10 $

% Determine input arguments.

x=flipud(x(:));
y=flipud(y(:));

[m,n] = size(y);

h=diff(x);

siz=size(u);
u=u(:);
p =[];

v = zeros(size(u,1),n*size(u,2));
q = length(u);
p = 1:q;

[ignore,k] = histc(u,x);
k(u<x(1) | ~isfinite(u)) = 1;
k(u>=x(m)) = m-1;

s = u - x(k);
for j = 1:n
  del = diff(y(:,j))./h;
  v(p,j) = y(k,j) + s.*del(k);
end

v=reshape(v,siz);

%if isempty(p)
%  p = 1:length(u);
%end
k = find(u<min(x) | u>max(x));
v(p(k),:) = nan;
