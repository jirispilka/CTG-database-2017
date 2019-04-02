function Y = fastmovav(X, w)

% Y = fastmovav(X, w)
% Returns the moving average for the columns of matrix X, given a running
% window w. The result is very fast and, possibly, slightly inaccurate. For
% more accurate (but slower) results, use accuratemovav instead.
%
% The first element of each column of Y is the sample mean of the first w
% elements of the corresponding column of X. Thus, if X is m-by-n, Y is
% (m-w+1)-by-n. The sliding window w must be a scalar greater than 1 and
% not greater than m.
%
% In order to use this function you need to compile its mex file first:
%     mex fastmovingaverage.c

%*************************************************************************%

% % Example
%
% % Some data
% T = 5000; N = 300; X = cumsum(randn(T, N)); w = ceil(T / 2);

% % Fast moving average
% tic, Y = fastmovav(X, w); toc
% % Elapsed time is 0.034531 seconds.
%
% % Check slow alternative with loops
% tic
% Z = zeros(size(X, 1) - w + 1, size(X, 2));
% for hh = 1:(size(X, 1) - w + 1)
%   Z(hh, :) = mean(X(hh:hh+w-1, :));
% end
% toc
% % Elapsed time is 43.973435 seconds.

% % Check accuracy
% all(abs(Y(:) - Z(:)) < 1e-12)

%*************************************************************************%
%                                                                         %
%            Author: Francesco Pozzi                                      %
%            E-Mail: francesco.pozzi@anu.edu.au                           %
%            Date: 11 December 2013                                       %
%                                                                         %
%*************************************************************************%

X = round(X'*10000);

chck = isnumeric(X) & isreal(X) & (ndims(X) == 2) & all(size(X) > 0);
chck = chck & isnumeric(w) & isreal(w) & isscalar(w);
if chck, chck = (w > 1) & (w <= size(X, 1)); end
if ~chck
  error('Check input: X must be a numeric, real, 2D array - w must be a positive integer.')
end

w = ceil(w);
Y = fastmovingaverage(X, w);

w1 = ceil(w/2)-1; 

% take care of boundary effect
ybegin = slidingAvg(X(1:w),w);

if rem(w,2) == 0
    yend = slidingAvg(X(end-w:end),w);
    yend = yend(end-w1:end);
else
    yend = slidingAvg(X(end-w+1:end),w);
    yend = yend(end-w1+1:end);
end

Y = [ybegin(1:w1); Y; yend];
Y = Y'/10000;

% while length(Y) > length(X)
%     Y = Y(1:end-1);
% end

