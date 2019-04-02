function out = slidingAvg(in, N, bPlot)
%% function STAN_slidingavg
%The function 'slidingavg' implements a one-dimensional filtering,
%applying a sliding window to a sequence.
%Such filtering replaces the center value in the window with the average 
%value of all the points within the window. 
%When the sliding window is exceeding the lower or upper boundaries 
%of the input vector INPUT_ARRAY, the average is computed among the 
%available points. 
%
%Input parameters:
% - in ..... input signal
% - N ..... scale of the window
% - bPlot.....variable if 1 will enable ploting of the function results
%
%Output parameters:
% - out ..... Filtered input
%
%Created in � 2002 by Michele Giugliano & Maura Arsiero
%<http://www.giugliano.info>
%
% Changed by v.ch.
%   ChangeLog on 9.04.2008:
%                - variable names
%                - complexness of the algorithm
%                - ploting introduced
%

if ~exist('bPlot','var')
    bPlot = false;
end

% If the input array is empty or N is non-positive
if (isempty(in)) || (N<=0)                                              
    disp(sprintf('SlidingAvg: (Error) empty input data or N null.'));     
    return;                                                               
end 

% If the number of neighbouring points over which the sliding
% average will be performed is '1', then no average actually occur and
% OUTPUT_ARRAY will be the copy of INPUT_ARRAY and the execution of the
% routine is stopped.
if (N==1)                                                              
    out = in;                                                             
    return;                                                               
end                                                               

% The length of the input data structure is acquired to later evaluate the
% 'mean' over the appropriate boundaries.
nx   = length(in);            

% If the number of neighbouring points over which the sliding
% average will be performed is large enough, then the average actually
% covers all the points of INPUT_ARRAY, for each index of OUTPUT_ARRAY
% and some CPU time can be gained by such an approach. 
if (N>=(2*(nx-1)))                                                     
    out = mean(in)*ones(size(in));                                        
    return;                                                               
end

% In all the other situations, the initialization of the output data
% structure is performed.
out = zeros(size(in));         
if rem(N,2)~=1                 
    m = N/2;                   
else                           
    m = (N-1)/2;               
end 

% For each element (i-th) contained in the input numerical array, a check
% must be performed: If not enough points are available on the left of the
% i-th element then we proceed to evaluate the mean from the first element
% to the (i + m)-th. If enough points are available on the left and on the
% right of the i-th element then we proceed to evaluate the mean on 2*m
% elements centered on the i-th position. If not enough points are
% available on the rigth of the i-th element. Then we proceed to evaluate
% the mean from the element (i - m)-th to the last one. If not enough
% points are available on the left and on the rigth of the  % i-th element
% then we proceed to evaluate the mean from the first element to the last.   
for i=1:nx,                                                 
    if ((i-m) < 1) && ((i+m) <= nx)                             
        out(i) = mean(in(1:i+m));                                 
    elseif ((i-m) >= 1) && ((i+m) <= nx)                        
        out(i) = mean(in(i-m:i+m));                              
    elseif ((i-m) >= 1) && ((i+m) > nx)                          
        out(i) = mean(in(i-m:nx));                               
    elseif ((i-m) < 1) && ((i+m) > nx)                          
        out(i) = mean(in(1:nx));                                  
    end 
end

if bPlot,
    figure;
    plot(in); hold on;
    plot(out,'r');
end
