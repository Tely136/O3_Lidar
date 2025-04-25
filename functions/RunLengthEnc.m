% https://www.mathworks.com/matlabcentral/answers/1663659-identify-consecutive-occurrence-of-1-in-a-binary-array
% created by Jan (2022)
function [b, n] = RunLengthEnc(x)
% See also: https://www.mathworks.com/matlabcentral/fileexchange/41813-runlength
d = [true; diff(x(:)) ~= 0];   % TRUE if values change
b = x(d);                      % Elements without repetitions
n = diff(find([d', true]));    % Number of repetitions
end