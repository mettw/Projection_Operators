% Author: MATLAB Copilot
% Date: 2024-08-02
% Function to test if a number is odd.
% Input:
%   - num: The number to be tested (integer).
% Output:
%   - out: Logical true (1) if the number is odd, false (0) otherwise.
% Example usage:
%   result = isodd(5); % result will be true (1)

function out = isodd(num)
    % Check if the number is odd by using the mod function
    out = mod(num, 2) == 1;
end
