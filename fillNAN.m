function output = fillNAN(xinput, yinput, xreference)
% fillNAN - Fill NaN values in yinput by interpolating based on xinput and
% xreference.
% Syntax:
%   youtput = fillNAN(xinput, yinput, xreference)
% Inputs:
%   xinput      - Vector of x values corresponding to yinput.
%   yinput      - Vector of y values with possible NaNs.
%   xreference  - Vector of x values where youtput is desired.
% Outputs:
%   youtput     - Vector of y values at xreference with NaNs filled.

% Ensure column vectors
xinput = xinput(:);
yinput = yinput(:);
xreference = xreference(:);
output = nan(length(xreference), 1);

% find the different segments
if length(xinput) > length(xreference)
    error('xinput should not be longer than xreference');
else
    j = 1;
    for i = 1:length(xreference)
        if xreference(i) == xinput(j)
            output(i) = yinput(j);
            if j < length(xinput)
                j = j + 1;
            end
        end
    end

end