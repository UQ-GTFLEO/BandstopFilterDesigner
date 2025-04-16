function [maxTransmissionYn] = single_frequency_heatmap(frequency, totalSpacing, varargin)
%SINGLE_FREQUENCY_HEATMAP Summary of this function goes here
%   Detailed explanation goes here

p = inputParser;
validBounds = @(x) isnumeric(x) && isscalar(x);
validPoints = @(x) isnumeric(x) && isscalar(x) && x == fix(x);
validSubstrate = @(x) isnumeric(x) && isscalar(x) && x >= 1;
validHeights = @(x) isnumeric(x) && isscalar(x) && x > 0;
addOptional(p, 'minYn', -10, validBounds);
addOptional(p, 'maxYn', 10, validBounds);
addOptional(p, 'numPoints', 1000, validPoints);
addOptional(p, 'substrateBackingEps', 3.55, validSubstrate);
addOptional(p, 'substrateBackingHeight', 0.508, validHeights);
parse(p, varargin{:});

if ~isscalar(totalSpacing) || ~(totalSpacing > 0)
    error("The value of 'totalSpacing' is invalid. It must satisfy the function: @(x)~isscalar(x)||~(x > 0).");
end

if ~isscalar(frequency) || ~(frequency > 0)
    error("The value of 'frequency' is invalid. It must satisfy the function: @(x)~isscalar(x)||~(x > 0).");
end

numPoints = p.Results.numPoints;
yns = linspace(p.Results.minYn, p.Results.maxYn, numPoints);


%% Define constants

% Scale fruencies
freq = frequency * 1e9;
S21s = zeros(numPoints, numPoints);
eps_r = p.Results.substrateBackingEps;


% Scale spacing
dSub = p.Results.substrateBackingHeight * 1e-3;
dAir = totalSpacing * 1e-3 - dSub;

% Base Case Air
lambdaAir = physconst('LightSpeed') ./ freq;
dWavelengthsAir = dAir ./ lambdaAir;
bdAir = 2 * pi * dWavelengthsAir;

% Substrate Case
lambdaSub = lambdaAir ./ sqrt(eps_r);
dSubWavelengths = dSub ./ lambdaSub;
bdSub = 2 * pi * dSubWavelengths;

% Free space and substrate impedances
Zd = 377/sqrt(eps_r);
Z0 = 377;

% Loop through the admittance values and populate the heatmap
for row = 1:numPoints
    for col = 1:numPoints
        
        Ys1 = (yns(col) * 1i) / Z0;
        Ys2 = (yns(row) * 1i ) / Z0;

        ABCD_Y1 = [1, 0; Ys1, 1];
        ABCD_Y2 = [1, 0; Ys2, 1];
        ABCD_TL_AIR = [cos(bdAir), 1i * Zd * sin(bdAir); 1i * sin(bdAir)/Zd, cos(bdAir)];
        ABCD_TL_SUB = [cos(bdSub), 1i * Zd * sin(bdSub); 1i * sin(bdSub)/Zd, cos(bdSub)];
        ABCD = ABCD_Y1 * ABCD_TL_SUB * ABCD_TL_AIR * ABCD_TL_SUB * ABCD_Y2;

        A = ABCD(1,1);
        B = ABCD(1,2);
        C = ABCD(2,1);
        D = ABCD(2,2);
        S21 = 2 / (A + B/Z0 + C * Z0 + D);
        S21s(row, col) = abs(S21);
    end
end

imagesc(yns, yns, S21s)
set(gca,'YDir','normal')
xlabel("Y1");
ylabel("Y2");
colorbar
colormap('jet')
axis square
set(gca,"FontSize",30)


% Get the exact maximums
sz = [numPoints numPoints];
[~, idx] = max(S21s, [], 'all');
[~, col] = ind2sub(sz, idx);
maxTransmissionYn = yns(col);


end