function [requiredImpedancesBelow, requiredImpedancesAbove, frequencies] = generate_required_impedance_curve(frequencies, totalSpacing, varargin)
%GENERATE_REQUIRED_IMPEDANCE_CURVE Vector of required impedance values.
%
% generate_required_impedance_curve(frequencies, totalSpacing) generates a
% vector of impedance values that optimally allow for transmission through
% a FSS for the given frequencies and spacing between layers that are
% specified.
%
% By default, plots are generated to represent transmission and of the
% returned impedance values. These can be disabled with
% generate_required_impedance_curve(frequencies, totalSpacing, 'enablePlot'
% 0).
%
% All FSSs are assumed to have a substrate backing.  By defult this is
% 0.508mm thick with a dielectric constant of 3.55.
%
% Example usage:
%   frequencies = [5:0.1:15];
%   totalSpacing = 2;
%   generate_required_impedance_curve(frequencies, ...
%       totalSpacing, ...
%       'enablePlot', 1, ...
%       'minYn', -10, ...
%       'maxYn', 10, ...
%       'numPoints', 1000, ...
%       'substrateBackingEps', 3.55, ...
%       'substrateBackingHeight', 0.508);

p = inputParser;
validPlotFlag = @(x) isnumeric(x) && isscalar(x) && (x == 0 || x == 1);
validBounds = @(x) isnumeric(x) && isscalar(x);
validPoints = @(x) isnumeric(x) && isscalar(x) && x == fix(x);
validSubstrate = @(x) isnumeric(x) && isscalar(x) && x >= 1;
validHeights = @(x) isnumeric(x) && isscalar(x) && x > 0;
addOptional(p, 'enablePlot', 1, validPlotFlag);
addOptional(p, 'minYn', -10, validBounds);
addOptional(p, 'maxYn', 10, validBounds);
addOptional(p, 'numPoints', 1000, validPoints);
addOptional(p, 'substrateBackingEps', 3.55, validSubstrate);
addOptional(p, 'substrateBackingHeight', 0.508, validHeights);
parse(p, varargin{:});

if ~isscalar(totalSpacing) || ~(totalSpacing > 0)
    error("The value of 'spacing' is invalid. It must satisfy the function: @(x)~isscalar(x)||~(x > 0).");
end

enablePlot = p.Results.enablePlot;
numPoints = p.Results.numPoints;
yns = linspace(p.Results.minYn, p.Results.maxYn, numPoints);

%% Define constants

% Scale fruencies
freqs = frequencies * 1e9;
bestYsBelow = zeros(size(freqs));
bestZsBelow = zeros(size(freqs));
bestYsAbove = zeros(size(freqs));
bestZsAbove = zeros(size(freqs));
S21s = zeros(numel(freqs), numPoints);
eps_r = p.Results.substrateBackingEps;


% Scale spacing
dSub = p.Results.substrateBackingHeight * 1e-3;
dAir = totalSpacing * 1e-3 - dSub;

% Base Case Air
lambdaAir = physconst('LightSpeed') ./ freqs;
dWavelengthsAir = dAir ./ lambdaAir;
bdsAir = 2 * pi * dWavelengthsAir;

% Substrate Case
lambdaSub = lambdaAir ./ sqrt(eps_r);
dSubWavelengths = dSub ./ lambdaSub;
bdsSub = 2 * pi * dSubWavelengths;

% Free space and substrate impedances
Zd = 377/sqrt(eps_r);
Z0 = 377;


%% Generate S21 parameters for each frequency and find ideal Y and Z.
for i = 1:numel(bdsAir)
    bdAir = bdsAir(i);
    bdSub = bdsSub(i);

    for point = 1:numPoints
        % FSS layer Admittance
        Ys = (yns(point) * 1i) / Z0;
        ABCD_Y = [1, 0; Ys, 1];

        % Air TL
        ABCD_TL_air = [cos(bdAir), 1i * Z0 * sin(bdAir); 1i * sin(bdAir)/Z0, cos(bdAir)];
        
        % Substrate TL
        ABCD_TL_sub = [cos(bdSub), 1i * Zd * sin(bdSub); 1i * sin(bdSub)/Zd, cos(bdSub)];

        ABCD = ABCD_Y * ABCD_TL_sub * ABCD_TL_air * ABCD_TL_sub * ABCD_Y;

        A = ABCD(1,1);
        B = ABCD(1,2);
        C = ABCD(2,1);
        D = ABCD(2,2);
        S21 = 2 / (A + B/Z0 + C * Z0 + D);
        %S21s(i, point) = 20 * log10(abs(S21));
        S21s(i, point) = abs(S21);
    end
    S21SingleFreq = S21s(i, :);
    below = S21SingleFreq(yns < 0);
    above = S21SingleFreq(~(yns < 0));
    [~, idxBelow] = max(below);
    bestYsBelow(i) = yns(idxBelow);
    bestZsBelow(i) = -1 ./ (bestYsBelow(i) / 377);
    [~, idxAbove] = max(above);
    bestYsAbove(i) = yns(idxAbove + numel(below));
    bestZsAbove(i) = -1 ./ (bestYsAbove(i) / 377);
end

%% Plot the surfaces and ideal impedances.
if enablePlot
    figure('Name', 'Tansmission vs Layer');
    subplot(2, 1, 1)
    hold on
    [Y, F] = meshgrid(yns, frequencies);
    S21s(S21s < 0.9) = 0;
    h = surf(Y, F, S21s);
    set(h, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
    zlim([0.1, 1])
    %colormap('jet')
    %view(3);
    subplot(2, 1, 2)
    hold on
    plot(frequencies, bestZsBelow);
    plot(frequencies, bestZsAbove);
end

requiredImpedancesBelow = bestZsBelow;
requiredImpedancesAbove = bestZsAbove;

end

