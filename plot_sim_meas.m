%function [eps_r] = get_exact_eps_r()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


indexat = @(expr, index) expr(index);
FreqMin = 0;
FreqMax = 40;
eps_r_layer = 3.55;
angMax = FreqMax * 2 * pi;
[LSym, C1Sym, C2Sym] = get_symbolic_impedances(0.9);

w_mesh = 4e-3;
w_patch = 3.2e-3;

period = 5;
scale_factor = period / 10;
ratio = 0.35;
dielectric_factor = (eps_r_layer * ratio) + 1 * (1 - ratio);

C1 = dielectric_factor * scale_factor * C1Sym(w_patch);
L = scale_factor * LSym(w_patch, w_mesh);
C2 = dielectric_factor * scale_factor * C2Sym(w_patch, w_mesh);

Z0 = 377;


Z = @(f) 1 / (1i * 2 * pi * f * C1 * 1e9)...
    + (1i * 2 * pi * f * C2 * 1e9 + 1 / (1i * 2 * pi * f * 1e9 * L))^-1;

C = @(f) 1 / Z(f);
ABCD_metal = @(f) [1, 0; C(f), 1];

% Base Case Air
lambdaAir = @(f) physconst('LightSpeed') ./ f;
dAir = 3.5 * 1e-3 - 2 * 0.508 * 1e-3;
dAirWavelengths = @(f) dAir ./ lambdaAir(f);
bdAir = @(f) 2 * pi * dAirWavelengths(f);
Z0 = 377;

% Substrate Case
dSub = 0.508 * 1e-3;
lambdaSub = @(f) lambdaAir(f) ./ sqrt(eps_r_layer);
dSubWavelengths = @(f) dSub ./ lambdaSub(f);
bdSub = @(f) 2 * pi * dSubWavelengths(f);
Zd = 377/sqrt(eps_r_layer);

ABCD_TL_sub = @(f) [cos(bdSub(f)), 1i * Zd * sin(bdSub(f)); 1i * sin(bdSub(f))/Zd, cos(bdSub(f))];
ABCD_TL_air = @(f) [cos(bdAir(f)), 1i * Z0 * sin(bdAir(f)); 1i * sin(bdAir(f))/Z0, cos(bdAir(f))];

ABCD = @(f) ABCD_metal(f) * ABCD_TL_sub(f * 10^9) * ABCD_TL_air(f * 10^9) * ABCD_TL_sub(f * 10^9) * ABCD_metal(f);

A = @(f) indexat(ABCD(f), 1);
B = @(f) indexat(ABCD(f), 3);
C = @(f) indexat(ABCD(f), 2);
D = @(f) indexat(ABCD(f), 4);

S21 = @(f) 20 * log10(abs(2 / (A(f) + B(f)/Z0 + C(f) * Z0 + D(f))));

figure;
hold on;
fplot(S21, [FreqMin, FreqMax], 'r', 'LineWidth', 2)
sParams = sparameters("HeightComparison/rogers4003c_0508_6.s4p");
plot_s = rfplot(sParams, 1, 3);
set(plot_s, 'Color', 'black', 'LineWidth', 2)
xlim([FreqMin FreqMax])
