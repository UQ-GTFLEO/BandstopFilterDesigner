%function [eps_r] = get_exact_eps_r()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


indexat = @(expr, index) expr(index);
FreqMin = 0;
FreqMax = 60;
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

% Substrate Case
dSub = 0.508 * 1e-3;
lambdaSub = @(f) lambdaAir(f) ./ sqrt(eps_r_layer);
dSubWavelengths = @(f) dSub ./ lambdaSub(f);
bdSub = @(f) 2 * pi * dSubWavelengths(f);
Zd = 377/sqrt(eps_r_layer);

ABCD_TL_sub = @(f) [cos(bdSub(f)), 1i * Zd * sin(bdSub(f)); 1i * sin(bdSub(f))/Zd, cos(bdSub(f))];
ABCD = @(f) ABCD_metal(f) * ABCD_TL_sub(f * 10^9);

A = @(f) indexat(ABCD(f), 1);
B = @(f) indexat(ABCD(f), 3);
C = @(f) indexat(ABCD(f), 2);
D = @(f) indexat(ABCD(f), 4);

S21 = @(f) 20 * log10(abs(2 / (A(f) + B(f)/Z0 + C(f) * Z0 + D(f))));

figure;
hold on;
fplot(S21, [FreqMin, FreqMax], 'r', 'LineWidth', 2)
sParams_sub = sparameters("rogers4003c_one_layer_for_sub_weight.s4p");
sParams_no_sub = sparameters("NO_SUB_one_layer.s4p");
plot_s_params_sub = rfplot(sParams_sub, 1, 3);
plot_s_params_no_sub = rfplot(sParams_no_sub, 1, 3);
set(plot_s_params_sub, 'Color', 'black', 'LineWidth', 2)
set(plot_s_params_no_sub, 'Color', 'blue', 'LineWidth', 2)

eps_r = 1;
%end