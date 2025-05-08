indexat = @(expr, index) expr(index);
FreqMin = 0;
FreqMax = 40;
eps_r = 3.55;
angMax = FreqMax * 2 * pi;
[LSym, C1Sym, C2Sym] = get_symbolic_impedances(0.9);

w_mesh = 4e-3;
w_patch = 3.2e-3;

period = 5;
scale_factor = period / 10;
ratio = 0.35;
dielectric_factor = (eps_r * ratio) + 1 * (1 - ratio);

C1 = dielectric_factor * scale_factor * C1Sym(w_patch);
L = scale_factor * LSym(w_patch, w_mesh);
C2 = dielectric_factor * scale_factor * C2Sym(w_patch, w_mesh);

Z_C1 = @(f) 1 / (1i * 2 * pi * f * C1);
Z_C2 = @(f) 1 / (1i * 2 * pi * f * C2);
Z_L = @(f) 1i * 2 * pi * f * L;

Z0 = 377;

A = 1;
B = 0;
Z = @(f) 1 / (1i * 2 * pi * f * C1 * 1e9)...
    + (1i * 2 * pi * f * C2 * 1e9 + 1 / (1i * 2 * pi * f * 1e9 * L))^-1;
C = @(f) 1 / Z(f);
D = 1;

S11 = @(f) 20 * log10(abs((A + B/Z0 - C(f) * Z0 - D) / (A + B/Z0 + C(f) * Z0 + D)));
S21 = @(f) 20 * log10(abs(2 / (A + B/Z0 + C(f) * Z0 + D)));


%subplot(1, 2, 1)
hold on
Yn_plot = @(f) imag(1 / Z(f)) * Z0;
fplot(Yn_plot, [FreqMin, FreqMax], 'b', 'LineWidth', 2);
ylim([-5, 5]);
xline(1e-9 / (2 * pi * sqrt(L * (C1 + C2))));
%xline(22.5, 'LineWidth',2);

% subplot(1, 2, 2)
% sParams = sparameters("NO_SUB_one_layer.s4p");
% fplot(S21, [FreqMin, FreqMax], 'b', 'LineWidth', 2);
% hold on 
% ylim([-60, 0])
% xline(1e-9 / (2 * pi * sqrt(L * (C1 + C2))));
% xlim([0 FreqMax])
% plot_s_params = rfplot(sParams, 1, 3);
% set(plot_s_params, 'Color', 'black', 'LineWidth', 2)

