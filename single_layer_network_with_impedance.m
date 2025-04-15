clear variables;
clc;
indexat = @(expr, index) expr(index);

FreqMin = 1;
FreqMax = 40;
eps_r = 3.55;
angMax = FreqMax * 2 * pi;
[LSym, C1Sym, C2Sym] = get_symbolic_impedances(1);

w_mesh = 4e-3;
w_patch = 3.2e-3;

period = 5;
scale_factor = period / 10;
dielectric_factor = (eps_r * 0.81 + 1) / 2;

C1 = dielectric_factor * scale_factor * C1Sym(w_patch);
L = scale_factor * LSym(w_patch, w_mesh);
C2 = dielectric_factor * scale_factor * C2Sym(w_patch, w_mesh);

Z_C1 = @(f) 1 / (1i * 2 * pi * f * C1);
Z_C2 = @(f) 1 / (1i * 2 * pi * f * C2);
Z_L = @(f) 1i * 2 * pi * f * L;

% BASE VALUES WE INPUT
ABCD_C1 = @(f) [1, Z_C1(f); 0, 1];
ABCD_C2 = @(f) [1, Z_C2(f); 0, 1];
Z_MAT_L = @(f) [Z_L(f), Z_L(f); Z_L(f), Z_L(f)];
ABCD_L = @(f) z2abcd(Z_MAT_L(f));

% CONVERT PAIR
Y_L = @(f) abcd2y(ABCD_L(f));
Y_C2 = @(f) abcd2y(ABCD_C2(f));
Y_PAIR = @(f) Y_L(f) + Y_C2(f);
ABCD_PAIR = @(f) y2abcd(Y_PAIR(f));
ABCD_LAYER = @(f) ABCD_C1(f) * ABCD_PAIR(f);

Z0 = 377;

A = 1;
B = 0;
Z = @(f) 1 / (1i * 2 * pi * f * C1 * 1e9)...
    + (1i * 2 * pi * f * C2 * 1e9 + 1 / (1i * 2 * pi * f * 1e9 * L))^-1;
C = @(f) 1 / Z(f);
D = 1;

S11 = @(f) 20 * log10(abs((A + B/Z0 - C(f) * Z0 - D) / (A + B/Z0 + C(f) * Z0 + D)));
S21 = @(f) 20 * log10(abs(2 / (A + B/Z0 + C(f) * Z0 + D)));

subplot(2, 1, 1)
fplot(S21, [FreqMin, FreqMax]);
hold on 
ylim([-60, 0])
xline(1e-9 / (2 * pi * sqrt(L * (C1 + C2))));
xlim([0 FreqMax])

subplot(2, 1, 2)
hold on
Yn_plot = @(f) imag(1 / Z(f)) * Z0;
fplot(Yn_plot, [0, 15], 'b');
fplot(Yn_plot, [30, 40], 'b');
ylim([-10, 10]);


