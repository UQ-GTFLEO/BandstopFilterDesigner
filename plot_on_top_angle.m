indexat = @(expr, index) expr(index);
FreqMin = 0;
FreqMax = 40;
eps_r = 3.55;
angMax = FreqMax * 2 * pi;
[LSym, C1Sym, C2Sym] = get_symbolic_impedances(1);

w_mesh = 4.4e-3;
w_patch = 4e-3;

period = 5;
scale_factor = period / 10;
ratio = 0.35;
dielectric_factor = (eps_r * ratio) + 1 * (1 - ratio);

% Scale spacing
dSub = 0.508 * 1e-3;
dAir = 3.2 * 1e-3 - (dSub * 2);

% Base Case Air
lambdaAir = @(f) physconst('LightSpeed') ./ f;
dWavelengthsAir = @(f) dAir ./ lambdaAir(f);
bdAir = @(f) 2 * pi * dWavelengthsAir(f);

% Substrate Case
lambdaSub = @(f) lambdaAir(f) ./ sqrt(eps_r);
dSubWavelengths = @(f) dSub ./ lambdaSub(f);
bdSub = @(f) 2 * pi * dSubWavelengths(f);

% Free space and substrate impedances
Zd = 377/sqrt(eps_r);
Z0 = 377;


ABCD_TL_air = @(f) [cos(bdAir(f)), 1i * Z0 * sin(bdAir(f)); 1i * sin(bdAir(f))/Z0, cos(bdAir(f))];
        
% Substrate TL
ABCD_TL_sub = @(f) [cos(bdSub(f)), 1i * Zd * sin(bdSub(f)); 1i * sin(bdSub(f))/Zd, cos(bdSub(f))];

C1 = dielectric_factor * scale_factor * C1Sym(w_patch);
L = scale_factor * LSym(w_patch, w_mesh);
C2 = dielectric_factor * scale_factor * C2Sym(w_patch, w_mesh);

Z_C1 = @(f) 1 / (1i * 2 * pi * f * C1);
Z_C2 = @(f) 1 / (1i * 2 * pi * f * C2);
Z_L = @(f) 1i * 2 * pi * f * L;

Z0 = 377;

A_Y = 1;
B_Y = 0;
Z = @(f) 1 / (1i * 2 * pi * f * C1 * 1e9)...
    + (1i * 2 * pi * f * C2 * 1e9 + 1 / (1i * 2 * pi * f * 1e9 * L))^-1;
C_Y = @(f) 1 / Z(f);
D_Y = 1;
ABCD_Y = @(f) [A_Y, B_Y; C_Y(f), D_Y];

ABCD = @(f) ABCD_Y(f) * ABCD_TL_sub(f  * 1e9) * ABCD_TL_air(f * 1e9) * ABCD_TL_sub(f * 1e9) * ABCD_Y(f);

A = @(f) indexat(ABCD(f), 1);
B = @(f) indexat(ABCD(f), 3);
C = @(f) indexat(ABCD(f), 2);
D = @(f) indexat(ABCD(f), 4);

S11 = @(f) 20 * log10(abs((A(f) + B(f)/Z0 - C(f) * Z0 - D(f)) / (A(f) + B(f)/Z0 + C(f) * Z0 + D(f))));
S21 = @(f) 20 * log10(abs(2 / (A(f) + B(f)/Z0 + C(f) * Z0 + D(f))));

fplot(S21, [0, FreqMax], 'Color', 'blue', 'LineWidth', 2, 'DisplayName', 'Model');



