clear variables;
clc;
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

bd = @(f) bdSub(f) + bdAir(f) + bdSub(f);

ABCD_TL_air = @(f) [cos(bdAir(f)), 1i * Z0 * sin(bdAir(f)); 1i * sin(bdAir(f))/Z0, cos(bdAir(f))];
        
% Substrate TL
ABCD_TL_sub = @(f) [cos(bdSub(f)), 1i * Zd * sin(bdSub(f)); 1i * sin(bdSub(f))/Zd, cos(bdSub(f))];

ABCD_TL = @(f) ABCD_TL_sub(f) * ABCD_TL_air(f) * ABCD_TL_sub(f);

C1 = dielectric_factor * scale_factor * C1Sym(w_patch);
L = scale_factor * LSym(w_patch, w_mesh);
C2 = dielectric_factor * scale_factor * C2Sym(w_patch, w_mesh);


figure;


% Choose file %
sParams = sparameters("./inline.s4p");
frequencies = sParams.Frequencies;
freqsGHz = linspace(0, 40, 1002);

km = -0.01;
ke1 = 0;
ke2 = 0.05;


%%%% INCDUCTIVE COUPLING AMOUNT %%%
Lm = L * km;
Cm1 = C1 * ke1;
Cm2 = C2 * ke2;

Z_C1 = @(f) 1/(1i * 2 * pi * f * (C1 - Cm1)); % Keep as is

Z_C2e = @(f) 1/(1i * 2 * pi * f * (C2 - Cm2));
Z_C2o = @(f) 1/(1i * 2 * pi * f * (C2 - Cm2 + 2 * Cm1));

Z_Le = @(f) 1i * 2 * pi * f * (L + Lm); % Keep as is
Z_Lo = @(f) 1i * 2 * pi * f * (L - Lm); % Keep as is

Y_Combined = @(f) abcd2y(ABCD_TL(f));

Y11_TL = @(f) indexat(Y_Combined(f), 1); %@(f) -1i * cot(bd(f)) / Zd;
Y12_TL = @(f) indexat(Y_Combined(f), 2); %1i * csc(bd(f)) / Zd;
Y_TL = @(f) -Y12_TL(f) * 2;
Y_TC = @(f) Y11_TL(f) + Y12_TL(f);

Y0 = 1 / 377;
Ye = @(f) 1 / ( Z_C1(f) + (Z_C2e(f) * Z_Le(f))/(Z_C2e(f) + Z_Le(f))) + Y_TC(f);
Yo = @(f) 1 / (Z_C1(f) + (Z_C2o(f) * Z_Lo(f))/(Z_C2o(f) + Z_Lo(f))) + Y_TC(f) + Y_TL(f);
S21 = @(f) 20 * log10(abs((Yo(f * 10^9) * Y0 - Ye(f * 10^9) * Y0) / ((Y0 + Ye(f * 10^9)) * (Y0 + Yo(f * 10^9)))));
S11 = @(f) 20 * log10(abs((Y0^2 - Ye(f * 10^9) * Yo(f * 10^9)) / ((Y0 + Ye(f * 10^9)) * (Y0 + Yo(f * 10^9)))));
S21_angle = @(f) mod(angle((Yo(f * 10^9) * Y0 - Ye(f * 10^9) * Y0) / ((Y0 + Ye(f * 10^9)) * (Y0 + Yo(f * 10^9)))) * 180/pi, -360);
S11_angle = @(f) mod(angle((Y0^2 - Ye(f * 10^9) * Yo(f * 10^9)) / ((Y0 + Ye(f * 10^9)) * (Y0 + Yo(f * 10^9)))) * 180/pi, -360);

fplot(S21, [0,  FreqMax] , 'b', 'LineWidth', 2);
hold on
ylim([-80, 0])
plot_s_params = rfplot(sParams, 1, 3);
set(plot_s_params, 'Color', 'black', 'LineWidth', 2)
legend('hide')
title("w_{m} = " + w_mesh + ", w_{p} = " + w_patch)
xlim([0, FreqMax])

% yyaxis right
% data = rfparam(sParams,1,3);
% wrap = @(x) mod(180*angle(x)/pi, -360);
% plot(sParams.Frequencies / 10^9, wrap(data), "--", "Color", 'black', 'LineWidth', 1);
% fplot(S21_angle, [0,  FreqMax] , "--", "Color", 'blue', 'LineWidth', 1);
% ylim([-360 0]) 

