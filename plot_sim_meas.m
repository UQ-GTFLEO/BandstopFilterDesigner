%function [eps_r] = get_exact_eps_r()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

indexat = @(expr, index) expr(index);
FreqMin = 0;
FreqMax = 40;
eps_r_layer = 3.55;
angMax = FreqMax * 2 * pi;
[LSym, C1Sym, C2Sym] = get_symbolic_impedances(1);

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

measuredFreq = [8.2	8.3	8.4	8.5	8.6	8.7	8.8	8.9	9	9.1	9.2	9.3	9.4	9.5	9.6	9.7	9.8	9.9	10	10.1	10.2	10.3	10.4	10.5	10.6	10.7	10.8	10.9	11	11.1	11.2	11.3	11.4	11.5	11.6	11.7	11.8	11.9	12	12.1	12.2	12.3	12.4	18	18.1	18.2	18.3	18.4	18.5	18.6	18.7	18.8	18.9	19	19.1	19.2	19.3	19.4	19.5	19.6	19.7	19.8	19.9	20	20.1	20.2	20.3	20.4	20.5	20.6	20.7	20.8	20.9	21	21.1	21.2	21.3	21.4	21.5	21.6	21.7	21.8	21.9	22	22.1	22.2	22.3	22.4	22.5	22.6	22.7	22.8	22.9	23	23.1	23.2	23.3	23.4	23.5	23.6	23.7	23.8	23.9	24	24.1	24.2	24.3	24.4	24.5	24.6	24.7	24.8	24.9	25	25.1	25.2	25.3	25.4	25.5	25.6	25.7	25.8	25.9	26	26.1	26.2	26.3	26.4	26.5];
measuredS21 = [-1.614	-1.837	-1.47	-1.342	-1.098	-0.612	-0.619	-0.161	-0.38	-0.075	-0.467	-0.484	-0.756	-0.885	-1.327	-1.112	-1.51	-1.023	-1.356	-0.762	-0.945	-0.362	-0.447	-0.082	-0.322	-0.208	-0.603	-0.779	-1.026	-1.162	-1.068	-1.052	-0.645	-0.73	-0.353	-0.375	-0.285	-0.327	-0.454	-0.326	-0.614	-0.426	-0.749	-10.336	-10.646	-11.387	-12.28	-13.5	-14.296	-14.697	-15.352	-15.548	-15.776	-15.736	-15.771	-15.682	-15.5	-15.502	-15.927	-16.452	-17.574	-19.235	-19.991	-19.996	-19.713	-19.899	-19.462	-20.164	-21.167	-21.027	-20.942	-19.532	-20.489	-21.524	-22.914	-23.544	-25.178	-25.189	-25.534	-27.198	-25.958	-24.956	-26.648	-28.642	-27.328	-26.909	-24.037	-23.097	-23.778	-25.993	-26.827	-26.996	-24.22	-22.55	-22.266	-22.892	-23.136	-23.567	-21.963	-19.448	-18.191	-17.506	-17.813	-18.435	-17.013	-15.347	-13.853	-13.33	-13.109	-13.568	-13.344	-12.243	-11.312	-10.406	-9.795	-9.421	-9.031	-8.155	-7.021	-6.041	-5.07	-4.675	-4.01	-3.754	-3.092	-2.631	-1.961	-1.154	-0.627];




figure;
hold on;
scatter(measuredFreq, measuredS21, 'b', 'x')
%fplot(S21, [FreqMin, FreqMax], 'r', 'LineWidth', 2)
sParams = sparameters("final_design3_4_5.s4p");
plot_s = rfplot(sParams, 1, 3);
set(plot_s, 'Color', 'black', 'LineWidth', 2)
xlim([FreqMin FreqMax])
