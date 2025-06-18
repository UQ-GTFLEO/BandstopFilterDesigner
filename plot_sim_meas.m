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

measuredFreq = [3.95	4.05	4.15	4.25	4.35	4.45	4.55	4.65	4.75	4.85	4.95	5.05	5.15	5.25	5.35	5.45	5.55	5.65	5.75	5.85	5.95	6.05	6.15	6.25	6.35	6.45	6.55	6.65	6.75	6.85	6.95	7.05	7.15	7.25	7.35	7.45	7.55	7.65	7.75	7.85	7.95	8.05	8.15	8.2	8.3	8.4	8.5	8.6	8.7	8.8	8.9	9	9.1	9.2	9.3	9.4	9.5	9.6	9.7	9.8	9.9	10	10.1	10.2	10.3	10.4	10.5	10.6	10.7	10.8	10.9	11	11.1	11.2	11.3	11.4	11.5	11.6	11.7	11.8	11.9	12	12.1	12.2	12.3	12.4	18	18.1	18.2	18.3	18.4	18.5	18.6	18.7	18.8	18.9	19	19.1	19.2	19.3	19.4	19.5	19.6	19.7	19.8	19.9	20	20.1	20.2	20.3	20.4	20.5	20.6	20.7	20.8	20.9	21	21.1	21.2	21.3	21.4	21.5	21.6	21.7	21.8	21.9	22	22.1	22.2	22.3	22.4	22.5	22.6	22.7	22.8	22.9	23	23.1	23.2	23.3	23.4	23.5	23.6	23.7	23.8	23.9	24	24.1	24.2	24.3	24.4	24.5	24.6	24.7	24.8	24.9	25	25.1	25.2	25.3	25.4	25.5	25.6	25.7	25.8	25.9	26	26.1	26.2	26.3	26.4	26.5];
measuredS21 = [-0.368	-0.23	-0.407	-0.312	-0.445	-0.228	-0.4	-0.193	-0.112	0.17	-0.253	-0.095	-0.474	-0.467	-0.794	-0.706	-1.047	-0.919	-0.954	-0.364	-0.558	-0.293	-0.578	-0.481	-0.714	-0.822	-1.199	-0.997	-1.613	-1.163	-1.681	-1.088	-1.358	-0.734	-1.215	-0.743	-1.063	-0.767	-0.861	-1.261	-0.995	-1.205	-0.889	-1.614	-1.837	-1.47	-1.342	-1.098	-0.612	-0.619	-0.161	-0.38	-0.075	-0.467	-0.484	-0.756	-0.885	-1.327	-1.112	-1.51	-1.023	-1.356	-0.762	-0.945	-0.362	-0.447	-0.082	-0.322	-0.208	-0.603	-0.779	-1.026	-1.162	-1.068	-1.052	-0.645	-0.73	-0.353	-0.375	-0.285	-0.327	-0.454	-0.326	-0.614	-0.426	-0.749	-9.105	-9.643	-10.621	-11.305	-11.395	-11.876	-12.798	-13.926	-14.991	-14.863	-14.569	-14.757	-15.258	-16.372	-17.589	-18.836	-19.279	-19.206	-19.933	-20.89	-21.877	-22.694	-24.399	-28.365	-30.799	-30.488	-29.857	-28.756	-29.24	-36.24	-34.081	-33.754	-33.256	-33.52	-35.54	-32.782	-34.29	-34.224	-32.312	-30.621	-29.152	-28.549	-28.234	-28.067	-28.428	-27.782	-27.922	-29.333	-31.226	-34.471	-36.884	-36.302	-35.61	-36.21	-33.932	-32.689	-32.57	-31.855	-28.234	-26.628	-25.536	-23.501	-21.009	-18.67	-18.781	-18.924	-17.983	-15.334	-13.105	-13.497	-13.47	-12.438	-9.99	-8.111	-9.034	-9.174	-8.244	-6.228	-4.391	-4.401	-4.408	-3.643	-2.19	-0.785	-0.899	-1.314];




figure;
hold on;
scatter(measuredFreq, measuredS21, 'b', 'x')
%fplot(S21, [FreqMin, FreqMax], 'r', 'LineWidth', 2)
sParams = sparameters("final_design3_4_5.s4p");
plot_s = rfplot(sParams, 1, 3);
set(plot_s, 'Color', 'black', 'LineWidth', 2)
xlim([FreqMin FreqMax])
