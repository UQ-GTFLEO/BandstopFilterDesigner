function [L, C1, C2] = get_symbolic_impedances(C1_reduction_factor)
% Retrieve all the required symbolic variables

% Physical Constants
% Define some physical constants that we will use throughout without change
mu0 = 4 * pi * 10^-7;
e0 = 8.8541878188e-12;
D = 10e-3;

c1_alpha = 2 * D * e0 / pi;
C1 = @(w) -c1_alpha * C1_reduction_factor * log(sin(pi * w / (2 * D)));

% SECOND COMPONENT PAIR - C2 and L
% * Run the code from the Predictive model.
[startSym, endSym] = getS11_symbolic();
S11 = @(wp, wm) (startSym(wm) + (endSym(wm) - startSym(wm)) * wp / wm);

% C2 Prediction
C2_Ratio_Start = @(wp) 0.2178 * (wp * 1e3)^2 - 0.9642 * (wp * 1e3) + 3.8577;
C2_Ratio_End = @(wp) 10 / sqrt(wp * 1e3);
C2_Ratio = @(wp, wm) C2_Ratio_Start(wp) + (5e-3 - wm) / (5e-3 - wp) * (C2_Ratio_End(wp) - C2_Ratio_Start(wp));
C2 = @(wp, wm) C1(wp) / C2_Ratio(wp, wm);

% Infer L
L = @(wp, wm) 1 / ((S11(wp, wm) * 2 * pi)^2 * C2(wp, wm));
end