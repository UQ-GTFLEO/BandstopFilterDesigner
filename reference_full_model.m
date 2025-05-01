%% OVERVIEW
% We have 3 components to be able to find.
% They are C1, C2 and L

%% Physical Constants
% Define some physical constants that we will use throughout without change
clear variables;
clc;
mu0 = 4 * pi * 10^-7;
e0 = 8.8541878188e-12;
D = 10e-3;


%% SUBJECT TO CHANGE
% * May change the averaged mode apporach to that of the geometric form
% from the other earlier papers
% * May replace the asssumption for C1 with the curves that we found of S21

%% FIRST COMPONENT - C1
% * We assume that the capacitance, C1, follows that of just a patch
% * Refer to: An Overview of Equivalent Circuit Modeling Techniques of
% Frequency Selective Surfaces and Metasurfaces, Costa 2014

c1_alpha = 2 * D * e0 / pi;
C1 = @(w) -c1_alpha * log(sin(pi * w / (2 * D)));
% subplot(3, 2, 1)
% fplot(C1);
% xlim([0 D]);
% title('C1 Averaged Prediction')
% xlabel('W PATCH')
% ylabel('C')



%% SECOND COMPONENT PAIR - C2 and L
% * Run the code from the Predictive model.
[startSym, endSym] = f_pred_s11();
S11_Sym = @(wp, wm) (startSym(wm) + (endSym(wm) - startSym(wm)) * wp / wm);


%% C2 Prediction
C2_Ratio_Start = @(wp) 0.2178 * (wp * 1e3)^2 - 0.9642 * (wp * 1e3) + 3.8577;
%C2_Ratio_Start = @(wp) 0.325 * (wp * 1e3);
C2_Ratio_End = @(wp) 10 / sqrt(wp * 1e3);
C2_Ratio = @(wp, wm) C2_Ratio_Start(wp) + (5e-3 - wm) / (5e-3 - wp) * (C2_Ratio_End(wp) - C2_Ratio_Start(wp));
%C2_Ratio = @(wp, wm) C2_Ratio_Start(wp) + (10e-3 - wm) / (10e-3 - wp) * (C2_Ratio_End(wp) - C2_Ratio_Start(wp));
C2 = @(wp, wm) C1(wp) / C2_Ratio(wp, wm);

%% Infer L
L_Sym = @(wp, wm) 1 / ((S11_Sym(wp, wm) * 2 * pi)^2 * C2(wp, wm));

%% Infer S21
S21_sym = @(wp, wm) 1 / (sqrt(L_Sym(wp, wm) * (C1(wp) + C2(wp, wm))) * 2 * pi);

%% Add a function of the Impedance
freq = 10e9;
ang_freq = freq * 2 * pi;
%Z_sym = @(wp, wm) imag(1/(1i * ang_freq * C1(wp)) + (1i * ang_freq * C2(wp, wm) + 1/(1i * ang_freq * L_Sym(wp, wm)))^-1);
Z_sym = @(wp, wm) imag((1 - ang_freq^2 * L_Sym(wp, wm) * (C1(wp) + C2(wp, wm))) / (1i * ang_freq * C1(wp) * (1 - ang_freq^2 * L_Sym(wp, wm) * C2(wp, wm))));



%% Plotting
for w_patch = 0.1e-3:0.1e-3:5e-3
    L1s = [];
    C2s = [];
    S21s = [];
    S11s = [];
    Zs = [];
    C2_ratios = [];
    w_mesh = (w_patch + 0.1e-3):0.1e-3:5e-3;
    s = size(w_mesh);
    s = s(:,2);
    for ws = w_mesh
        L_tmp = L_Sym(w_patch, ws) * 1e9;
        L1s = [L1s, L_tmp];

        s11_tmp = S11_Sym(w_patch, ws);
        S11s = [S11s, s11_tmp];

        C2_Ratios_tmp = C2_Ratio(w_patch, ws);
        C2_ratios = [C2_ratios, C2_Ratios_tmp];

        C2_tmp = C2(w_patch, ws) * 1e15;
        C2s = [C2s C2_tmp];

        s21_tmp = S21_sym(w_patch, ws);
        S21s = [S21s s21_tmp];

        Z_tmp = Z_sym(w_patch, ws);
        Zs = [Zs, Z_tmp];

    end
    ys = linspace(w_patch, w_patch, s);
    %subplot(3, 2, 1)
    subplot(2, 3, 1)
    plot3(w_mesh, ys, L1s);
    title('L Result')
    xlabel('W MESH')
    ylabel('W PATCH')
    hold on
    %subplot(3, 2, 2)
    subplot(2, 3, 2)
    plot3(w_mesh, ys, C2s);
    title('C2 Result')
    xlabel('W MESH')
    ylabel('W PATCH')
    hold on
    %subplot(3, 2, 3)
    subplot(2, 3, 3)
    plot3(w_mesh, ys, C2_ratios);
    title('C2 Ratio')
    xlabel('W MESH')
    ylabel('W PATCH')
    hold on
    %subplot(3, 2, 4)
    subplot(2, 3, 4)
    plot3(w_mesh, ys, S21s);
    title('S21 Prediction')
    xlabel('W MESH')
    ylabel('W PATCH')
    hold on
    %subplot(3, 2, 5)
    subplot(2, 3, 5)
    plot3(w_mesh, ys, S11s);
    title('S11 Prediction')
    xlabel('W MESH')
    ylabel('W PATCH')
    hold on
    %subplot(3, 2, 6)
    subplot(2, 3, 6)
    plot3(w_mesh, ys, Zs);
    title('Impedance Prediction')
    xlabel('W MESH')
    ylabel('W PATCH')
    hold on

end








