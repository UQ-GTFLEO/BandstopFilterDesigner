clear variables;
clc;

total_samples = 1002;
FreqMin = 0;
FreqMax = 40;
eps_r = 2.2;
d = 0.508 * 10^-3;
%d = 0.787 * 10^-3;
%d = 1.575 * 10^-3;
period = 5;
w_meshes    = [ones(1, 24) * 2.5, ones(1, 19) * 2, ones(1, 14) * 1.5, ones(1, 9), ones(1, 4) * 0.5];
w_patches   = [linspace(0.1, 2.4, 24), linspace(0.1, 1.9, 19), linspace(0.1, 1.4, 14), linspace(0.1, 0.9, 9), linspace(0.1, 0.4, 4)];

figure('Name', "Three Layer - Coupling and No Coupling")

for tmp = 1:70
    filename = "./unit_cell/508_unit_cell_transmitarray_" + tmp + ".s4p";
    %filename = "./unit_cell/787_unit_cell_transmitarray_" + tmp + ".s4p";
    %filename = "./unit_cell/1575_unit_cell_transmitarray_" + tmp + ".s4p";
    sParams = sparameters(filename);
    frequencies = sParams.Frequencies;
    freqsGHz = linspace(FreqMin, FreqMax, total_samples);
    km = 0.4;
    w_mesh_outside = w_meshes(tmp) * 2e-3;
    w_mesh_inside = w_meshes(tmp) * 2e-3;
    w_patch_outside = w_patches(tmp) * 2e-3;
    w_patch_inside = w_patches(tmp) * 2e-3;
    
    [S21_no_coupling ,S21_coupling] = three_layer_prediction_model_comparison(d, eps_r, km, period, w_mesh_outside, w_mesh_inside, w_patch_outside, w_patch_inside);
    S21_no_coupling_db = @(f) 10 * log(abs(S21_no_coupling(f)));
    S21_no_coupling_angle = @(f) mod(angle(S21_no_coupling(f)) * 180 / pi, -360);
    S21_coupling_db = @(f) 10 * log(abs(S21_coupling(f)));
    S21_coupling_angle = @(f) mod(angle(S21_coupling(f)) * 180 / pi, -360);

    S21_no_coupling_dbs = zeros(total_samples - 1, 1);
    S21_no_coupling_angles = zeros(total_samples - 1, 1);
    S21_coupling_dbs = zeros(total_samples - 1, 1);
    S21_coupling_angles = zeros(total_samples - 1, 1);
    freqsGHz = freqsGHz(2:end);
    for i = 1:(total_samples - 1)
        freq = freqsGHz(i);
        S21_no_coupling_dbs(i) = S21_no_coupling_db(freq);
        S21_no_coupling_angles(i) = S21_no_coupling_angle(freq);
        S21_coupling_dbs(i) = S21_coupling_db(freq);
        S21_coupling_angles(i) = S21_coupling_angle(freq);
    end

    subplot(7, 10, tmp)
    %plot(freqsGHz, S21_coupling_dbs, 'b', 'LineWidth', 2);
    hold on
    %plot(freqsGHz, S21_no_coupling_dbs, 'g', 'LineWidth', 2);
    
    ylim([-80, 0])
    plot_s_params = rfplot(sParams, 1, 3);
    set(plot_s_params, 'Color', 'black', 'LineWidth', 2)

    yyaxis right
    %plot(freqsGHz, S21_coupling_angles, "--", "Color", 'blue', 'LineWidth', 1);
    %plot(freqsGHz, S21_no_coupling_angles, "--", "Color", 'green', 'LineWidth', 1);
    data = rfparam(sParams,1,3);
    wrap = @(x) mod(180*angle(x)/pi, -360);
    plot(sParams.Frequencies / 10^9,wrap(data), "--", "Color", 'black', 'LineWidth', 1);
    ylim([-360, 0])
    legend('hide')
    title("w_{m} = " + w_meshes(tmp) + ", w_{p} = " + w_patches(tmp))

end