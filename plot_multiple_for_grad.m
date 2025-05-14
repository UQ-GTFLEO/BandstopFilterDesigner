%figure;

w_meshes    = [ones(1, 24) * 2.5, ones(1, 19) * 2, ones(1, 14) * 1.5, ones(1, 9), ones(1, 4) * 0.5];
w_patches   = [linspace(0.1, 2.4, 24), linspace(0.1, 1.9, 19), linspace(0.1, 1.4, 14), linspace(0.1, 0.9, 9), linspace(0.1, 0.4, 4)];

w_meshes = [2.2];
w_patches = [2];

for idx = 1
    w_m = w_meshes(idx) * 1e-3;
    w_p = w_patches(idx) * 1e-3;
    period = 5;
    add_dim_plot(w_m, w_p, period)

    pause;
end