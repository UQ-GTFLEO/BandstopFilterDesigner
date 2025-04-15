clear variables;
clc;
current_cs = 5;
freq = 20e9;

red_factor = 1;
Z_sym = scale_sym_vars(current_cs, freq, 1, 4.38, red_factor);
zero_dims = [];
w_patchs = 0.1e-3:0.1e-3:5e-3;

for w_patch = w_patchs
    Yn1 = [];
    Yn2 = [];
    Z1s = [];
    Z2s = [];

    w_mesh = (w_patch + 0.1e-3):0.1e-3:5e-3;
    s = size(w_mesh);
    s = s(:,2);
    for ws = w_mesh
        Z_tmp = Z_sym(w_patch, ws);
        Z1s = [Z1s, Z_tmp];
    end
    adj_patch = 1e3 * w_patch * current_cs / 10;
    ys = linspace(adj_patch, adj_patch, s);
    adj_mesh = 1e3 * w_mesh * current_cs / 10;

    scatter3(adj_mesh, ys, Z1s, 'b', 'x');
    title('Normalised Impedance at 20GHz')
    xlabel('w_{m}')
    ylabel('w_{p}')
    hold on

end



