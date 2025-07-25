function [startSym,endSym] = getS11_symbolic()

mu0 = 4 * pi * 10^-7;
e0 = 8.8541878188e-12;
D = 10 * 10^-3;
w_mesh1 = 0.1;
w_meshn = 8.5;
step = 0.1;
n = (w_meshn - w_mesh1) / step + 1;
w_mesh = 10^-3 * linspace(w_mesh1, w_meshn, n);
L = -(D * mu0) / (2 * pi) * log(sin((pi * w_mesh) / (2 * D)));
 
s11_res = 1e9 *[29.95	29.95	29.9	29.85	29.85	29.75	29.65	29.65	29.6	29.45	29.35	29.2	29.1	29	28.9	28.75	28.65	28.55	28.5	28.35	28.25	28.15	28.05	27.9	27.8	27.7	27.65	27.45	27.4	27.4	27.3	27.2	27.15	27.15	27.15	27.15	27.15	27.15	27.2	27.25	27.3	27.35	27.4	27.5	27.6	27.7	27.8	27.95	28.05	28.2	28.35	28.5	28.75	29	29.2	29.45	29.7	29.9	30.2	30.45	30.8	31.15	31.5	31.85	32.2	32.65	33	33.35	33.7	34.05	34.5	34.9	35.25	35.7	36.1	36.45	36.75	37.05	37.2	37.35	37.85	38.05	38.25	38.4	38.5];
s11_angular = 2 * pi * s11_res;

C = 1 ./ (L .* s11_angular.^2);

C_restricted = C(1:70);
w_restricted = w_mesh(1:70);
alpha = 2 * D * e0 / pi;

f = @(x) -alpha * log(sin(pi * (D - x) / (2 * D))) * 10^15;
g = @(params, x) params(1)* f(params(2) * x + params(3)) + params(4);
C_data = C_restricted * 10^15;
initial_guess = [1, 1, 0, 0];
[params, ~] = lsqcurvefit(g, initial_guess, w_restricted, C_data, [], []);

p1Start = [3.7e-3, 27.15e9];
p1End = [0.4e-3, 24.65e9];
i = 53;
p2Start = [w_mesh(i), s11_res(i)];
p2End = [5e-3, 40e9];

deltaStart = p2Start - p1Start;
deltaEnd = p2End - p1End;

l_sym = @(x) -(D * mu0) / (2 * pi) * log(sin((pi * x) / (2 * D)));
s11_sym = @(x) 1./sqrt( g(params, x) * 10^-15 .* l_sym(x)) / (2 * pi);
scaleX = deltaEnd(1)/deltaStart(1);
scaleY = deltaEnd(2)/deltaStart(2);
shiftX = p1End(1) - scaleX * p1Start(1);
shiftY = p1End(2) - scaleY * p1Start(2);
startSym = s11_sym;
endSym = @(x) scaleY * s11_sym((x - shiftX) / scaleX) + shiftY;

end