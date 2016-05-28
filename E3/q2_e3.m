v_max = 5;
v_min = -5;

L_t = 30;
L_x = 12;
L_v = (v_max - v_min);

N_x = 64;
N_v = 64;
N_t = 200;

x = linspace(0, L_x, N_x + 1);
v = linspace(v_min, v_max, N_v + 1);
t = linspace(0, L_t, N_t);

dx = L_x / N_x;
dv = L_v / N_v;
dt = L_t / N_t;

f_0 = @(x, v) (1 + 0.01 * cos(2 * pi * x / L_x)) * ... 
    (1 / sqrt(2 * pi)) * exp(-1 * v^2 / 2);

f = zeros(N_t, N_x, N_v);

for i = 1:N_x
    for j = 1:N_v
        f(0, i, j) = f_0(x(i), v(j));
    end
end

dA = diag(2 * ones(1, N_x)); % diagonal matrix
dAp1 = diag(-1 * ones(1, N_x - 1), 1); % super-diagonal matrix
dAm1 = diag(-1 * ones(1, N_x - 1), -1); % sub-diagonal matrix
A = (dA + dAp1 + dAm1);
A(1 ,N_x) = -1;
A(N_x, 1) = -1;

A(1, :) = zeros(1, N_x);
A(1, 1) = 1;
A; % Don't forget to multiply the rhs by dx^2 and to set rhs(1)=0.

rhs = zeros(N_x, 1);











