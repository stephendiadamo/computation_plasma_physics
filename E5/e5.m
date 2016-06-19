eps = 0.01;
T = 20;
k = 0.5;
Lx = 2 * pi / k;
Lv = 10;
d = 3;

% Initial condition of f at time 0
f_0 = @(x, v) (1 + eps * cos(k * x)) * (1 / sqrt(2 * pi)) * exp(-v^2 / 2);

Nv = 32;
Nx = 32;
Nk = 5e4;

dt = 1e-1;
dx = Lx / Nx;
dv = Lv / Nv;

Nt = floor(T / dt) + 1;

t = linspace(0, T, Nt);
x = linspace(0, Lx, Nx);
v = linspace(-Lv / 2, Lv / 2, Nv);

% Initialize f at time 0
f = zeros(Nx, Nv);
for i = 1:Nx
    for j = 1:Nv
        f(i, j) = f_0(x(i), v(j));
    end
end








