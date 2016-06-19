eps = 0.01;
T = 20;
k = 0.5;
Lx = 2 * pi / k;
Lv = 10;
d = 3;

% Initial condition of f at time 0
f_0 = @(x, v) (1 + eps * cos(k * x)) * (1 / sqrt(2 * pi)) * exp(-v^2 / 2);
% Importance sampling distribution 
g_0 = @(x, v) (1 / Lx) * (1 / sqrt(2 * pi)) * exp(-v^2 / 2); 

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

kx_ind = (1:Nx) - Nx / 2 - 1;
kx = 2 * pi / Lx * kx_ind;
[kkx, vv] = ndgrid(kx, v);
X_shift = exp(-1j * vv .* dt .* kkx);

% Initialize f at time 0
f = zeros(Nx, Nv);
for i = 1:Nx
    for j = 1:Nv
        f(i, j) = f_0(x(i), v(j));
    end
end

knot_vector = zeros(d + 2 , 1);
for i = 1:d + 1
    knot_vector(i) = (-(d + 1) / 2 + (i - 1)) * dx;
end
knot_vector(d + 2) = dx * (d + 1) / 2;
% Part (b) Plotting the bsplines and checking their normalizations
% disp(knot_vector);
bspline(knot_vector);

Sd = bspline(knot_vector);

% for n = 1:Nt
%     % Part (a) Solve Poisson equation via spectral method to determine 
%     % the electric field at the grid points
%     density = dv * sum(f, 2);
%     rho = fft(1 - density);
%     E = fftshift(rho) ./ (1j * kx');
%     E(Nx / 2 + 1) = 0; % setting median to zero
%     E = ifft(fftshift(E), 'symmetric');
% end

