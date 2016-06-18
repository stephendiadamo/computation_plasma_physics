scrsz = get(0,'ScreenSize');

% White noise stuff

eps = 0.001;
v_max = 10;
v_min = -10;

f_0 = @(x, v, v0) 0.5 * ((1 / sqrt(2 * pi)) * ( exp(-(v + v0)^2 / 2) + ...
    exp(-(v - v0)^2 / 2))) * (1 + eps * cos(0.2 * x));

L_x = 2 * pi / 0.2;
L_t = 50;
L_v = (v_max - v_min);

N_x = 256;
N_v = 256;

dt = 1e-1; % no CFL condition for spectral-splitting method
dx = L_x / N_x;
dv = L_v / N_v;

N_t = floor(L_t / dt) + 1;

x = linspace(0, L_x, N_x);
v = linspace(v_min, v_max, N_v);
t = linspace(0, L_t, N_t);

% Diagnostic arrays
E_energy = zeros(N_t, 1);

f = zeros(N_x, N_v);
phi_periodic = zeros(N_x + 1, 1);

% v0 = 1.3;
% v0 = 2.4;
v0 = 3.0;

% Initialize f
for i = 1:N_x
    for j = 1:N_v
        f(i, j) = f_0(x(i), v(j), v0);
    end
end

density = dv * sum(f, 2);

dA = diag(2 * ones(1, N_x)); % diagonal matrix
dAp1 = diag(-1 * ones(1, N_x - 1), 1); % super-diagonal matrix
dAm1 = diag(-1 * ones(1, N_x - 1), -1); % sub-diagonal matrix
A = (dA + dAp1 + dAm1);
A(1 ,N_x) = -1;
A(N_x, 1) = -1;

A(1, :) = zeros(1, N_x);
A(1, 1) = 1;
A; 

rhs = (1 - density) * dx^2;
rhs(1) = 0;
phi = A \ rhs;

phi_periodic(1:N_x) = phi;
phi_periodic(N_x + 1) = phi(1);
phi_periodic_plus = circshift(phi_periodic, [-1 0]);
phi_periodic_minus = circshift(phi_periodic, [1 0]);
E_periodic = -(phi_periodic_plus - phi_periodic_minus) / 2 / dx;
E = E_periodic(1:N_x);

E_energy(1) = dx * sum(E .* E) / 2;

kx_ind = (1:N_x) - N_x / 2 - 1;
kx = 2 * pi / L_x * kx_ind;
[kkx, vv] = ndgrid(kx, v);
X_shift = exp(-1j * vv .* dt .* kkx);

% method = 0 for part a (E via finite difference)
% methid = 1 for part d (E via spectral) 
method = 1;

for c = 2:N_t
    if method == 0
        rhs = (1 - density) * dx^2;
        rhs(1) = 0;
        phi = A \ rhs;
   
        phi_plus = circshift(phi, [-1 0]);
        phi_minus = circshift(phi, [1 0]);
    
        E = -(phi_plus - phi_minus) / 2 / dx;
    
        E_periodic(1:N_x) = E; 
        E_periodic(N_x + 1) = E(1);

        phi_periodic(1:N_x) = phi;
        phi_periodic(N_x + 1) = phi(1);
    elseif method == 1
        rho = fft(1 - density);
        E = fftshift(rho) ./ (1j * kx');
        E(N_x / 2 + 1) = 0; % setting mean to zero

        E = ifft(fftshift(E), 'symmetric');
        E_periodic(1:N_x) = E;
        E_periodic(N_x + 1) = E(1);
    else
        break;
    end

    E_energy(c) = dx * sum(E .* E) / 2;

    kv_ind = (1:N_v) - N_v / 2 - 1;
    kv = 2 * pi / L_v * kv_ind;
    [EE, kkv] = ndgrid(E, kv);
    V_shift = exp(1j * EE .* dt .* kkv);
    
    kx_ind = (1:N_x) - N_x / 2 - 1;
    kx = 2 * pi / L_x * kx_ind;
    [kkx, vv] = ndgrid(kx, v);
    X_shift = exp(-1j * vv .* dt .* kkx);
    
    Fv = fft(f, [], 2);
    Fv = V_shift .* fftshift(Fv, 2);
    f = ifft(fftshift(Fv, 2), [], 2, 'symmetric');
    
    Fx = fft(f, [], 1);
    Fx = X_shift .* fftshift(Fx, 1);
    f = ifft(fftshift(Fx, 1), [], 1, 'symmetric');
    
    density = dv * sum(f, 2);
end

% wi = 0.0011;
% wi = 0.2258;
wi = 0.2845;

E_ana = 2e-6 * exp(2 * wi * t);

figure('Position',[100 100 scrsz(3)*0.6 scrsz(4)*0.6])
semilogy(t, E_energy, t, E_ana, ':','linewidth',3)
set(gca,'fontsize',16)
xlabel('t')
legend('numerics','analytic growth rate','location','best')
