scrsz = get(0,'ScreenSize');

% White noise stuff

eps = 0.001;
v_max = 10;
v_min = -10;

f_0 = @(x, v) (1 / sqrt(2 * pi)) * exp(-v^2 / 2) * (1 + eps * randn(1));

L_x = 2 * pi / 0.015;
L_t = 600;
L_v = (v_max - v_min);

N_x = 64;
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

% Initialize f
for i = 1:N_x
    for j = 1:N_v
        f(i, j) = f_0(x(i), v(j));
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

Et = zeros(N_x, N_t);
Et(:, 1) = E;

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
    Et(:, c) = E;
    
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

N_t_FT = 2^(nextpow2(N_t) - 1);
t_FT = (0:N_t_FT - 1) * dt;
T_FT = t_FT(end);
om_ind = (1:N_t_FT) - (N_t_FT / 2) - 1; %omega mode numbers
om = 2 * pi / T_FT * om_ind;

[xx2, tt] = ndgrid(x, t_FT);
[kkx, tt2] = ndgrid(kx, t_FT);
[kkx2, oo] = ndgrid(kx, om);

Et_FT = Et(:, 1:N_t_FT);

Et_hat = fftshift( fft(Et_FT, [], 1), 1);
Et_til = fftshift( fft(Et_hat, [], 2), 2);

figure('Position',[100 100 scrsz(3)*0.9 scrsz(4)*0.6])
subplot(1,2,1)
pcolor(xx2,tt,Et_FT)
set(gca,'fontsize',16)
xlabel('x')
ylabel('t')
title('E(x,t)')
colorbar
grid off
shading interp
axis tight

subplot(1,2,2)
pcolor(kkx,tt2,log10(abs(Et_hat)))
set(gca,'fontsize',16)
xlabel('k_x')
ylabel('t')
title('log10|E(k_x,t)|')
colorbar
grid off
shading interp
axis tight

figure('Position',[100 100 scrsz(3)*0.6 scrsz(4)*0.6])
pcolor(kkx2,oo,log10(abs(Et_til)))
set(gca,'fontsize',16)
xlabel('k_x')
ylabel('\omega')
title('log10(|E(k_x,\omega)|)')
colorbar
grid off
shading interp
ylim([0 2])
hold on
plot(kx,1+3/2*kx.^2,'k:','linewidth',0.5)
legend('white noise simulation','\omega = \omega_p[1 + 3/2 (\lambda_Dk_x)^2]')
hold off
