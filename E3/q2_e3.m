scrsz = get(0,'ScreenSize');

v_max = 5;
v_min = -5;

L_t = 200;
L_x = 12;
L_v = (v_max - v_min);

N_x = 128;
N_v = 128;

dt = 1e-1; % no CFL condition for spectral-splitting method
dx = L_x / N_x;
dv = L_v / N_v;

N_t = floor(L_t / dt) + 1;

x_b = linspace(0, L_x, N_x + 1);
v_b = linspace(v_min, v_max, N_v + 1);
[xx_b,vv_b] = ndgrid(x_b, v_b);

x = x_b(1:N_x);
v = v_b(1:N_v);
t = linspace(0, L_t, N_t);

% Diagnostic arrays
momentum = zeros(N_t, 1);
P_energy = zeros(N_t, 1);
E_energy = zeros(N_t, 1);
max_f = zeros(N_t, 1);
mass_f = zeros(N_t, 1);
L2_f = zeros(N_t, 1);

f_0 = @(x, v) (1 + 0.01 * cos(2 * pi * x / L_x)) * ... 
    (1 / sqrt(2 * pi)) * exp(-1 * v^2 / 2);

f = zeros(N_x, N_v);
phi_periodic = zeros(N_x + 1, 1);

% Initialize f
for i = 1:N_x
    for j = 1:N_v
        f(i, j) = f_0(x(i), v(j));
    end
end

max_f(1) = max(max(f));
mass_f(1) = dx * dv * sum(sum(f));
L2_f(1) = sqrt(dx * dv * sum(sum(f.^2)));

u0 = 0;
w0 = 0;
for i = 1:N_x
    u0 = u0 + dx * dv * sum(v .* f(i,:)); % initial momentum
    w0 = w0 + 1/2 * dx * dv * sum((v.^2) .* f(i,:)); % initial plasma energy
end

momentum(1) = u0;
P_energy(1) = w0;

f_periodic = zeros(N_x + 1, N_v + 1);
f_periodic(1:N_x, 1:N_v) = f(:, :);
f_periodic(N_x + 1, 1:N_v) = f(1, :);
f_periodic(1:N_x, N_v + 1) = f(:, 1);

% figure
% surf(xx_b,vv_b,f_periodic)
% set(gca,'fontsize',16)
% axis tight
% grid off
% shading interp
% colorbar
% %view([0 90])
% xlabel('x')
% ylabel('v')
% title('Initial condition')
%pause

density = dv * sum(f, 2);

dA = diag(2 * ones(1, N_x)); % diagonal matrix
dAp1 = diag(-1 * ones(1, N_x - 1), 1); % super-diagonal matrix
dAm1 = diag(-1 * ones(1, N_x - 1), -1); % sub-diagonal matrix
A = (dA + dAp1 + dAm1);
A(1 ,N_x) = -1;
A(N_x, 1) = -1;

A(1, :) = zeros(1, N_x);
A(1, 1) = 1;
A; % Don't forget to multiply the rhs by dx^2 and to set rhs(1)=0.

rhs = (1 - density) * dx^2;
rhs(1) = 0;
phi = A \ rhs;

% Solving E from phi (E = d_phi / dx)
phi_periodic(1:N_x) = phi;
phi_periodic(N_x + 1) = phi(1);
phi_periodic_plus = circshift(phi_periodic, [-1 0]);
phi_periodic_minus = circshift(phi_periodic, [1 0]);
E_periodic = -(phi_periodic_plus - phi_periodic_minus) / 2 / dx;
E = E_periodic(1:N_x);

% Conservation of electric field energy
E_energy(1) = dx * sum(E .* E) / 2;

% method = 0 for part a (splitting), and 1 for part d (spectral) 
method = 1;

% creating shifts for FFT
kv_ind = (1:N_v) - N_v / 2 - 1;
kv = 2 * pi / L_v * kv_ind;
[EE, kkv] = ndgrid(E, kv);
V_shift = exp(1j * EE .* dt .* kkv);

kx_ind = (1:N_x) - N_x / 2 - 1;
kx = 2 * pi / L_x * kx_ind;
[kkx, vv] = ndgrid(kx, v);
X_shift = exp(-1j * vv .* dt .* kkx);

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
        E = fftshift(rho)./(1j * kx');
        E(N_x / 2 + 1) = 0; % setting median to zero
        E = ifft(fftshift(E), 'symmetric');
        E_periodic(1:N_x) = E;
        E_periodic(N_x + 1) = E(1);
    else
        break;
    end
    
    E_energy(c) = dx * sum(E .* E) / 2;
    
    Fv = fft(f, [], 2);
    Fv = V_shift .* fftshift(Fv, 2);
    f = ifft(fftshift(Fv, 2), [], 2, 'symmetric');
    
    Fx = fft(f, [], 1);
    Fx = X_shift .* fftshift(Fx, 1);
    f = ifft(fftshift(Fx, 1), [], 1, 'symmetric');
    
    density = dv * sum(f, 2);
    
    max_f(c) = max(max(f));
    mass_f(c) = dx * dv * sum(sum(f));
    L2_f(c) = sqrt(dx * dv * sum(sum(f.^2)));

    u0 = 0;
    w0 = 0;
    for i = 1:N_x
        u0 = u0 + dx*dv*sum( v.*f(i,:) ); % initial momentum
        w0 = w0 + 1/2*dx*dv*sum( (v.^2).*f(i,:) ); % initial plasma energy
    end

    momentum(c) = u0;
    P_energy(c) = w0;
end

figure
plot(t,max_f,t,mass_f,t,L2_f,'linewidth',2)
set(gca,'fontsize',13)
legend('maxf','massf','L2f')
xlabel('t')
title('Conservation properties')

figure('Position',[800 100 scrsz(3) scrsz(4)])

subplot(2, 2, 1)
plot(t,momentum)
set(gca,'fontsize',13)
title('momentum')
xlabel('t')

subplot(2, 2, 2)
plot(t,P_energy)
set(gca,'fontsize',13)
title('plasma energy')
xlabel('t')

subplot(2, 2, 3)
semilogy(t, E_energy)
set(gca,'fontsize',13)
title('field energy')
xlabel('t')

subplot(2, 2, 4)
plot(t,P_energy + E_energy)
set(gca,'fontsize',13)
title('total energy')
xlabel('t')
