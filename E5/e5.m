scrsz = get(0,'ScreenSize');
eps = 0.001;
T = 20;
k = 0.5;
Lx = 2 * pi / k;
Lv = 10;
d = 3;

k_vals = [0.5, 0.4, 0.3, 0.2];
w_r = [1.4156, 1.2850, 1.1598, 1.0640];
w_i = [-0.1533, -0.0661, -0.0126, -5.51e-5];
r = [0.3677, 0.4247, 0.6368, 1.1297];
p = [0.5362, 0.3358, 0.1143, 0.0013];
k_index = 1;

% Landau Damping:
E_analytical = @(x, t, k_ind) (2 * eps / k_vals(k_ind)) * r(k_ind) * exp(w_i(k_ind) * t) ...
                       * cos(w_r(k_ind) * t - p(k_ind)) * sin(k_vals(k_ind) * x);
                   
% Initial condition of f at time 0
f_0 = @(x, v) (1 + eps * cos(k * x)) * (1 / sqrt(2 * pi)) * exp((-v^2) / 2);
% Importance sampling distribution 
g_0 = @(v) (1 / Lx) * (1 / sqrt(2 * pi)) * exp(-v^2 / 2); 
% Maxwellian 
M = @(v) exp(-v^2 / 2) / sqrt(2 * pi);

Nv = 32;
Nx = 32;
% Nk = 5e4;
Nk = 4000;

dt = 0.1;
dx = Lx / Nx;
dv = Lv / Nv;

Nt = floor(T / dt) + 1;

% Diagnostics
momentum = zeros(Nt, 1);
P_energy = zeros(Nt, 1);
E_energy = zeros(Nt, 1);

t = linspace(0, T, Nt);
xj = linspace(0, Lx, Nx);
vj = linspace(-Lv / 2, Lv / 2, Nv);

knot_vector = zeros(d + 2 , 1);
for i = 1:d + 1
    knot_vector(i) = (-(d + 1) / 2 + (i - 1)) * dx;
end
knot_vector(d + 2) = dx * (d + 1) / 2;

% Part b
% Plotting the bsplines and checking their normalizations
Sd = bspline(knot_vector);
S = @(x) ppval(Sd, x) / dx;

% Normalization
spline_normalization = integral(S, knot_vector(1), knot_vector(end));
disp('-- Spline normaliztion --');
disp(spline_normalization);

% Initializing particles
xk = zeros(Nt, Nk);
vk = zeros(Nt, Nk);
% Part c
% Generate random positions for particles according to g_0
xk(1, :) = Lx * rand(Nk, 1);
vk(1, :) = randn(Nk, 1);
wk = zeros(Nk, 1);

% Compute the weights 
for marker = 1:Nk 
    wk(marker) = f_0(xk(1, marker), vk(1, marker)) / g_0(vk(1, marker));
end

% Part d
% The right side of equation (3) is the expectation value for 
% the position random variable Sd(x - X) * f(t, X, V)/g(t, X, V). 
% The estimator for X is 
% 1 / Nk * sum( Sd(x - xk) * f(t, xk, vk) / g(t, xk, vk) ).

% Part e 
% Refer to notes for derivation seen in Tutorial
nj = zeros(Nt, Nx);
nj(1, :) = get_density(1, xk, vk, xj, Nx, Nk, Lx, S, wk, '', 0);

% Part f
% Calculating the electric field
Ej = zeros(Nt, Nx);
kx_ind = (1:Nx) - Nx / 2 - 1;
kx = 2 * pi / Lx * kx_ind;

% Part j
% Analyical solution for electric field 
E_ana = zeros(Nx, Nt);                   
for i = 1:Nx
    for j = 1:Nt
        E_ana(i, j) = E_analytical(xj(i), t(j), k_index);
    end
end
% Analytical electric field energy w.r.t. time
E_ana_energy = dx * sum(E_ana .* E_ana) / 2;

density = 1 - nj(1, :)';
rho = fft(density);
E = fftshift(rho) ./ (1j * kx');
E(Nx / 2 + 1) = 0; % setting median to zero
Ej(1, :) = ifft(fftshift(E), 'symmetric');

E_energy(1) = dx * sum(E .* E) / 2;

Eh = zeros(Nt, Nk);
for marker = 1:Nk 
    s = 0;
    for j = 1:Nx
        if (abs(xj(j) - xk(1, marker)) <= knot_vector(end))
            s = s + Ej(1, j) * S(xk(1, marker) - xj(j));
        end
    end
    Eh(1, marker) = dx * s;
end

vk_half = zeros(Nt, Nk);
for time = 1:Nt-1    
    for marker = 1:Nk
        vk_half(time, marker) = vk(time, marker) + dt * Eh(time, marker);
        
        % Advance xk
        xk(time + 1, marker) = xk(time, marker) + dt * vk_half(time, marker);
    end
    
    % Update density
    nj(time + 1, :) = get_density(time, xk, vk, xj, Nx, Nk, Lx, S, wk, 'verlet', vk_half);
        
    % Update Ej(t + 1, x)
    density = 1 - nj(time + 1, :)';
    rho = fft(density);
    E = fftshift(rho) ./ (1j * kx');
    E(Nx / 2 + 1) = 0; % setting median to zero
    Ej(time + 1, :) = ifft(fftshift(E), 'symmetric');
    
    E_energy(time + 1) = dx * sum(E .* E) / 2;
    
    % Update Eh(time + 1, marker)    
    % TODO: This doesn't handle the boundaries
    for marker = 1:Nk 
        s = 0;
        for j = 1:Nx
            if (abs(xj(j) - xk(time + 1, marker)) <= knot_vector(end))
                s = s + Ej(time + 1, j) * S(xk(time + 1, marker) - xj(j));
            end
        end
        Eh(time + 1, marker) = dx * s;
        % Advance vk              
        vk(time + 1, marker) = vk_half(time, marker) + dt * Eh(time + 1, marker);
    end  
end

f = zeros(Nt, Nx, Nv);
for i = 1:Nx
    for j = 1:Nv
        f(1, i, j) = f_0(xj(i), vj(j));
    end
end

x_b = linspace(0,Lx,Nx+1);
v_b = linspace(-Lv/2,Lv/2,Nv+1);
[xx_b,vv_b] = ndgrid(x_b,v_b);
f_b = zeros(Nx+1,Nv+1);
f0_b = zeros(Nx+1, Nv+1);

for i = 1:Nx+1
    for j = 1:Nv+1
        f0_b(i, j) = f_0(x_b(i), v_b(j));
    end
end

for time = 1:Nt
    [N, Xedges, Yedges, Xbin, Ybin] = histcounts2(xk(time, :), vk(time, :), xj, vj);
    coord_bin = zeros(Nk, 2);
    
    % Clearly a bug if I have to do this!
    Xbin = Xbin + 1;
    Ybin = Ybin + 1;
    
    for marker = 1:Nk
        coord_bin(marker, 1) = Xbin(marker);
        coord_bin(marker, 2) = Ybin(marker);
    end

    f = (1 / Nk) * accumarray(coord_bin, wk, [32, 32]);

    f(1, 1) = 0;
    f(32, 32) = 0;
    
    f_b(1:Nx,1:Nv) = f;
    f_b(Nx+1,:) = f_b(1,:);
    f_b(:,Nv+1) = f_b(:,1);
    
    u0 = 0;
    w0 = 0;
    for i = 1:Nx
        u0 = u0 + dx*dv*sum( vj .* f(i,:) ); % initial momentum
        w0 = w0 + 1/2*dx*dv*sum((vj.^2) .* f(i,:) ); % initial plasma energy
    end

    momentum(time) = u0;
    P_energy(time) = w0;
    
%     if (mod(time - 1, 40) == 0)
%         figure        
%         surf(xx_b,vv_b,f_b - f0_b)        
%         title(['t = ' num2str(t(time),'%5.3f')])
%         grid off
%         shading interp
%         colorbar
%         view([0 90])
%         xlim([0 Lx])
%         ylim([-Lv/2 Lv/2])
%         xlabel('x')
%         ylabel('v')
%         drawnow
%    end   
end

figure('Position',[800 100 scrsz(3) scrsz(4)])

subplot(2, 2, 1)
plot(t, momentum)
set(gca,'fontsize',13)
title('momentum')
xlabel('t')

subplot(2, 2, 2)
plot(t, P_energy)
set(gca,'fontsize',13)
title('plasma energy')
xlabel('t')

subplot(2, 2, 3)
semilogy(t, abs(E_energy))
set(gca,'fontsize',13)
title('field energy')
xlabel('t')

subplot(2, 2, 4)
plot(t,abs(P_energy) + abs(E_energy))
set(gca,'fontsize',13)
title('total energy')
xlabel('t')

figure('Position',[100 100 scrsz(3)*0.6 scrsz(4)*0.6])
semilogy(t, abs(E_energy), t, E_ana_energy);
set(gca,'fontsize', 13);
legend('numerics','analytic');
title('field energy');
xlabel('t');

display('done');
