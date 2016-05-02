phi_f = @(t) t .* cos(t);
dphi_f = @(t) cos(t) - t .* sin(t);
rho_f = @(t) 2 .* sin(t) + t .* cos(t);

a = 0;
b = 2 * pi;
N = 256;
L = b - a;
h = L / N;
x = linspace(a, b, N + 1)';
phi_exact = phi_f(x);

alpha = phi_f(a);
gamma = dphi_f(b);

% Part A
clear rho A phi;
rho = rho_f(x(2:N + 1)); 
rho(1) = rho(1) + alpha / h^2;
rho(N) = rho(N) + 2 * gamma / h;
phi_a = q2_finite_difference_solver_a(a, b, alpha, N, rho);
gamma_approx_a = (phi_a(N + 1) - phi_a(N - 1)) / (2 * h);
display(gamma_approx_a);

% Part B
clear rho A phi;
rho = rho_f(x(2:N + 1)); 
rho(1) = rho(1) + alpha / h^2;
rho(N) = rho(N) + gamma / h;
phi_b = q2_finite_difference_solver_b(a, b, alpha, N, rho);
gamma_approx_b = (phi_b(N + 1) - phi_b(N)) / h;
display(gamma_approx_b);

plot(x, phi_exact, 'b', x, phi_a, 'ro-', x ,phi_b, 'g*-');
xlim([a b]);
xlabel('x');
ylabel('\phi(x)');
legend('exact', 'numerical a)', 'numerical b)', 'location', 'northwest');

% Diagnostics 
% L1_a = h * sum(abs(phi_exact - phi_a));
% L2_a = sqrt(h * sum(abs(phi_exact - phi_a).^2));
% Linf_a = max(abs(phi_exact - phi_a));
% 
% L1_b = h * sum(abs(phi_exact - phi_b));
% L2_b = sqrt(h * sum(abs(phi_exact - phi_b).^2));
% Linf_b = max(abs(phi_exact - phi_b));
% 
% save(['q2_a_err' num2str(N)], 'h', 'L1_a', 'L2_a', 'Linf_a');
% save(['q2_b_err' num2str(N)], 'h', 'L1_b', 'L2_b', 'Linf_b');

NVec = [8 16 32 64 128 256];
hVec = zeros(6, 1);

L1Vec_a = zeros(6, 1);
L2Vec_a = zeros(6, 1);
LInfVec_a = zeros(6, 1);

L1Vec_b = zeros(6, 1);
L2Vec_b = zeros(6, 1);
LInfVec_b = zeros(6, 1);

for i = 1:6
    err_a = load(['q2_a_err' num2str(NVec(i)) '.mat']);
    err_b = load(['q2_b_err' num2str(NVec(i)) '.mat']);
    
    hVec(i) = err_a.h;
    
    L1Vec_a(i) = err_a.L1_a;
    L2Vec_a(i) = err_a.L2_a;
    LInfVec_a(i) = err_a.Linf_a;
    
    L1Vec_b(i) = err_b.L1_b;
    L2Vec_b(i) = err_b.L2_b;
    LInfVec_b(i) = err_b.Linf_b;
end

figure
loglog(hVec, L1Vec_a, 'o-b', hVec, L2Vec_a, 's-r', hVec, LInfVec_a, '+-g', ...
    hVec(1:3), 1e-1 * hVec(1:3).^2, '-k');
legend('L^1 error', 'L^2 error', 'L^\infty error', 'f(h) = c * h^2', ...
    'location', 'northwest');
xlabel('h');
ylabel('Errors'); 
title('Case a');

figure 
loglog(hVec, L1Vec_b, 'o-b', hVec, L2Vec_b, 's-r', hVec, LInfVec_b, '+-g', ...
    hVec(1:3), hVec(1:3), '-k');
legend('L^1 error', 'L^2 error', 'L^\infty error', 'f(h) = c * h', ...
    'location', 'northwest');
xlabel('h');
ylabel('Errors'); 
title('Case b');
