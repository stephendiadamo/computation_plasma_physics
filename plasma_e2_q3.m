phi_f = @(t) sin(2 * t);
rho_f = @(t) 4 * sin(2 .* t);

a = 0;
b = 2 * pi;
N = 8;
L = b - a;
h = L / N;
x = linspace(a, b, N + 1)';

phi_exact = phi_f(x);
rho = rho_f(x(1:N));

dA = diag(2 * ones(1, N));
dAp_1 = diag(-1 * ones(1, N - 1), 1);
dAm_1 = diag(-1 * ones(1, N - 1), -1);
A = (dA + dAp_1 + dAm_1);
A(1, N) = -1;
A(N, 1) = -1;
A = A / h^2;

A_1 = A;
A_1(1, :) = zeros(1, N);
A_1(1, 1) = 1;

rho_1 = rho;
rho_1(1) = 0;

A_2 = zeros(N + 1, N + 1);
A_2(1:N, 1:N) = A;
A_2(N + 1, 1:N) = h;
A_2(1:N, N + 1) = 1;

rho_2 = zeros(N + 1,1);
rho_2(1:N) = rho_f(x(1:N));

phi_1 = A_1 \ rho_1;
phi_2 = A_2 \ rho_2;

% solver_err_1 = max(abs(A_1 * phi_1 - rho_1));
% solver_err_2 = max(abs(A_2 * phi_2 - rho_2));
% display(solver_err_1);
% display(solver_err_2);

phi_num_1 = [phi_1; phi_1(1)];
phi_num_2 = [phi_2(1:N); phi_2(1)];

% plot(x, phi_exact, 'b', x, phi_num_1, 'ro-', x, phi_num_2, 'g*--');
% xlim([a b]);
% xlabel('x');
% ylabel('\phi(x)');
% legend('exact', 'numerical 1', 'numerical 2', 'location', 'northwest');

% Diagnostics 
% L1_a = h * sum(abs(phi_exact - phi_num_1));
% L2_a = sqrt(h * sum(abs(phi_exact - phi_num_1).^2));
% Linf_a = max(abs(phi_exact - phi_num_1));
% 
% L1_b = h * sum(abs(phi_exact - phi_num_2));
% L2_b = sqrt(h * sum(abs(phi_exact - phi_num_2).^2));
% Linf_b = max(abs(phi_exact - phi_num_2));
% 
% save(['q3_a_err' num2str(N)], 'h', 'L1_a', 'L2_a', 'Linf_a');
% save(['q3_b_err' num2str(N)], 'h', 'L1_b', 'L2_b', 'Linf_b');

NVec = [8 16 32 64 128 256];
hVec = zeros(6, 1);

L1Vec_a = zeros(6, 1);
L2Vec_a = zeros(6, 1);
LInfVec_a = zeros(6, 1);

L1Vec_b = zeros(6, 1);
L2Vec_b = zeros(6, 1);
LInfVec_b = zeros(6, 1);

for i = 1:6
    err_a = load(['q3_a_err' num2str(NVec(i)) '.mat']);
    err_b = load(['q3_b_err' num2str(NVec(i)) '.mat']);
    
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
title('Possibility 1');

figure 
loglog(hVec, L1Vec_b, 'o-b', hVec, L2Vec_b, 's-r', hVec, LInfVec_b, '+-g', ...
    hVec(1:3), 1e-1 * hVec(1:3).^2, '-k');
legend('L^1 error', 'L^2 error', 'L^\infty error', 'f(h) = c * h^2', ...
    'location', 'northwest');
xlabel('h');
ylabel('Errors'); 
title('Possibility 2');
