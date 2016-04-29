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

solver_err_1 = max(abs(A_1 * phi_1 - rho_1));
solver_err_2 = max(abs(A_2 * phi_2 - rho_2));

phi_num_1 = [phi_1; phi_1(1)];
phi_num_2 = [phi_2(1:N); phi_2(1)];

plot(x,phi_exact,'b',x,phi_num_1,'ro-',x,phi_num_1,'g*--');
xlim([a b]);
xlabel('x');
ylabel('\phi(x)');
legend('exact','numerical 1','numerical 2','location','northwest');
