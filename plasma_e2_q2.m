a = 0;
b = 2 * pi;
N = 32;
L = b - a;
h = L / N;
x = linspace(a, b, N + 1)';

phi_f = @(t) t .* cos(t);
dphi_f = @(t) cos(t) - t .* sin(t);
rho_f = @(t) 2.*sin(t) + t .* cos(t);

alpha = phi_f(a);
gamma = dphi_f(b);

% Part A
rho = rho_f(x(2:N + 1)); % exact rho at the grid points
rho(1) = rho(1) + alpha / h^2; % add boundary term
rho(N) = rho(N) + 2 * gamma / h; % add boundary term

dA = diag(2 * ones(1, N));
dAp1 = diag(-1 * ones(1, N - 1), 1);
dAm1 = diag([-1 * ones(1, N - 2) -2], -1);
A = (dA + dAp1 + dAm1) / h^2;

phi = A \ rho;
phi_num_a = [alpha; phi];

% Part B
clear rho A phi

rho = rho_f(x(2:N + 1));
rho(1) = rho(1) + alpha / h^2;
rho(N) = rho(N) + gamma / h;

dA = diag([ 2 * ones(1, N-1) 1]);
dAp1 = diag(-1 * ones(1, N - 1), 1); % super-diagonal matrix
dAm1 = diag(-1 * ones(1, N - 1), -1); % sub-diagonal matrix
A = (dA + dAp1 + dAm1);
A = A/h^2;

phi = A \ rho;
solver_err = max(abs(A * phi - rho));

phi_exact = phi_f(x);
phi_num_b = [alpha; phi];

plot(x, phi_exact, 'b', x, phi_num_a, 'ro-', x ,phi_num_b, 'g*-');
xlim([a b]);
xlabel('x');
ylabel('\phi(x)');
legend('exact', 'numerical a)', 'numerical b)', 'location', 'northwest');
