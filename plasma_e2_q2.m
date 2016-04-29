phi_f = @(t) t .* cos(t);
dphi_f = @(t) cos(t) - t .* sin(t);
rho_f = @(t) 2 .* sin(t) + t .* cos(t);

a = 0;
b = 2 * pi;
N = 8;
L = b - a;
h = L / N;
x = linspace(a, b, N + 1)';
phi_exact = phi_f(x);

alpha = phi_f(a);
gamma = dphi_f(b);
rho = rho_f(x(2:N + 1)); 
rho(1) = rho(1) + alpha / h^2;

% Part A
rho(N) = rho(N) + 2 * gamma / h;
phi_a = q2_finite_difference_solver_a(a, b, alpha, N, rho);

% Part B
rho(N) = rho(N) + gamma / h;
phi_b = q2_finite_difference_solver_b(a, b, alpha, N, rho);

plot(x, phi_exact, 'b', x, phi_a, 'ro-', x ,phi_b, 'g*-');
xlim([a b]);
xlabel('x');
ylabel('\phi(x)');
legend('exact', 'numerical a)', 'numerical b)', 'location', 'northwest');




