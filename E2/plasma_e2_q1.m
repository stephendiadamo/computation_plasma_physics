phi_f = @(t) t .* cos(t);
rho_f = @(t) 2 .* sin(t) + t .* cos(t);

a = 0;
b = 2 * pi;
N = 256;
L = b - a;
h = L / N;
x = linspace(a, b, N + 1)';

alpha = phi_f(a);
beta = phi_f(b);

rho = rho_f(x(2:N));
rho(1) = rho(1) + alpha / h^2;
rho(N - 1) = rho(N - 1) + beta / h^2;

phi = q1_finite_difference_solver(a, b, a, b, N, rho);
phi_exact = phi_f(x);

plot(x, phi_exact, 'b', x, phi, 'ro-');
xlim([a b]);
xlabel('x');
ylabel('\phi(x)');
legend('exact', 'numerical', 'location', 'northwest');

% Diagnostics 
% L1 = h * sum(abs(phi_exact - phi));
% L2 = sqrt(h * sum(abs(phi_exact - phi).^2));
% Linf = max(abs(phi_exact - phi));
% save(['q1_err' num2str(N)], 'h', 'L1', 'L2', 'Linf'); 

NVec = [8 16 32 64 128 256];
hVec = zeros(6, 1);
L1Vec = zeros(6, 1);
L2Vec = zeros(6, 1);
LInfVec = zeros(6, 1);

for i = 1:6
    err = load(['q1_err' num2str(NVec(i)) '.mat']);

    hVec(i) = err.h;
    L1Vec(i) = err.L1;
    L2Vec(i) = err.L2;
    LInfVec(i) = err.Linf;
end

loglog(hVec, L1Vec, 'o-b', hVec, L2Vec, 's-r', hVec, LInfVec, '+-g', ...
    hVec(1:3), 1e-1 * hVec(1:3).^2, '-k');
legend('L^1 error', 'L^2 error', 'L^\infty error', 'f(h) = 0.1 * h^2', ...
    'location', 'northwest');
xlabel('h');
ylabel('Errors'); 
