function phi = finite_difference_solver(a, b, alpha, beta, N, rho)
L = b - a;
h = L / N;

dA = diag(2 * ones(1, N - 1)); 
dAp1 = diag(-1 * ones(1, N - 2), 1);
dAm1 = diag(-1 * ones(1, N - 2), -1);

A = dA + dAp1 + dAm1;
A = A / h^2;

phi = A \ rho;
phi = [alpha; phi; beta];
end