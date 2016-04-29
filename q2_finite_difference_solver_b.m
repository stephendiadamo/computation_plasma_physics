function phi = q2_finite_difference_solver_b(a, b, alpha, N, rho)
L = b - a;
h = L / N;

dA = diag([ 2 * ones(1, N - 1) 1]);
dAp1 = diag(-1 * ones(1, N - 1), 1);
dAm1 = diag(-1 * ones(1, N - 1), -1);
A = (dA + dAp1 + dAm1) / h^2;

phi = A \ rho;
phi = [alpha; phi]; 
end