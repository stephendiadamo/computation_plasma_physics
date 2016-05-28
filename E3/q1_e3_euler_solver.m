% Euler solver
function u = q1_e3_euler_solver(L, T, a, sigma, N_x, N_t)

u = zeros(N_t, N_x);
h = L / N_x;
delta_t = T / N_t;

x = linspace(0, L, N_x);
a_minus = min(a, 0);
a_plus = max(a, 0);

u_0 = @(x, sigma, L) (1 / sqrt(2 * pi * sigma^2)) * exp(-1 * ...
    (x - (L / 2))^2 / (2 * sigma^2));

euler = @(u_n_j, u_n_j_p1, u_n_j_m1) u_n_j - (delta_t / h) * (a_minus * ...
    (u_n_j_p1 - u_n_j) + a_plus * (u_n_j - u_n_j_m1));

for i = 1:N_x
    u(1, i) = u_0(x(i), sigma, L);
end

for n = 1:N_t - 1
    for j = 1:N_x
        if (j == 1)
            u(n + 1, j) = euler(u(n, j), u(n, j + 1), u(n, N_x));
        elseif (j == N_x)
            u(n + 1, j) = euler(u(n, j), u(n, 1), u(n, j - 1));
        else
            u(n + 1, j) = euler(u(n, j), u(n, j + 1), u(n, j - 1));
        end
    end
end

end
