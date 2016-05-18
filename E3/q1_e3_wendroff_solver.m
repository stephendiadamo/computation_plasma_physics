% Wendroff solver
function u = q1_e3_wendroff_solver(L, T, a, sigma, N_x, N_t)

u = zeros(N_t, N_x);
h = L / N_x;
delta_t = T / N_t;

x = linspace(0, L, N_x);

u_0 = @(x, sigma, L) (1 / sqrt(2 * pi * sigma^2)) * exp(-1 * ...
    (x - (L / 2))^2 / (2 * sigma^2));

lax_wendroff = @(u_n_j, u_n_j_p1, u_n_j_m1) u_n_j - ... 
    (a * delta_t / (2 * h)) * (u_n_j_p1 - u_n_j_m1) + ...
    (a^2 * delta_t^2 / (2 * h^2)) * (u_n_j_p1 - 2 * u_n_j + u_n_j_m1);

for i = 1:N_x
    u(1, i) = u_0(x(i), sigma, L);
end

for n = 1:N_t - 1
    for j = 1:N_x
        if (j == 1)
            u(n + 1, j) = lax_wendroff(u(n, j), u(n, j + 1), u(n, N_x));
        elseif (j == N_x)
            u(n + 1, j) = lax_wendroff(u(n, j), u(n, 1), u(n, j - 1));
        else
            u(n + 1, j) = lax_wendroff(u(n, j), u(n, j + 1), u(n, j - 1));
        end
    end
end

end