function u_b = q1_e3_spectral_solver(L, T, a, sigma, N_x, N_t)

delta_t = T / N_t;

x = linspace(0, L, N_x);

u_b = zeros(N_t, N_x + 1);

u_0 = @(x, sigma, L) (1 / sqrt(2 * pi * sigma^2)) * exp(-1 * ...
    (x - (L / 2)).^2 / (2 * sigma^2));

u = u_0(x, sigma, L);

% Fix k indicies 
k = fftshift((1:N_x) - (N_x / 2) - 1);
k_vector = k * -2j * pi / L;
U_shift = exp(a * delta_t * k_vector);

for n = 2:N_t
    U = fft(u);
    U = U_shift .* U;
    u = ifft(U, 'symmetric');
    u_b(n, 1:N_x) = real(u);
    u_b(n, N_x + 1) = real(u(1));
end
end
