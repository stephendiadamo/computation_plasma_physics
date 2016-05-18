L = 1;
T = 1.5;
a = 2; 
sigma = 0.1;
N_x = 100;
N_t = 10000;

x = linspace(0, L, N_x);
t = linspace(0, T, N_t);

u_0 = @(x) (1 / sqrt(2 * pi * sigma^2)) * exp(-1 * ...
    (x - (L / 2))^2 / (2 * sigma^2));

u_exact = zeros(N_t, N_x);
for n = 1:N_t
    for j = 1:N_x
        u_exact(n, j) = u_0(x(j) - a * t(n));
    end
end

u1 = q1_e3_a_solver(L, T, a, sigma, N_x, N_t);
u2 = q1_e3_b_solver(L, T, a, sigma, N_x, N_t);

offset_1 = 10;
offset_2 = 10;
t_1 = 1;
t_2 = 1000;
num_frames = 30;

x = x(offset_1:N_x - offset_2);
u_exact = u_exact(1:end, offset_1:N_x - offset_2);
u1 = u1(1:end, offset_1:N_x - offset_2);
u2 = u2(1:end, offset_1:N_x - offset_2);

% make_movie(x, u_exact, u1, u2, num_frames);

L2_exact = zeros(N_t, 1);
L2_euler = zeros(N_t, 1);
L2_wendroff = zeros(N_t, 1);

u_exact_max = zeros(N_t, 1);
u1_max = zeros(N_t, 1);
u2_max = zeros(N_t, 1);

for n = 1:N_t
    L2_exact(n) = sqrt(sum(u_exact(n, 1:end).^2));
    L2_euler(n) = sqrt(sum(u1(n, 1:end).^2));
    L2_wendroff(n) = sqrt(sum(u2(n, 1:end).^2));
    
    u_exact_max(n) = max(u_exact(n, 1:end));
    u1_max(n) = max(u1(n, 1:end));
    u2_max(n) = max(u2(n, 1:end));
end

display('--- Conservation of L^2 norm ---');
str_1 = ['u_exact: ', num2str(max(L2_exact) - min(L2_exact))];
str_2 = ['u_euler: ', num2str(max(L2_euler) - min(L2_euler))];
str_3 = ['u_wendroff: ', num2str(max(L2_wendroff) - min(L2_wendroff))];
disp(str_1);
disp(str_2);
disp(str_3);

display('--- Conservation of maximums ---');
str_1 = ['u_exact: ', num2str(max(u_exact_max) - min(u_exact_max))];
str_2 = ['u_euler: ', num2str(max(u1_max) - min(u1_max))];
str_3 = sprintf('u_wendroff: %f', max(u2_max) - min(u2_max));
disp(str_1);
disp(str_2);
disp(str_3);

% Conservation of L^2 
% - u_exact has conservation and then loses it
% - u1 has fairly good convervation
% - u2 also has fairly good conservation
% - results improve as N_x increases
plot(1:N_t, L2_exact, 1:N_t, L2_euler, 1:N_t, L2_wendroff);
title('L^2 norm conservation');
legend('Exact', 'Euler', 'Wendroff');
figure;

plot(1:N_t, u_exact_max, 1:N_t, u1_max, 1:N_t, u2_max);
title('Maximum of u conservation');
legend('Exact', 'Euler', 'Wendroff');
