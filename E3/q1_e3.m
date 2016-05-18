L = 1;
T = 1.5;
a = 2; 
sigma = 0.1;
N_x = 1000;
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

u_euler = q1_e3_euler_solver(L, T, a, sigma, N_x, N_t);
u_wendroff = q1_e3_wendroff_solver(L, T, a, sigma, N_x, N_t);

offset_1 = 10;
offset_2 = 10;
t_1 = 1;
t_2 = 1000;
num_frames = 30;

x = x(offset_1:N_x - offset_2);
u_exact = u_exact(1:end, offset_1:N_x - offset_2);
u_euler = u_euler(1:end, offset_1:N_x - offset_2);
u_wendroff = u_wendroff(1:end, offset_1:N_x - offset_2);

% make_movie(x, u_exact, u_euler, u_wendroff, num_frames);

mass_exact = zeros(N_t, 1);
mass_euler = zeros(N_t, 1);
mass_wendroff = zeros(N_t, 1);

L2_exact = zeros(N_t, 1);
L2_euler = zeros(N_t, 1);
L2_wendroff = zeros(N_t, 1);

u_exact_max = zeros(N_t, 1);
u_euler_max = zeros(N_t, 1);
u_wandroff_max = zeros(N_t, 1);

for n = 1:N_t
    mass_exact(n) = sum(u_exact(n, 1:end));
    mass_euler(n) = sum(u_euler(n, 1:end));
    mass_wendroff(n) = sum(u_wendroff(n, 1:end));
    
    L2_exact(n) = sqrt(sum(u_exact(n, 1:end).^2));
    L2_euler(n) = sqrt(sum(u_euler(n, 1:end).^2));
    L2_wendroff(n) = sqrt(sum(u_wendroff(n, 1:end).^2));
    
    u_exact_max(n) = max(u_exact(n, 1:end));
    u_euler_max(n) = max(u_euler(n, 1:end));
    u_wandroff_max(n) = max(u_wendroff(n, 1:end));
end

% If conservation holds, these should all be close to 0
display('--- Conservation of Mass ---');
str_1 = ['u_exact: ', num2str(max(mass_exact) - min(mass_exact))];
str_2 = ['u_euler: ', num2str(max(mass_euler) - min(mass_euler))];
str_3 = ['u_wendroff: ', num2str(max(mass_wendroff) - min(mass_wendroff))];
disp(str_1);
disp(str_2);
disp(str_3);

display('--- Conservation of L^2 norm ---');
str_1 = ['u_exact: ', num2str(max(L2_exact) - min(L2_exact))];
str_2 = ['u_euler: ', num2str(max(L2_euler) - min(L2_euler))];
str_3 = ['u_wendroff: ', num2str(max(L2_wendroff) - min(L2_wendroff))];
disp(str_1);
disp(str_2);
disp(str_3);

display('--- Conservation of maximums ---');
str_1 = ['u_exact: ', num2str(max(u_exact_max) - min(u_exact_max))];
str_2 = ['u_euler: ', num2str(max(u_euler_max) - min(u_euler_max))];
str_3 = sprintf('u_wendroff: %f', max(u_wandroff_max) - min(u_wandroff_max));
disp(str_1);
disp(str_2);
disp(str_3);

% Conservation of Mass
% - u_exact has conservation and then loses it
% - u_euler has fairly good convervation
% - u_wendroff also has fairly good conservation
% - results improve as N_x increases
plot(t, mass_exact, t, mass_euler, t, mass_wendroff);
title('Mass conservation');
legend('Exact', 'Euler', 'Wendroff');
figure;

% Conservation of L^2 
% - u_exact has conservation and then loses it exponentially fast (why?)
% - u_euler has fairly good convervation with periodic dips
% - u_wendroff also has fairly good conservation 
% - results improve as N_x increases
plot(t, L2_exact, t, L2_euler, t, L2_wendroff);
title('L^2 norm conservation');
legend('Exact', 'Euler', 'Wendroff');
figure;

% Conservation of max(u)
% - u_exact behaves similarly to the other converservations 
%   such that is starts conserved, but quickly loses conservation.
%   This should be investigated further
% - u_euler slowly loses convservation, but this improves with bigger N_x, N_t
% - u_wendroff is quite conserved, but at the bounds loses some conservation.
%   This is again improved as N_x and N_t are increased. 
plot(t, u_exact_max, t, u1_max, t, u_wandroff_max);
title('Maximum of u conservation');
legend('Exact', 'Euler', 'Wendroff');
