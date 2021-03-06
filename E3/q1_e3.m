L = 1;
T = 1.5;
a = 2; 
sigma = 0.1;
N_x = 500;
N_t = 10000;
h = L / N_x;

x = linspace(0, L, N_x);
t = linspace(0, T, N_t);

% CFL Stability
cfl = a * (T / N_t) / (L / N_x);
disp('--CFL--');
disp(cfl);
if (cfl <= 1)
    display('Parameters meet CFL condition: running program');
else
    display('Adjust parameters to meet CFL condition');
    display('Program ending');
    return;
end

u_0 = @(x) (1 / sqrt(2 * pi * sigma^2)) * exp(-1 * ...
    (x - (L / 2))^2 / (2 * sigma^2));

u_exact = zeros(N_t, N_x);
for n = 1:N_t
    for j = 1:N_x
        u_exact(n, j) = u_0(mod(x(j) - a * t(n), L));
    end
end

u_euler = q1_e3_euler_solver(L, T, a, sigma, N_x, N_t);
u_wendroff = q1_e3_wendroff_solver(L, T, a, sigma, N_x, N_t);
u_fourier = q1_e3_spectral_solver(L, T, a, sigma, N_x, N_t);

offset_1 = 10;
offset_2 = 10;
num_frames = 30;

x = x(offset_1:N_x - offset_2);
u_exact = u_exact(1:end, offset_1:N_x - offset_2);
u_euler = u_euler(1:end, offset_1:N_x - offset_2);
u_wendroff = u_wendroff(1:end, offset_1:N_x - offset_2);
u_fourier = u_fourier(1:end, offset_1:N_x - offset_2);

% make_movie(x, u_exact, u_euler, u_fourier, num_frames);

mass_exact = zeros(N_t, 1);
mass_euler = zeros(N_t, 1);
mass_wendroff = zeros(N_t, 1);
mass_fourier = zeros(N_t, 1);

L2_exact = zeros(N_t, 1);
L2_euler = zeros(N_t, 1);
L2_wendroff = zeros(N_t, 1);
L2_fourier = zeros(N_t, 1);

u_exact_max = zeros(N_t, 1);
u_euler_max = zeros(N_t, 1);
u_wandroff_max = zeros(N_t, 1);
u_fourier_max = zeros(N_t, 1);

for n = 1:N_t
    mass_exact(n) = sum(u_exact(n, 1:end));
    mass_euler(n) = sum(u_euler(n, 1:end));
    mass_wendroff(n) = sum(u_wendroff(n, 1:end));
    mass_fourier(n) = sum(u_fourier(n, 1:end));
    
    L2_exact(n) = sqrt(sum(u_exact(n, 1:end).^2));
    L2_euler(n) = sqrt(sum(u_euler(n, 1:end).^2));
    L2_wendroff(n) = sqrt(sum(u_wendroff(n, 1:end).^2));
    L2_fourier(n) = sqrt(sum(u_fourier(n, 1:end).^2));
    
    u_exact_max(n) = max(u_exact(n, 1:end));
    u_euler_max(n) = max(u_euler(n, 1:end));
    u_wandroff_max(n) = max(u_wendroff(n, 1:end));
    u_fourier_max(n) = max(u_fourier(n, 1:end));
end

% If conservation holds, these should all be close to 0
display('--- Conservation of Mass ---');
str_1 = ['u_exact: ', num2str(max(mass_exact) - min(mass_exact))];
str_2 = ['u_euler: ', num2str(max(mass_euler) - min(mass_euler))];
str_3 = ['u_wendroff: ', num2str(max(mass_wendroff) - min(mass_wendroff))];
str_4 = ['u_fourier: ', num2str(max(mass_fourier) - min(mass_fourier))];
disp(str_1);
disp(str_2);
disp(str_3);
disp(str_4);

display('--- Conservation of L^2 norm ---');
str_1 = ['u_exact: ', num2str(max(L2_exact) - min(L2_exact))];
str_2 = ['u_euler: ', num2str(max(L2_euler) - min(L2_euler))];
str_3 = ['u_wendroff: ', num2str(max(L2_wendroff) - min(L2_wendroff))];
str_4 = ['u_fourier: ', num2str(max(L2_fourier) - min(L2_fourier))];
disp(str_1);
disp(str_2);
disp(str_3);
disp(str_4);

display('--- Conservation of maximums ---');
str_1 = ['u_exact: ', num2str(max(u_exact_max) - min(u_exact_max))];
str_2 = ['u_euler: ', num2str(max(u_euler_max) - min(u_euler_max))];
str_3 = sprintf('u_wendroff: %f', max(u_wandroff_max) - min(u_wandroff_max));
str_4 = sprintf('u_fourier: %f', max(u_fourier_max) - min(u_fourier_max));
disp(str_1);
disp(str_2);
disp(str_3);
disp(str_4);

% Conservation of Mass
% - u_exact has conservation and then loses it
% - u_euler has fairly good convervation
% - u_wendroff also has fairly good conservation
% - results improve as N_x increases
plot(t, mass_exact, t, mass_euler, t, mass_wendroff, t, mass_fourier);
title('Mass conservation');
legend('Exact', 'Euler', 'Wendroff', 'Fourier');
figure;

% Conservation of L^2 
% - u_exact has conservation and then loses it exponentially fast (why?)
% - u_euler has fairly good convervation with periodic dips
% - u_wendroff also has fairly good conservation 
% - results improve as N_x increases
plot(t, L2_exact, t, L2_euler, t, L2_wendroff, t, L2_fourier);
title('L^2 norm conservation');
legend('Exact', 'Euler', 'Wendroff', 'Fourier');
figure;

% Conservation of max(u)
% - u_exact behaves similarly to the other converservations 
%   such that is starts conserved, but quickly loses conservation.
%   This should be investigated further
% - u_euler slowly loses convservation, but this improves with bigger N_x, N_t
% - u_wendroff is quite conserved, but at the bounds loses some conservation.
%   This is again improved as N_x and N_t are increased. 
plot(t, u_exact_max, t, u_euler_max, t, u_wandroff_max, t, u_fourier_max);
title('Maximum of u conservation');
legend('Exact', 'Euler', 'Wendroff', 'Fourier');

% Convergence tests
% L1_euler = max(h * sum(abs(u_exact - u_euler)));
% L2_euler = max(sqrt(h * sum(abs(u_exact - u_euler).^2)));
% Linf_euler = max(max(abs(u_exact - u_euler)));
% save(['euler_err' num2str(N_x)], 'h', 'L1_euler', 'L2_euler', 'Linf_euler'); 
% 
% L1_wendroff = max(h * sum(abs(u_exact - u_wendroff)));
% L2_wendroff = max(sqrt(h * sum(abs(u_exact - u_wendroff).^2)));
% Linf_wendroff= max(max(abs(u_exact - u_wendroff)));
% save(['wendroff_err' num2str(N_x)], 'h', 'L1_wendroff', 'L2_wendroff', 'Linf_wendroff'); 

% NVec = [50 100 200 500 1000 2000];
% hVec = zeros(6, 1);
% 
% L1Vec_euler = zeros(6, 1);
% L2Vec_euler = zeros(6, 1);
% LInfVec_euler = zeros(6, 1);
% 
% L1Vec_wendroff = zeros(6, 1);
% L2Vec_wendroff = zeros(6, 1);
% LInfVec_wendroff = zeros(6, 1);
% 
% for i = 1:6
%     err_a = load(['euler_err' num2str(NVec(i)) '.mat']);
%     err_b = load(['wendroff_err' num2str(NVec(i)) '.mat']);
%     
%     hVec(i) = err_a.h;
%     
%     L1Vec_euler(i) = err_a.L1_euler;
%     L2Vec_euler(i) = err_a.L2_euler;
%     LInfVec_euler(i) = err_a.Linf_euler;
%     
%     L1Vec_wendroff(i) = err_b.L1_wendroff;
%     L2Vec_wendroff(i) = err_b.L2_wendroff;
%     LInfVec_wendroff(i) = err_b.Linf_wendroff;
% end
% 
% figure
% loglog(hVec, L1Vec_euler, 'o-b', hVec, L2Vec_euler, 's-r', hVec, LInfVec_euler, '+-g', ...
%     hVec(1:3), 1e-1 * hVec(1:3), '-k');
% legend('L^1 error', 'L^2 error', 'L^\infty error', 'f(h) = c * h', ...
%     'location', 'northwest');
% xlabel('h');
% ylabel('Errors'); 
% title('Case Euler');
% 
% figure 
% loglog(hVec, L1Vec_wendroff, 'o-b', hVec, L2Vec_wendroff, 's-r', hVec, LInfVec_wendroff, '+-g', ...
%     hVec(1:3), hVec(1:3).^2, '-k');
% legend('L^1 error', 'L^2 error', 'L^\infty error', 'f(h) = c * h^2', ...
%     'location', 'northwest');
% xlabel('h');
% ylabel('Errors'); 
% title('Case Wendroff');
