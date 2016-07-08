function nj_t = get_density(time, xk, vk, xj, Nx, Nk, Lx, S, M, wk, g_0)

nj_t = zeros(Nx, 1);
[Nj, xj, Xbin] = histcounts(xk(1, :), xj);
for j = 1:Nx
    ind1 = Xbin == j;
    if (j == Nx - 1)
        ind2 = Xbin == 1;
        ind3 = Xbin == Nx - 2;
        ind4 = Xbin == Nx - 3;       
        particles_near_j = [xk(time, ind1) (xk(time, ind2) + Lx) xk(time, ind3) xk(time, ind4)];
    elseif (j == 1 || j == Nx)
        ind1 = Xbin == 1;
        ind2 = Xbin == 2;
        ind3 = Xbin == Nx - 1;
        ind4 = Xbin == Nx - 2;        
        if (j == Nx)            
            particles_near_j = [(xk(time, ind1) + Lx) (xk(time, ind2) + Lx) xk(time, ind3) xk(time, ind4)];    
        else
            particles_near_j = [xk(time, ind1) xk(time, ind2) (xk(time, ind3) - Lx) (xk(time, ind4) - Lx)];    
        end
    elseif (j == 2)
        ind2 = Xbin == 3;
        ind3 = Xbin == 1;
        ind4 = Xbin == Nx - 1; 
        particles_near_j = [xk(time, ind1) xk(time, ind2) xk(time, ind3) (xk(time, ind4) - Lx)];
    else 
        ind2 = Xbin == j + 1;
        ind3 = Xbin == j - 1;
        ind4 = Xbin == j - 2;
        particles_near_j = [xk(time, ind1) xk(time, ind2) xk(time, ind3) xk(time, ind4)];
    end
    
    % TODO: Find a better way to do this with information already given
    k_pos = [find(ind1) find(ind2) find(ind3) find(ind4)];
    
    num_particles_near_j = size(particles_near_j);
    num_particles_near_j = num_particles_near_j(2);
    s = 0;    
    for particle_index = 1:num_particles_near_j        
        s = s + S(xj(j) - particles_near_j(particle_index)) * (wk(k_pos(particle_index)) ...
        - M(vk(time, k_pos(particle_index))) / g_0(vk(time, k_pos(particle_index))));    
    end
    nj_t(j) = 1 + (time / Nk) * s; 
end
