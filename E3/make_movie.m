function make_movie(x, u_exact, u1, u2, num_frames)

axes tight manual

frames(num_frames) = struct('cdata',[],'colormap',[]);
t = 1;

for i = 1:num_frames
    plot(x, u_exact(t, 1:end), ...
         x, u1(t, 1:end), ...
         x, u2(t, 1:end));
    legend('Exact', 'Euler', 'Lax Wendroff');
    drawnow
    frames(i) = getframe;
    t = t + 200;
end

movie(frames, 5, 2);
end