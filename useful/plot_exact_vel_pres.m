function fig = plot_exact_vel_pres(problem_data, geo, npts, output_filename)
    arguments
        problem_data
        geo
        npts = 11
        output_filename = 'exact_solution'
    end

    % Time evaluation points
    t = {linspace(0,1,npts)};
    lables_t = geo.tgeo.map(t);

    % Spatial evaluation points (match original sampling)
    xv = linspace(0,1,30);
    xp = linspace(0,1,50);
    [Xv, Yv] = meshgrid(xv, xv);
    [Xp, Yp] = meshgrid(xp, xp);

    % Create video writer
    v = VideoWriter([output_filename '.mp4'], 'MPEG-4');
    v.FrameRate = 1;
    v.Quality = 100;
    open(v);

    % Evaluate exact solutions over time
    exact_velx = [];
    exact_presx = [];
    for i = 1:numel(t{:})
        tt = lables_t(i);

        vel = problem_data.velex(Xv, Yv, tt); % 2 × Nx × Ny
        pres = problem_data.presex(Xp, Yp, tt); % Nx × Ny

        exact_velx = cat(4, exact_velx, vel); % shape: [2, Nx, Ny, T]
        exact_presx = cat(3, exact_presx, pres); % shape: [Nx, Ny, T]
    end

    % Compute global color range for pressure
    min_pres = min(exact_presx(:));
    max_pres = max(exact_presx(:));

    % Plot and write frames
    fig = figure();
    fig.WindowState = 'maximized';

    for i = 1:numel(t{:})
        clf;

        % Velocity magnitude and vectors
        subplot(1,2,1)
        contourf(Xv, Yv, squeeze(sqrt(exact_velx(1,:,:,i).^2 + exact_velx(2,:,:,i).^2)));
        hold on;
        quiver(Xv, Yv, squeeze(exact_velx(1,:,:,i)), squeeze(exact_velx(2,:,:,i)), 'k');
        colorbar("southoutside"); view(2);
        title(sprintf('Exact Velocity at t=%.2f', lables_t(i)), 'FontSize', 10);
        axis tight; axis equal;
        hold off;

        % Pressure surface
        subplot(1,2,2)
        s = surf(Xp, Yp, exact_presx(:,:,i));
        s.FaceAlpha = 0.5; s.EdgeColor = 'none';
        colorbar("southoutside"); view(3);
        axis tight; axis square;
        title(sprintf('Exact Pressure at t=%.2f', lables_t(i)), 'FontSize', 10);
        % fig.CurrentAxes.CLim = [min_pres, max_pres];
        % fig.CurrentAxes.ZLim = [min_pres, max_pres];

        pause(1);
        frame = getframe(fig);
        writeVideo(v, frame);
    end

    close(v);
end
