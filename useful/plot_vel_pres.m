function fig=plot_vel_pres(vel, pres, space, geo, npts, output_filename)
    arguments
        vel
        pres
        space
        geo
        npts = 11
        output_filename = 'solution'
    end
xv={linspace(0,1,30)}; % punti di valutazione in spazio
xp={linspace(0,1,50)}; % punti di valutazione in spazio
t={linspace(0,1,npts)}; % punti di valutazione in tempo
lables_t = geo.tgeo.map(t);
vel=reshape(vel,space.spv.ndof,space.spt_vel.ndof);
pres=reshape(pres,space.spp.ndof,space.spt_pres.ndof);

evelt = [];
eprest = [];
for sdof = 1:space.spv.ndof
    [evel, ~] = sp_eval (vel(sdof,:), space.spt_vel, geo.tgeo, t);
    evelt = cat(1,evelt,evel');
end
for sdof = 1:space.spp.ndof
    [epres, ~] = sp_eval (pres(sdof,:), space.spt_pres, geo.tgeo, t);
    eprest = cat(1,eprest,epres');
end

evelx = [];
epresx = [];
for tpnt = 1: numel(t{:})
    [evel, Fvel] = sp_eval (evelt(:,tpnt), space.spv, geo.xgeo, repmat(xv, [1, space.spv.ncomp]));
    [epres, Fpres] = sp_eval (eprest(:,tpnt), space.spp, geo.xgeo, repmat(xp,[1,numel(space.spp.knots)]));
    evelx = cat(4,evelx,evel);
    epresx = cat(3,epresx,epres);
end

min_epres = min(epresx(:));
max_epres = max(epresx(:));

% Create video writer object
v = VideoWriter([output_filename '.mp4'], 'MPEG-4'); % Use .mp4 if supported on your system
v.FrameRate = 1;                      % Set frames per second
v.Quality = 100;                      % Set highest possible quality
open(v);

fig=figure();
fig.WindowState = 'maximized';

for it=1:numel(t{:})
    subplot(1,2,1)
    contourf(squeeze(Fvel(1,:,:)), squeeze(Fvel(2,:,:)), ... 
             sqrt(squeeze(evelx(1,:,:,it)).^2 + squeeze(evelx(2,:,:,it)).^2));
    % s = surf(squeeze(Fvel(1,:,:)),squeeze(Fvel(2,:,:)),sqrt(squeeze(evelx(1,:,:,it)).^2 + squeeze(evelx(2,:,:,it)).^2));
    hold on;
    % s.FaceAlpha =0.5; s.EdgeColor = "none";
    colorbar("southoutside"); view(2);    
    quiver(squeeze(Fvel(1,:,:)),squeeze(Fvel(2,:,:)),squeeze(evelx(1,:,:,it)),squeeze(evelx(2,:,:,it)),'Color','black'); 
    title(sprintf('Velocity at t=%.2f', lables_t(it)), 'FontSize', 10);
    axis tight; axis equal;
    hold off;

    subplot(1,2,2)
    % contourf(squeeze(Fpres(1,:,:,1)),squeeze(Fpres(2,:,:,1)),epresx(:,:,it));
    s = surf(squeeze(Fpres(1,:,:,1)),squeeze(Fpres(2,:,:,1)),epresx(:,:,it));
    s.FaceAlpha =0.5; s.EdgeColor = "none";
    colorbar("southoutside"); view(3);
    axis tight; axis square;
    title(sprintf('Pressure at t=%.2f', lables_t(it)), 'FontSize', 10);
    % fig.CurrentAxes.CLim=[min_epres,max_epres];
    % fig.CurrentAxes.ZLim=[min_epres,max_epres];

    pause(2);
    % Capture frame
    frame = getframe(fig);
    writeVideo(v, frame);    
end
close(v)
end