% HEAT_ST_SURF: <description>
%
%   CALL: 
%
%  INPUT:
%
% OUTPUT:
%
%
% ProjectName - STIGA
% Copyright (C) 2025 Alen Kushova
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% See <https://www.gnu.org/licenses/> for more details.

function fig = heat_st_surf(u, exu, space, geo, npts)
    arguments
        u
        exu 
        space
        geo
        npts = 11
    end
xp={linspace(0,1,50)}; % punti di valutazione in spazio
t={linspace(0,1,npts)}; % punti di valutazione in tempo
[Xp, Yp] = meshgrid(xp{:}, xp{:}); % punti sulla meshgrid 

lables_t = geo.tgeo.map(t);
pres = reshape(u,space.xsp_trial.ndof,space.tsp_trial.ndof);

eprest = [];
for sdof = 1:space.xsp_trial.ndof
    [epres, ~] = sp_eval (pres(sdof,:), space.tsp_trial, geo.tgeo, t);
    eprest = cat(1,eprest,epres');
end
epresx = [];
exact_presx = [];
for tpnt = 1: numel(t{:})
    tt = lables_t(tpnt);
    [epres, Fpres] = sp_eval (eprest(:,tpnt), space.xsp_trial, geo.xgeo, repmat(xp,[1,numel(space.xsp_trial.knots)]));
    epresex = exu(Xp, Yp, tt);
    exact_presx = cat(3, exact_presx, epresex); 
    epresx = cat(3,epresx,epres);
end



min_epres = min(epresx(:));
max_epres = max(epresx(:));
min_epresex = min(exact_presx(:));
max_epresex = max(exact_presx(:));
minimo = min(min_epres,min_epresex);
massimo= max(max_epres,max_epresex);

fig=figure();
fig.WindowState = 'maximized';

for it=1:numel(t{:})
    subplot(1,2,1)
    surf(squeeze(Fpres(1,:,:,1)),squeeze(Fpres(2,:,:,1)),epresx(:,:,it),'EdgeColor','none');
    colorbar("southoutside"); view(3);
    axis tight;
    title(sprintf('Solution at t=%.2f', lables_t(it)), 'FontSize', 10);
    fig.CurrentAxes.CLim=[minimo,massimo];
    fig.CurrentAxes.ZLim=[minimo,massimo];

    subplot(1,2,2)
    surf(Xp,Yp,exact_presx(:,:,it),'EdgeColor','none');
    colorbar("southoutside"); view(3);
    axis tight;
    title(sprintf('Exact solution at t=%.2f', lables_t(it)), 'FontSize', 10);
    fig.CurrentAxes.CLim=[minimo,massimo];
    fig.CurrentAxes.ZLim=[minimo,massimo];

    pause(2);
end

end
