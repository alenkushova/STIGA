% HEAT_ST_PLOT: <description>
%
%   CALL: 
%
%  INPUT:
%
% OUTPUT:
%
%
% ProjectName - STIGA
% Copyright (C) 2025 Alen Kushova, Gabriele Loli
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

function fig = heat_st_plot(u, exu, space, geo, npts)
    arguments
        u
        exu 
        space
        geo
        npts = 11
    end
xp={linspace(0,1,50)}; % punti di valutazione in spazio
t={linspace(0,1,npts)}; % punti di valutazione in tempo
lables_t = geo.tgeo.map(t);
u=reshape(u,space.xsp_trial.ndof,space.tsp_trial.ndof);

esolt = [];
for sdof = 1:space.xsp_trial.ndof
    [evs, ~] = sp_eval (u(sdof,:), space.tsp_trial, geo.tgeo, t);
    esolt = cat(1,esolt,evs');
end

esolx = [];
euex  = [];
for tpnt = 1: numel(t{:})
    tt = lables_t(tpnt);
    euex = cat(3,euex, exu(xp{:}, tt));
    [evs, F] = sp_eval (esolt(:,tpnt), space.xsp_trial, geo.xgeo, repmat(xp,[1,numel(space.xsp_trial.knots)]));
    esolx = cat(3,esolx,evs);
end

min_eu = min(esolx(:));
max_eu = max(esolx(:));
min_eexu = min(euex(:));
max_eexu = max(euex(:));
min_esol = min(min_eu,min_eexu);
max_esol = max(max_eu,max_eexu);

fig=figure();
% fig.WindowState = 'maximized';

for it=1:numel(t{:})
plot(squeeze(F(1,:,1)),esolx(:,:,it),'LineWidth',1.5);
hold on 
plot(squeeze(F(1,:,1)),euex(:,:,it),'-.','LineWidth',1.5);
%axis tight; axis square;
title('Heat solution', 'FontSize', 10);
hold off 
ylim([min_esol,max_esol]);
pause(2)

end

end