function problemData = computeProblemData2d(sym_velexx, sym_velexy, sym_presex, viscosity)

syms x y t;

velexx=matlabFunction(sym_velexx,'Vars',[x,y,t]);
problemData.velexx=@(x,y,t) velex(x,y,t)+0*x+0*y+0*t;
velexy=matlabFunction(sym_velexy,'Vars',[x,y,t]);
problemData.velexy=@(x,y,t) velex(x,y,t)+0*x+0*y+0*t;
problemData.velex=@(x,y,t) cat (1, ...
    reshape (velexx(x,y,t), [1, size(x)]), ...
    reshape (velexy(x,y,t), [1, size(x)]));


presex=matlabFunction(sym_presex,'Vars',[x,y,t]);
problemData.presex=@(x,y,t) presex(x,y,t)+0*x+0*y+0*t;

% gradiente velocità componente x
sym_dx_velexx=diff(sym_velexx,x);
dx_velexx=matlabFunction(sym_dx_velexx,'Vars',[x,y,t]);
dx_velexx=@(x,y,t) dx_velexx(x,y,t)+0*x+0*y+0*t;
sym_dy_velexx=diff(sym_velexx,y);
dy_velexx=matlabFunction(sym_dy_velexx,'Vars',[x,y,t]);
dy_velexx=@(x,y,t) dy_velexx(x,y,t)+0*x+0*y+0*t;
sym_dt_velexx=diff(sym_velexx,t);
dt_velexx=matlabFunction(sym_dt_velexx,'Vars',[x,y,t]);
dt_velexx=@(x,y,t) dt_velexx(x,y,t)+0*x+0*y+0*t;
% problemData.dt_velexx=dt_velexx;

problemData.grad_velexx=@(x,y,t) cat (1, ...
    reshape (dx_velexx(x,y,t), [1, size(x)]), ...
    reshape (dy_velexx(x,y,t), [1, size(x)]));

% gradiente velocità componente y
sym_dx_velexy=diff(sym_velexy,x);
dx_velexy=matlabFunction(sym_dx_velexy,'Vars',[x,y,t]);
dx_velexy=@(x,y,t) dx_velexy(x,y,t)+0*x+0*y+0*t;
sym_dy_velexy=diff(sym_velexy,y);
dy_velexy=matlabFunction(sym_dy_velexy,'Vars',[x,y,t]);
dy_velexy=@(x,y,t) dy_velexy(x,y,t)+0*x+0*y+0*t;
sym_dt_velexy=diff(sym_velexy,t);
dt_velexy=matlabFunction(sym_dt_velexy,'Vars',[x,y,t]);
dt_velexy=@(x,y,t) dt_velexy(x,y,t)+0*x+0*y+0*t;
% problemData.dt_velexy=dt_velexy;

problemData.grad_velexy=@(x,y,t) cat (1, ...
    reshape (dx_velexy(x,y,t), [1, size(x)]), ...
    reshape (dy_velexy(x,y,t), [1, size(x)]));

problemData.grad_velex = @(x,y,t) cat(1, ...
    reshape (problemData.grad_velexx(x,y,t) ,[1, 2, size(x)]), ...
    reshape (problemData.grad_velexy(x,y,t) ,[1, 2, size(x)]));

problemData.dt_velex = @(x,y,t) cat(1, ...
    reshape (dt_velexx(x,y,t) ,[1, size(x)]), ...
    reshape (dt_velexy(x,y,t), [1, size(x)]));

% gradiente presione
sym_dx_presex=diff(sym_presex,x);
dx_presex=matlabFunction(sym_dx_presex,'Vars',[x,y,t]);
dx_presex=@(x,y,t) dx_presex(x,y,t)+0*x+0*y+0*t;
sym_dy_presex=diff(sym_presex,y);
dy_presex=matlabFunction(sym_dy_presex,'Vars',[x,y,t]);
dy_presex=@(x,y,t) dy_presex(x,y,t)+0*x+0*y+0*t;
sym_dt_presex=diff(sym_presex,t);
dt_presex=matlabFunction(sym_dt_presex,'Vars',[x,y,t]);
dt_presex=@(x,y,t) dt_presex(x,y,t)+0*x+0*y+0*t;
% problemData.dt_presex=dt_presex;

problemData.grad_presex=@(x,y,t) cat (1, ...
    reshape (dx_presex(x,y,t), [1, size(x)]), ...
    reshape (dy_presex(x,y,t), [1, size(x)]),...
    reshape (dt_presex(x,y,t), [1, size(x)]));

sym_dxx_velexx=diff(sym_dx_velexx,x);
dxx_velexx=matlabFunction(sym_dxx_velexx,'Vars',[x,y,t]);
dxx_velexx=@(x,y,t) dxx_velexx(x,y,t)+0*x+0*y+0*t;
sym_dyy_velexx=diff(sym_dy_velexx,y);
dyy_velexx=matlabFunction(sym_dyy_velexx,'Vars',[x,y,t]);
dyy_velexx=@(x,y,t) dyy_velexx(x,y,t)+0*x+0*y+0*t;


sym_dxx_velexy=diff(sym_dx_velexy,x);
dxx_velexy=matlabFunction(sym_dxx_velexy,'Vars',[x,y,t]);
dxx_velexy=@(x,y,t) dxx_velexy(x,y,t)+0*x+0*y+0*t;
sym_dyy_velexy=diff(sym_dy_velexy,y);
dyy_velexy=matlabFunction(sym_dyy_velexy,'Vars',[x,y,t]);
dyy_velexy=@(x,y,t) dyy_velexy(x,y,t)+0*x+0*y+0*t;

problemData.viscosity=viscosity;

f1 =@(x,y,t) dt_velexx(x,y,t) -viscosity*(dxx_velexx(x,y,t) + dyy_velexx(x,y,t)) + dx_presex(x,y,t);
f2 =@(x,y,t) dt_velexy(x,y,t) -viscosity*(dxx_velexy(x,y,t) + dyy_velexy(x,y,t)) + dy_presex(x,y,t);

problemData.f=@(x,y,t) cat (1, ...
    reshape (f1(x,y,t), [1, size(x)]), ...
    reshape (f2(x,y,t), [1, size(x)]));
end