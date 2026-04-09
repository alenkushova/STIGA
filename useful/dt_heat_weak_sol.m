function u = dt_heat_weak_sol(x,t,n)
u = 0;
for k = 1:2:n
u = u - 4*k*pi*(exp(-(k*pi)^2*t)).*sin(k*pi*x);
end
end