function u = grad_heat_weak_sol(x,t,n)
u = 0;
for k = 1:2:n
u = u + 4*(exp(-(k*pi)^2*t)).*cos(k*pi*x);
end
end