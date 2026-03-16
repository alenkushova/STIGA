function [u, grad_u, d1, d2, f] = Fourier_expansion(M)
syms U(x,t) GRADU(x,t) F(x,t) d1(x,t) d2(x,t)
U(x,t)= 0 ; GRADU(x,t)=0; F(x,t)=0;
GRADU(x,t)=  (cat (1, ...
             reshape (0*x , [1, size(x)]), ...
             reshape (0*t , [1, size(t)])));
D1(x,t) = 0*x;
D2(x,t) = 0*t;
for k = 1 : M 
    syms ek(x) pk(x) fk(t) uk(t) 
    omgk  = k*pi;
    norm  = int((sin(omgk*x))^2,[0 1]);
    ek(x) = sin(omgk*x)/norm;
    pk(x) = omgk*cos(omgk*x)/norm;
    fk(t) = exp(1i*omgk^2*t)/k;
    uk(t) = -1i*t*fk(t);
    U(x,t) = U(x,t) + uk(t).*ek(x);
    F(x,t) = F(x,t) + fk(t).*ek(x);
    GRADU(x,t) = GRADU(x,t) +(cat (1, ...
             reshape (uk(t).*pk(x), [1, size(x)]), ...
             reshape (exp(1i*omgk^2*t)*(-1i + t*omgk^2)/k.*ek(x) , [1, size(t)])));
    D1(x,t) = D1(x,t) + uk(t).*pk(x);
    D2(x,t) = D2(x,t) + exp(1i*omgk^2*t)*(-1i + t*omgk^2)/k.*ek(x);
end
u      = matlabFunction(U);
grad_u = matlabFunction(GRADU);
d1     = matlabFunction(D1); 
d2     = matlabFunction(D2); 
f      = matlabFunction(F);
end


