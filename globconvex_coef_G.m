function G = globconvex_coef_G(Nq,s0)
% compute the coefficient F in the globally convergent method using
% Laguerre functions. 
% Nq: number of Laguerre's coefficients used in the expansion. 
% s0: the lower bound of the pseudo frequency


smax = 200; 
ds = 0.00005; 
Ns = round(smax-s0)/ds + 1;
s = linspace(s0,smax,Ns);



% compute fn: 
f = zeros(Nq,Ns);
for n = 1:Nq
    f(n,:) = laguerre_explicit(n-1,s,s0);
end

% compute vn: 
v = zeros(Nq,Ns);
for n = 1:Nq
    v(n,:) = laguerre_integral_explicit(n-1,s,s0);
end

% compute F: 
G = zeros(Nq,Nq);
for k = 1:Nq
    for n = 1:Nq
        G(k,n) = sum(2*s.*f(k,:).*f(n,:) - 2*f(k,:).*v(n,:))*ds; 
    end
end
    




