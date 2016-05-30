function qn = laguerre_coefficient(Idx,q,s,s0)
% % compute qn(x) = int_{s0}^infty q(s) f_n(s) ds, where f_n(s) = L_n(s-s0)
% q can be a COLUMN vector or a matrix whose columns are functions of s.
% Idx: a vector of indices of the Laguerre's expansion

if size(s,2)==1
    s = s';
end

smax = 1000; 
ds2 = 0.0001; 
ds = s(2) - s(1);

Nx = size(q,2);

s2 = s(end):ds2:smax;

Nq = length(Idx);
qn = zeros(Nq,Nx);

for idx = 1:Nq
    n = Idx(idx);
    fn = laguerre(n,s,s0);
    fn2 = laguerre(n,s2,s0);

    qn(idx,:) = fn*q(:,:)*ds + q(end,:)*s(end)^2*sum(fn2./(s2.^2))*ds2;
end

