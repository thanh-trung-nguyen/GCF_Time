function v = compute_v_from_laguerre_coef(qn,s,s0)
% compute v from the Laguerre's coefficient
% each column of qn is the function w.r.t. s


[Nq,Nx] = size(qn);

v = zeros(length(s),Nx);

for k = 0:Nq-1    
    vn = laguerre_integral_explicit(k,s,s0);
    
    v = v - vn'*qn(k+1,:);
end
